/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "superlu_bridge.h"
#include "superlu_hooks.h"
#include <slu_ddefs.h>
#include <setjmp.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

_Static_assert(sizeof(int_t)==sizeof(int),"bridge requires locked 32-bit SuperLU indices");
_Static_assert(SUPERLU_MAJOR_VERSION==7 && SUPERLU_MINOR_VERSION==0 && SUPERLU_PATCH_VERSION==0,
               "review bridge allocation guards when changing SuperLU");

/* Every allocation, including temporary arrays and headers, belongs to a
   factor arena. Cleanup therefore also covers upstream early-return paths. */
typedef union cs_block cs_block;
union cs_block {
    max_align_t alignment;
    struct {cs_block *previous,*next;struct cs_lu *owner;size_t bytes;} h;
};
struct cs_lu {
    cs_block *head;
    size_t used,peak,budget;
    jmp_buf failure;
    int code,n,ready;
    char error[256];
    SuperMatrix l,u;
    int *perm_c,*perm_r;
    double rcond;
};
#if defined(_MSC_VER)
static __declspec(thread) cs_lu *active;
#else
static _Thread_local cs_lu *active;
#endif
static void fail(int code,const char *message) {
    /* An unguarded call is an internal programming error, never input-driven. */
    if(!active) abort();
    active->code=code;
    snprintf(active->error,sizeof(active->error),"%s",message);
    longjmp(active->failure,1);
}
void *cs_superlu_malloc(size_t bytes) {
    cs_block *block;
    if(!active) abort();
    if(bytes>SIZE_MAX-sizeof(cs_block) || bytes+sizeof(cs_block)>active->budget-active->used)
        fail(CS_LU_MEMORY,"SuperLU allocation exceeds factor memory budget");
    block=(cs_block *)malloc(sizeof(cs_block)+bytes);
    if(!block) fail(CS_LU_MEMORY,"SuperLU allocation failed");
    block->h.previous=NULL;block->h.next=active->head;block->h.owner=active;block->h.bytes=sizeof(cs_block)+bytes;
    if(active->head)active->head->h.previous=block;
    active->head=block;active->used+=block->h.bytes;
    if(active->used>active->peak)active->peak=active->used;
    return block+1;
}
void cs_superlu_free(void *address) {
    cs_block *block;cs_lu *owner;
    if(!address)return;
    block=(cs_block *)address-1;owner=block->h.owner;
    if(block->h.previous)block->h.previous->h.next=block->h.next;
    else owner->head=block->h.next;
    if(block->h.next)block->h.next->h.previous=block->h.previous;
    owner->used-=block->h.bytes;free(block);
}
void cs_superlu_abort(const char *message) {fail(CS_LU_NUMERICAL,message);}
static void release_until(cs_lu *factor,cs_block *stop) {
    while(factor->head!=stop)cs_superlu_free(factor->head+1);
}
cs_lu *cs_lu_create(size_t budget) {
    cs_lu *factor;
    if(budget<sizeof(cs_lu))return NULL;
    factor=(cs_lu *)calloc(1,sizeof(cs_lu));
    if(factor) {factor->used=factor->peak=sizeof(cs_lu);factor->budget=budget<(size_t)INT_MAX?budget:(size_t)INT_MAX;}
    return factor;
}
int cs_lu_factor(cs_lu *factor,int n,int nnz,int *colptr,int *rows,double *values,double norm) {
    cs_lu *previous=active;
    SuperLUStat_t stat;
    superlu_options_t options;
    NCformat storage;
    DNformat rhs_storage;
    SuperMatrix a,b;
    int info=0;
    /* Upstream's workspace and initial fill estimates use signed int
       arithmetic even with size_t allocators. Guard those products first. */
    const uint64_t panel=(uint64_t)sp_ienv(1),rowblk=(uint64_t)sp_ienv(4);
    const uint64_t maxsuper=(uint64_t)SUPERLU_MAX(sp_ienv(3),sp_ienv(7));
    const uint64_t extent=n>0?(uint64_t)n:0;
    const uint64_t iwork=(2*panel+2+NO_MARKER)*extent*sizeof(int);
    const uint64_t temp=SUPERLU_MAX(extent,(maxsuper+rowblk)*panel);
    const uint64_t dwork=(extent*panel+temp)*sizeof(double);
    if(n<=0 || nnz<0 || iwork>INT_MAX || dwork>INT_MAX ||
       (uint64_t)nnz*(uint64_t)sp_ienv(6)>INT_MAX) {
        snprintf(factor->error,sizeof(factor->error),"SuperLU signed workspace extent exceeded");
        return CS_LU_ARGUMENT;
    }
    active=factor;
    if(setjmp(factor->failure)) {
        release_until(factor,NULL);active=previous;return factor->code;
    }
    factor->n=n;factor->perm_c=int32Malloc(n);factor->perm_r=int32Malloc(n);
    storage.nnz=nnz;storage.nzval=values;storage.rowind=rows;storage.colptr=colptr;
    a.Stype=SLU_NC;a.Dtype=SLU_D;a.Mtype=SLU_GE;a.nrow=n;a.ncol=n;a.Store=&storage;
    /* The zero-column RHS still gets a real address: dgstrs forms &Bmat[0] outside its
       per-column loops, and with NULL that is pointer arithmetic on a null pointer
       (UBSan, T25). Nothing is read or written through it (no columns; the locked
       build uses SuperLU's internal CBLAS, so no vendor path offsets it). */
    double no_columns=0.0;
    rhs_storage.lda=n;rhs_storage.nzval=&no_columns;
    b.Stype=SLU_DN;b.Dtype=SLU_D;b.Mtype=SLU_GE;b.nrow=n;b.ncol=0;b.Store=&rhs_storage;
    set_default_options(&options);
    options.ColPerm=COLAMD;options.DiagPivotThresh=1.0;options.SymmetricMode=NO;
    options.Equil=NO;options.IterRefine=NOREFINE;options.PrintStat=NO;
    StatInit(&stat);
    /* dgssv uses gstrf without equilibration/refinement, like scipy splu.
       A zero-column RHS builds factors once; later solves call dgstrs. */
    dgssv(&options,&a,factor->perm_c,factor->perm_r,&factor->l,&factor->u,&b,&stat,&info);
    if(info) {
        char message[100];snprintf(message,sizeof(message),"SuperLU factorization failed: info=%d (singular pivot if <=n)",info);
        fail(CS_LU_NUMERICAL,message);
    }
    {char norm_kind='1';dgscon(&norm_kind,&factor->l,&factor->u,norm,&factor->rcond,&stat,&info);}
    if(info)fail(CS_LU_NUMERICAL,"SuperLU condition estimate failed");
    StatFree(&stat);factor->ready=1;active=previous;return CS_LU_OK;
}
int cs_lu_solve(cs_lu *factor,int columns,double *rhs,size_t budget) {
    cs_lu *previous=active;cs_block *checkpoint=factor->head;
    SuperLUStat_t stat;DNformat storage;SuperMatrix b;int info=0;
    if(!factor->ready || columns<=0 || columns>INT_MAX/factor->n) {
        snprintf(factor->error,sizeof(factor->error),"invalid SuperLU solve extent");return CS_LU_ARGUMENT;
    }
    if(budget<factor->used) {
        snprintf(factor->error,sizeof(factor->error),"SuperLU solve exceeds factor memory budget");return CS_LU_MEMORY;
    }
    factor->budget=budget<(size_t)INT_MAX?budget:(size_t)INT_MAX;factor->peak=factor->used;active=factor;
    if(setjmp(factor->failure)) {
        release_until(factor,checkpoint);active=previous;return factor->code;
    }
    storage.lda=factor->n;storage.nzval=rhs;
    b.Stype=SLU_DN;b.Dtype=SLU_D;b.Mtype=SLU_GE;b.nrow=factor->n;b.ncol=columns;b.Store=&storage;
    StatInit(&stat);
    dgstrs(NOTRANS,&factor->l,&factor->u,factor->perm_c,factor->perm_r,&b,&stat,&info);
    if(info)fail(CS_LU_NUMERICAL,"SuperLU triangular solve failed");
    StatFree(&stat);active=previous;return CS_LU_OK;
}
size_t cs_lu_bytes(const cs_lu *factor) {return factor->used;}
size_t cs_lu_peak(const cs_lu *factor) {return factor->peak;}
double cs_lu_rcond(const cs_lu *factor) {return factor->rcond;}
const char *cs_lu_error(const cs_lu *factor) {return factor->error;}
void cs_lu_destroy(cs_lu *factor) {if(factor) {release_until(factor,NULL);free(factor);}}
