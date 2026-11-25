/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CLEANSEAMS_SUPERLU_BRIDGE_H
#define CLEANSEAMS_SUPERLU_BRIDGE_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct cs_lu cs_lu;
/* All setjmp/longjmp frames are C frames. No C++ object is bypassed. */
enum { CS_LU_OK=0, CS_LU_MEMORY=1, CS_LU_NUMERICAL=2, CS_LU_ARGUMENT=3 };
cs_lu *cs_lu_create(size_t budget);
int cs_lu_factor(cs_lu *,int n,int nnz,int *colptr,int *rows,double *values,double norm);
int cs_lu_solve(cs_lu *,int columns,double *column_major_rhs,size_t budget);
size_t cs_lu_bytes(const cs_lu *);
size_t cs_lu_peak(const cs_lu *);
double cs_lu_rcond(const cs_lu *);
const char *cs_lu_error(const cs_lu *);
void cs_lu_destroy(cs_lu *);
#ifdef __cplusplus
}
#endif
#endif
