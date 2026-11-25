/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CLEANSEAMS_SUPERLU_HOOKS_H
#define CLEANSEAMS_SUPERLU_HOOKS_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
void *cs_superlu_malloc(size_t);
void cs_superlu_free(void *);
void cs_superlu_abort(const char *);
#ifdef __cplusplus
}
#endif
#define USER_MALLOC cs_superlu_malloc
#define USER_FREE cs_superlu_free
#define USER_ABORT cs_superlu_abort
#endif
