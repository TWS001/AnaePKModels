#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void PropEleveldFinal(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void PropSchnider(void *, void *, void *, void *, void *);
extern void PropSchnider2(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"PropEleveldFinal", (DL_FUNC) &PropEleveldFinal, 9},
    {"PropSchnider",     (DL_FUNC) &PropSchnider,     5},
    {"PropSchnider2",    (DL_FUNC) &PropSchnider2,    5},
    {NULL, NULL, 0}
};

void R_init_AnaePKModels(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
