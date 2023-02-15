#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
  Code extracted from https://cran.r-project.org/src/contrib/Archive/condSURV/ - 2022-04-29
*/

/* .C calls */
extern void SurvBeranKernel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

// static const R_CMethodDef CEntries[] = {
//   {"SurvBeranKernel", (DL_FUNC) &SurvBeranKernel, 10},
//   {NULL, NULL, 0}
// };
