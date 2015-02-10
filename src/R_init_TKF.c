#include "TKF91.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {
  /* TKF91.c */
  CALLMETHOD_DEF(TKF91LikelihoodFunctionWrapper, 7),
  CALLMETHOD_DEF(TKF91LikelihoodFunction1DMain, 6),
  CALLMETHOD_DEF(TKF91LikelihoodFunction2DMainNM, 5),

  /* TKF92.c */

  {NULL, NULL, 0}
};

void R_init_TKF(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  return;
}

