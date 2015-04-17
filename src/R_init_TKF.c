#include "TKF91.h"
#include "TKF92.h"
#include "TKF92HG.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {
  /* TKF91.c */
  CALLMETHOD_DEF(TKF91LikelihoodFunctionWrapper, 7),
  CALLMETHOD_DEF(TKF91LikelihoodFunction1DMain, 6),
  CALLMETHOD_DEF(TKF91LikelihoodFunction2DMainNM, 5),

  /* TKF92.c */
  CALLMETHOD_DEF(TKF92LikelihoodFunctionWrapper, 8),
  CALLMETHOD_DEF(TKF92LikelihoodFunction1DMain, 7),
  CALLMETHOD_DEF(TKF92LikelihoodFunction3DMainNM, 5),

  /* TKF92HG.c */
  CALLMETHOD_DEF(TKF92HGLikelihoodFunctionWrapper, 10),
  CALLMETHOD_DEF(TKF92HGLikelihoodFunction1DMain, 9),
  CALLMETHOD_DEF(TKF92HGLikelihoodFunction5DMainNM, 5),
  {NULL, NULL, 0}
};

void R_init_TKF(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  return;
}

