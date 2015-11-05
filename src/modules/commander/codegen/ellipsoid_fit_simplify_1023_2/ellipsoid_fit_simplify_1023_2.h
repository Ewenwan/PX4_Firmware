/*
 * File: ellipsoid_fit_simplify_1023_2.h
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 24-Oct-2015 14:16:57
 */

#ifndef __ELLIPSOID_FIT_SIMPLIFY_1023_2_H__
#define __ELLIPSOID_FIT_SIMPLIFY_1023_2_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "ellipsoid_fit_simplify_1023_2_types.h"

/* Function Declarations */
extern void ellipsoid_fit_simplify_1023_2(const emxArray_real32_T *x, const
  emxArray_real32_T *y, const emxArray_real32_T *z, float N, float center[3],
  float radii[3]);
extern void ellipsoid_fit_simplify_1023_2_initialize(void);
extern void ellipsoid_fit_simplify_1023_2_terminate(void);
extern emxArray_real32_T *emxCreateND_real32_T(int numDimensions, int *size);
extern emxArray_real32_T *emxCreateWrapperND_real32_T(float *data, int
  numDimensions, int *size);
extern emxArray_real32_T *emxCreateWrapper_real32_T(float *data, int rows, int
  cols);
extern emxArray_real32_T *emxCreate_real32_T(int rows, int cols);
extern void emxDestroyArray_real32_T(emxArray_real32_T *emxArray);
extern void emxInitArray_real32_T(emxArray_real32_T **pEmxArray, int
  numDimensions);

#endif

/*
 * File trailer for ellipsoid_fit_simplify_1023_2.h
 *
 * [EOF]
 */
