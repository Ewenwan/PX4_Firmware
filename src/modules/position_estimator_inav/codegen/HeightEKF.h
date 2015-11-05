/*
 * File: HeightEKF.h
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 15-Oct-2015 20:51:56
 */

#ifndef __HEIGHTEKF_H__
#define __HEIGHTEKF_H__

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "HeightEKF_types.h"

/* Function Declarations */
extern void HeightEKF(float x_apo[2], float P_apo[4], const unsigned char zFlag
                      [2], float dt, const float z[2], float r_baro, float r_acc,
                      float xa_apo[2], float Pa_apo[4]);

#endif

/*
 * File trailer for HeightEKF.h
 *
 * [EOF]
 */
