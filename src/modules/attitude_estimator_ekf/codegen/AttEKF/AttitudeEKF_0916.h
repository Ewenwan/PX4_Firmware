/*
 * File: AttitudeEKF_0916.h
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 21-Oct-2015 14:29:05
 */

#ifndef __ATTITUDEEKF_0916_H__
#define __ATTITUDEEKF_0916_H__

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "AttitudeEKF_0916_types.h"

/* Function Declarations */
extern void AttitudeEKF_0916(unsigned char approx_order, unsigned char use_J,
  const float X_apo_k[12], const float P_apo_k[144], const unsigned char zFlag[3],
  float dt, const float z[9], const float q_vector[4], const float r_vector[3],
  const float J[9], float xa_apo[12], float Pa_apo[144], float Rot_matrix[9],
  float eulerAngles[3]);
extern void AttitudeEKF_0916_initialize(void);
extern void AttitudeEKF_0916_terminate(void);

#endif

/*
 * File trailer for AttitudeEKF_0916.h
 *
 * [EOF]
 */
