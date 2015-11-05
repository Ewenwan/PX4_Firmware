/*
 * File: AttitudeEKF_0916.c
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 21-Oct-2015 14:29:05
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "AttitudeEKF_0916.h"

/* Function Declarations */
static void b_mrdivide(float A[108], const float B[81]);
static void c_mrdivide(float A[72], const float B[36]);
static void cross(const float a[3], const float b[3], float c[3]);
static void diag(const float v[12], float d[144]);
static void eye(double I[144]);
static void inv(const float x[9], float y[9]);
static void mpower(const float a[9], float c[9]);
static void mrdivide(const float A[36], const float B[9], float y[36]);
static float norm(const float x[3]);
static void rdivide(const float x[3], float y, float z[3]);
static float rt_atan2f_snf(float u0, float u1);

/* Function Definitions */

/*
 * Arguments    : float A[108]
 *                const float B[81]
 * Return Type  : void
 */
static void b_mrdivide(float A[108], const float B[81])
{
  static float b_A[81];
  signed char ipiv[9];
  int k;
  int j;
  int c;
  int kBcol;
  int ix;
  float temp;
  float s;
  int i;
  int jp;
  int jAcol;
  memcpy(&b_A[0], &B[0], 81U * sizeof(float));
  for (k = 0; k < 9; k++) {
    ipiv[k] = (signed char)(1 + k);
  }

  for (j = 0; j < 8; j++) {
    c = j * 10;
    kBcol = 0;
    ix = c;
    temp = (real32_T)fabs(b_A[c]);
    for (k = 2; k <= 9 - j; k++) {
      ix++;
      s = (real32_T)fabs(b_A[ix]);
      if (s > temp) {
        kBcol = k - 1;
        temp = s;
      }
    }

    if (b_A[c + kBcol] != 0.0F) {
      if (kBcol != 0) {
        ipiv[j] = (signed char)((j + kBcol) + 1);
        ix = j;
        kBcol += j;
        for (k = 0; k < 9; k++) {
          temp = b_A[ix];
          b_A[ix] = b_A[kBcol];
          b_A[kBcol] = temp;
          ix += 9;
          kBcol += 9;
        }
      }

      k = (c - j) + 9;
      for (i = c + 1; i + 1 <= k; i++) {
        b_A[i] /= b_A[c];
      }
    }

    jp = c;
    jAcol = c + 9;
    for (kBcol = 1; kBcol <= 8 - j; kBcol++) {
      temp = b_A[jAcol];
      if (b_A[jAcol] != 0.0F) {
        ix = c + 1;
        k = (jp - j) + 18;
        for (i = 10 + jp; i + 1 <= k; i++) {
          b_A[i] += b_A[ix] * -temp;
          ix++;
        }
      }

      jAcol += 9;
      jp += 9;
    }
  }

  for (j = 0; j < 9; j++) {
    jp = 12 * j;
    jAcol = 9 * j;
    for (k = 1; k <= j; k++) {
      kBcol = 12 * (k - 1);
      if (b_A[(k + jAcol) - 1] != 0.0F) {
        for (i = 0; i < 12; i++) {
          A[i + jp] -= b_A[(k + jAcol) - 1] * A[i + kBcol];
        }
      }
    }

    temp = 1.0F / b_A[j + jAcol];
    for (i = 0; i < 12; i++) {
      A[i + jp] *= temp;
    }
  }

  for (j = 8; j > -1; j += -1) {
    jp = 12 * j;
    jAcol = 9 * j - 1;
    for (k = j + 2; k < 10; k++) {
      kBcol = 12 * (k - 1);
      if (b_A[k + jAcol] != 0.0F) {
        for (i = 0; i < 12; i++) {
          A[i + jp] -= b_A[k + jAcol] * A[i + kBcol];
        }
      }
    }
  }

  for (kBcol = 7; kBcol > -1; kBcol += -1) {
    if (ipiv[kBcol] != kBcol + 1) {
      jp = ipiv[kBcol] - 1;
      for (jAcol = 0; jAcol < 12; jAcol++) {
        temp = A[jAcol + 12 * kBcol];
        A[jAcol + 12 * kBcol] = A[jAcol + 12 * jp];
        A[jAcol + 12 * jp] = temp;
      }
    }
  }
}

/*
 * Arguments    : float A[72]
 *                const float B[36]
 * Return Type  : void
 */
static void c_mrdivide(float A[72], const float B[36])
{
  float b_A[36];
  signed char ipiv[6];
  int k;
  int j;
  int c;
  int kBcol;
  int ix;
  float temp;
  float s;
  int i;
  int jp;
  int jAcol;
  memcpy(&b_A[0], &B[0], 36U * sizeof(float));
  for (k = 0; k < 6; k++) {
    ipiv[k] = (signed char)(1 + k);
  }

  for (j = 0; j < 5; j++) {
    c = j * 7;
    kBcol = 0;
    ix = c;
    temp = (real32_T)fabs(b_A[c]);
    for (k = 2; k <= 6 - j; k++) {
      ix++;
      s = (real32_T)fabs(b_A[ix]);
      if (s > temp) {
        kBcol = k - 1;
        temp = s;
      }
    }

    if (b_A[c + kBcol] != 0.0F) {
      if (kBcol != 0) {
        ipiv[j] = (signed char)((j + kBcol) + 1);
        ix = j;
        kBcol += j;
        for (k = 0; k < 6; k++) {
          temp = b_A[ix];
          b_A[ix] = b_A[kBcol];
          b_A[kBcol] = temp;
          ix += 6;
          kBcol += 6;
        }
      }

      k = (c - j) + 6;
      for (i = c + 1; i + 1 <= k; i++) {
        b_A[i] /= b_A[c];
      }
    }

    jp = c;
    jAcol = c + 6;
    for (kBcol = 1; kBcol <= 5 - j; kBcol++) {
      temp = b_A[jAcol];
      if (b_A[jAcol] != 0.0F) {
        ix = c + 1;
        k = (jp - j) + 12;
        for (i = 7 + jp; i + 1 <= k; i++) {
          b_A[i] += b_A[ix] * -temp;
          ix++;
        }
      }

      jAcol += 6;
      jp += 6;
    }
  }

  for (j = 0; j < 6; j++) {
    jp = 12 * j;
    jAcol = 6 * j;
    for (k = 1; k <= j; k++) {
      kBcol = 12 * (k - 1);
      if (b_A[(k + jAcol) - 1] != 0.0F) {
        for (i = 0; i < 12; i++) {
          A[i + jp] -= b_A[(k + jAcol) - 1] * A[i + kBcol];
        }
      }
    }

    temp = 1.0F / b_A[j + jAcol];
    for (i = 0; i < 12; i++) {
      A[i + jp] *= temp;
    }
  }

  for (j = 5; j > -1; j += -1) {
    jp = 12 * j;
    jAcol = 6 * j - 1;
    for (k = j + 2; k < 7; k++) {
      kBcol = 12 * (k - 1);
      if (b_A[k + jAcol] != 0.0F) {
        for (i = 0; i < 12; i++) {
          A[i + jp] -= b_A[k + jAcol] * A[i + kBcol];
        }
      }
    }
  }

  for (kBcol = 4; kBcol > -1; kBcol += -1) {
    if (ipiv[kBcol] != kBcol + 1) {
      jp = ipiv[kBcol] - 1;
      for (jAcol = 0; jAcol < 12; jAcol++) {
        temp = A[jAcol + 12 * kBcol];
        A[jAcol + 12 * kBcol] = A[jAcol + 12 * jp];
        A[jAcol + 12 * jp] = temp;
      }
    }
  }
}

/*
 * Arguments    : const float a[3]
 *                const float b[3]
 *                float c[3]
 * Return Type  : void
 */
static void cross(const float a[3], const float b[3], float c[3])
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

/*
 * Arguments    : const float v[12]
 *                float d[144]
 * Return Type  : void
 */
static void diag(const float v[12], float d[144])
{
  int j;
  memset(&d[0], 0, 144U * sizeof(float));
  for (j = 0; j < 12; j++) {
    d[j + 12 * j] = v[j];
  }
}

/*
 * Arguments    : double I[144]
 * Return Type  : void
 */
static void eye(double I[144])
{
  int k;
  memset(&I[0], 0, 144U * sizeof(double));
  for (k = 0; k < 12; k++) {
    I[k + 12 * k] = 1.0;
  }
}

/*
 * Arguments    : const float x[9]
 *                float y[9]
 * Return Type  : void
 */
static void inv(const float x[9], float y[9])
{
  float b_x[9];
  int p1;
  int p2;
  int p3;
  float absx11;
  float absx21;
  float absx31;
  int itmp;
  float b_y;
  for (p1 = 0; p1 < 9; p1++) {
    b_x[p1] = x[p1];
  }

  p1 = 0;
  p2 = 3;
  p3 = 6;
  absx11 = (real32_T)fabs(x[0]);
  absx21 = (real32_T)fabs(x[1]);
  absx31 = (real32_T)fabs(x[2]);
  if ((absx21 > absx11) && (absx21 > absx31)) {
    p1 = 3;
    p2 = 0;
    b_x[0] = x[1];
    b_x[1] = x[0];
    b_x[3] = x[4];
    b_x[4] = x[3];
    b_x[6] = x[7];
    b_x[7] = x[6];
  } else {
    if (absx31 > absx11) {
      p1 = 6;
      p3 = 0;
      b_x[0] = x[2];
      b_x[2] = x[0];
      b_x[3] = x[5];
      b_x[5] = x[3];
      b_x[6] = x[8];
      b_x[8] = x[6];
    }
  }

  absx21 = b_x[1] / b_x[0];
  b_x[1] /= b_x[0];
  absx11 = b_x[2] / b_x[0];
  b_x[2] /= b_x[0];
  b_x[4] -= absx21 * b_x[3];
  b_x[5] -= absx11 * b_x[3];
  b_x[7] -= absx21 * b_x[6];
  b_x[8] -= absx11 * b_x[6];
  if ((real32_T)fabs(b_x[5]) > (real32_T)fabs(b_x[4])) {
    itmp = p2;
    p2 = p3;
    p3 = itmp;
    b_x[1] = absx11;
    b_x[2] = absx21;
    absx11 = b_x[4];
    b_x[4] = b_x[5];
    b_x[5] = absx11;
    absx11 = b_x[7];
    b_x[7] = b_x[8];
    b_x[8] = absx11;
  }

  absx31 = b_x[5];
  b_y = b_x[4];
  absx21 = b_x[5] / b_x[4];
  b_x[8] -= absx21 * b_x[7];
  absx11 = (absx21 * b_x[1] - b_x[2]) / b_x[8];
  absx21 = -(b_x[1] + b_x[7] * absx11) / b_x[4];
  y[p1] = ((1.0F - b_x[3] * absx21) - b_x[6] * absx11) / b_x[0];
  y[p1 + 1] = absx21;
  y[p1 + 2] = absx11;
  absx11 = -(absx31 / b_y) / b_x[8];
  absx21 = (1.0F - b_x[7] * absx11) / b_x[4];
  y[p2] = -(b_x[3] * absx21 + b_x[6] * absx11) / b_x[0];
  y[p2 + 1] = absx21;
  y[p2 + 2] = absx11;
  absx11 = 1.0F / b_x[8];
  absx21 = -b_x[7] * absx11 / b_x[4];
  y[p3] = -(b_x[3] * absx21 + b_x[6] * absx11) / b_x[0];
  y[p3 + 1] = absx21;
  y[p3 + 2] = absx11;
}

/*
 * Arguments    : const float a[9]
 *                float c[9]
 * Return Type  : void
 */
static void mpower(const float a[9], float c[9])
{
  int i2;
  int i3;
  int i4;
  for (i2 = 0; i2 < 3; i2++) {
    for (i3 = 0; i3 < 3; i3++) {
      c[i2 + 3 * i3] = 0.0F;
      for (i4 = 0; i4 < 3; i4++) {
        c[i2 + 3 * i3] += a[i2 + 3 * i4] * a[i4 + 3 * i3];
      }
    }
  }
}

/*
 * Arguments    : const float A[36]
 *                const float B[9]
 *                float y[36]
 * Return Type  : void
 */
static void mrdivide(const float A[36], const float B[9], float y[36])
{
  float b_A[9];
  int rtemp;
  int r1;
  int r2;
  int r3;
  float maxval;
  float a21;
  for (rtemp = 0; rtemp < 9; rtemp++) {
    b_A[rtemp] = B[rtemp];
  }

  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = (real32_T)fabs(B[0]);
  a21 = (real32_T)fabs(B[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if ((real32_T)fabs(B[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  b_A[r2] = B[r2] / B[r1];
  b_A[r3] /= b_A[r1];
  b_A[3 + r2] -= b_A[r2] * b_A[3 + r1];
  b_A[3 + r3] -= b_A[r3] * b_A[3 + r1];
  b_A[6 + r2] -= b_A[r2] * b_A[6 + r1];
  b_A[6 + r3] -= b_A[r3] * b_A[6 + r1];
  if ((real32_T)fabs(b_A[3 + r3]) > (real32_T)fabs(b_A[3 + r2])) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }

  b_A[3 + r3] /= b_A[3 + r2];
  b_A[6 + r3] -= b_A[3 + r3] * b_A[6 + r2];
  for (rtemp = 0; rtemp < 12; rtemp++) {
    y[rtemp + 12 * r1] = A[rtemp] / b_A[r1];
    y[rtemp + 12 * r2] = A[12 + rtemp] - y[rtemp + 12 * r1] * b_A[3 + r1];
    y[rtemp + 12 * r3] = A[24 + rtemp] - y[rtemp + 12 * r1] * b_A[6 + r1];
    y[rtemp + 12 * r2] /= b_A[3 + r2];
    y[rtemp + 12 * r3] -= y[rtemp + 12 * r2] * b_A[6 + r2];
    y[rtemp + 12 * r3] /= b_A[6 + r3];
    y[rtemp + 12 * r2] -= y[rtemp + 12 * r3] * b_A[3 + r3];
    y[rtemp + 12 * r1] -= y[rtemp + 12 * r3] * b_A[r3];
    y[rtemp + 12 * r1] -= y[rtemp + 12 * r2] * b_A[r2];
  }
}

/*
 * Arguments    : const float x[3]
 * Return Type  : float
 */
static float norm(const float x[3])
{
  float y;
  float scale;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  scale = 1.17549435E-38F;
  for (k = 0; k < 3; k++) {
    absxk = (real32_T)fabs(x[k]);
    if (absxk > scale) {
      t = scale / absxk;
      y = 1.0F + y * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * (real32_T)sqrt(y);
}

/*
 * Arguments    : const float x[3]
 *                float y
 *                float z[3]
 * Return Type  : void
 */
static void rdivide(const float x[3], float y, float z[3])
{
  int i;
  for (i = 0; i < 3; i++) {
    z[i] = x[i] / y;
  }
}

/*
 * Arguments    : float u0
 *                float u1
 * Return Type  : float
 */
static float rt_atan2f_snf(float u0, float u1)
{
  float y;
  int b_u0;
  int b_u1;
  if (rtIsNaNF(u0) || rtIsNaNF(u1)) {
    y = ((real32_T)rtNaN);
  } else if (rtIsInfF(u0) && rtIsInfF(u1)) {
    if (u0 > 0.0F) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0F) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = (real32_T)atan2((float)b_u0, (float)b_u1);
  } else if (u1 == 0.0F) {
    if (u0 > 0.0F) {
      y = RT_PIF / 2.0F;
    } else if (u0 < 0.0F) {
      y = -(float)(RT_PIF / 2.0F);
    } else {
      y = 0.0F;
    }
  } else {
    y = (real32_T)atan2(u0, u1);
  }

  return y;
}

/*
 * %
 *  Arguments:
 *  approx_order: 0-1st  1-2nd
 *  use_J: set to true if you have the inertia matrix J for your quadrotor
 *  xa_apo_k: old state vectotr
 *  zFlag: if sensor measurement is available [gyro, acc, mag]
 *  dt: dt in s
 *  z: measurements [gyro, acc, mag]
 *  q_vector = [q_rotSpeed, q_rotAcc, q_acc, q_mag]';
 *  process noise: gyro, gyro acceleration, acceleration,magnetometer
 *  r_vector = [r_gyro, r_accel, r_mag]';
 *  measurement noise: gyro,accel,mag
 *  J: moment of inertia matrix
 * Arguments    : unsigned char approx_order
 *                unsigned char use_J
 *                const float X_apo_k[12]
 *                const float P_apo_k[144]
 *                const unsigned char zFlag[3]
 *                float dt
 *                const float z[9]
 *                const float q_vector[4]
 *                const float r_vector[3]
 *                const float J[9]
 *                float xa_apo[12]
 *                float Pa_apo[144]
 *                float Rot_matrix[9]
 *                float eulerAngles[3]
 * Return Type  : void
 */
void AttitudeEKF_0916(unsigned char approx_order, unsigned char use_J, const
                      float X_apo_k[12], const float P_apo_k[144], const
                      unsigned char zFlag[3], float dt, const float z[9], const
                      float q_vector[4], const float r_vector[3], const float J
                      [9], float xa_apo[12], float Pa_apo[144], float
                      Rot_matrix[9], float eulerAngles[3])
{
  float Ji[9];
  float b_q_vector[12];
  static float Q[144];
  float b_X_apo_k[3];
  float b_J[3];
  float c_X_apo_k[3];
  int i;
  int i0;
  float zek[3];
  float d_X_apo_k[3];
  float wak[3];
  float y;
  float b_Ji[9];
  float fv0[9];
  static const signed char iv0[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  float e_X_apo_k[3];
  float f_X_apo_k[3];
  float muk[3];
  float x_apr[12];
  static double dv0[144];
  float fv1[9];
  float fv2[9];
  static float A_lin[144];
  static const signed char iv1[36] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  static float b_A_lin[144];
  int i1;
  static float P_apr[144];
  static float a[108];
  static const signed char b_a[108] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  static float S_k[81];
  static const signed char b[108] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  float b_r_vector[9];
  static float K_k[108];
  float c_a[36];
  static const signed char d_a[36] = { 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  static const signed char b_b[36] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  float c_r_vector[3];
  float b_P_apr[36];
  float b_K_k[36];
  static float e_a[72];
  static const signed char f_a[72] = { 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0 };

  float b_S_k[36];
  static const signed char c_b[72] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0 };

  float d_r_vector[6];
  float c_S_k[6];
  float c_K_k[72];
  static const signed char g_a[72] = { 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1 };

  static const signed char d_b[72] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1 };

  float b_z[6];

  /*  Output: */
  /*  xa_apo: updated state vectotr */
  /*  Pa_apo: updated state covariance matrix */
  /*  Rot_matrix: rotation matrix */
  /*  eulerAngles: euler angles */
  /* % init */
  inv(J, Ji);

  /* % copy the states */
  /*  x  body angular rate */
  /*  y  body angular rate */
  /*  z  body angular rate */
  /*  x  body angular acceleration */
  /*  y  body angular acceleration */
  /*  z  body angular acceleration */
  /*  x  component gravity vector */
  /*  y  component gravity vector */
  /*  z  component gravity vector */
  /*  x  component magnetic field vector */
  /*  y  component magnetic field vector */
  /*  z  component magnetic field vector */
  /* % copy the Q & R */
  /*  process covariance matrix  */
  b_q_vector[0] = q_vector[0];
  b_q_vector[1] = q_vector[0];
  b_q_vector[2] = q_vector[0];
  b_q_vector[3] = q_vector[1];
  b_q_vector[4] = q_vector[1];
  b_q_vector[5] = q_vector[1];
  b_q_vector[6] = q_vector[2];
  b_q_vector[7] = q_vector[2];
  b_q_vector[8] = q_vector[2];
  b_q_vector[9] = q_vector[3];
  b_q_vector[10] = q_vector[3];
  b_q_vector[11] = q_vector[3];
  diag(b_q_vector, Q);

  /* % prediction */
  /*  compute the apriori state estimate from the previous aposteriori estimate */
  /*  body angular accelerations */
  if (use_J == 1) {
    b_X_apo_k[0] = X_apo_k[3];
    b_X_apo_k[1] = X_apo_k[4];
    b_X_apo_k[2] = X_apo_k[5];
    c_X_apo_k[0] = X_apo_k[3];
    c_X_apo_k[1] = X_apo_k[4];
    c_X_apo_k[2] = X_apo_k[5];
    for (i = 0; i < 3; i++) {
      b_J[i] = 0.0F;
      for (i0 = 0; i0 < 3; i0++) {
        b_J[i] += J[i + 3 * i0] * c_X_apo_k[i0];
      }
    }

    cross(b_X_apo_k, b_J, zek);
    for (i = 0; i < 3; i++) {
      c_X_apo_k[i] = -zek[i];
    }

    d_X_apo_k[0] = X_apo_k[3];
    d_X_apo_k[1] = X_apo_k[4];
    d_X_apo_k[2] = X_apo_k[5];
    for (i = 0; i < 3; i++) {
      b_X_apo_k[i] = 0.0F;
      for (i0 = 0; i0 < 3; i0++) {
        b_X_apo_k[i] += Ji[i + 3 * i0] * c_X_apo_k[i0];
      }

      wak[i] = d_X_apo_k[i] + b_X_apo_k[i] * dt;
    }
  } else {
    wak[0] = X_apo_k[3];
    wak[1] = X_apo_k[4];
    wak[2] = X_apo_k[5];
  }

  /* body angular rates */
  /* derivative of the prediction rotation matrix */
  Ji[0] = 0.0F;
  Ji[1] = -X_apo_k[2];
  Ji[2] = X_apo_k[1];
  Ji[3] = X_apo_k[2];
  Ji[4] = 0.0F;
  Ji[5] = -X_apo_k[0];
  Ji[6] = -X_apo_k[1];
  Ji[7] = X_apo_k[0];
  Ji[8] = 0.0F;

  /* prediction of the earth z vector */
  if (approx_order == 1) {
    /* e^(Odt)=I+dt*O+dt^2/2!O^2 */
    y = dt * dt / 2.0F;
    mpower(Ji, b_Ji);
    for (i = 0; i < 3; i++) {
      for (i0 = 0; i0 < 3; i0++) {
        fv0[i0 + 3 * i] = ((float)iv0[i0 + 3 * i] + Ji[i0 + 3 * i] * dt) + y *
          b_Ji[i0 + 3 * i];
      }
    }

    e_X_apo_k[0] = X_apo_k[6];
    e_X_apo_k[1] = X_apo_k[7];
    e_X_apo_k[2] = X_apo_k[8];
    for (i = 0; i < 3; i++) {
      zek[i] = 0.0F;
      for (i0 = 0; i0 < 3; i0++) {
        zek[i] += fv0[i + 3 * i0] * e_X_apo_k[i0];
      }
    }
  } else {
    /*  so we do a first order approximation of the exponential map */
    for (i = 0; i < 3; i++) {
      for (i0 = 0; i0 < 3; i0++) {
        b_Ji[i0 + 3 * i] = Ji[i0 + 3 * i] * dt + (float)iv0[i0 + 3 * i];
      }
    }

    e_X_apo_k[0] = X_apo_k[6];
    e_X_apo_k[1] = X_apo_k[7];
    e_X_apo_k[2] = X_apo_k[8];
    for (i = 0; i < 3; i++) {
      zek[i] = 0.0F;
      for (i0 = 0; i0 < 3; i0++) {
        zek[i] += b_Ji[i + 3 * i0] * e_X_apo_k[i0];
      }
    }

    /* zek =expm2(O*dt)*[zex;zey;zez]; not working because use double */
  }

  /* prediction of the magnetic vector */
  if (approx_order == 1) {
    /* e^(Odt)=I+dt*O+dt^2/2!O^2 */
    y = dt * dt / 2.0F;
    mpower(Ji, b_Ji);
    for (i = 0; i < 3; i++) {
      for (i0 = 0; i0 < 3; i0++) {
        fv0[i0 + 3 * i] = ((float)iv0[i0 + 3 * i] + Ji[i0 + 3 * i] * dt) + y *
          b_Ji[i0 + 3 * i];
      }
    }

    f_X_apo_k[0] = X_apo_k[9];
    f_X_apo_k[1] = X_apo_k[10];
    f_X_apo_k[2] = X_apo_k[11];
    for (i = 0; i < 3; i++) {
      muk[i] = 0.0F;
      for (i0 = 0; i0 < 3; i0++) {
        muk[i] += fv0[i + 3 * i0] * f_X_apo_k[i0];
      }
    }
  } else {
    /*  so we do a first order approximation of the exponential map */
    for (i = 0; i < 3; i++) {
      for (i0 = 0; i0 < 3; i0++) {
        b_Ji[i0 + 3 * i] = Ji[i0 + 3 * i] * dt + (float)iv0[i0 + 3 * i];
      }
    }

    f_X_apo_k[0] = X_apo_k[9];
    f_X_apo_k[1] = X_apo_k[10];
    f_X_apo_k[2] = X_apo_k[11];
    for (i = 0; i < 3; i++) {
      muk[i] = 0.0F;
      for (i0 = 0; i0 < 3; i0++) {
        muk[i] += b_Ji[i + 3 * i0] * f_X_apo_k[i0];
      }
    }

    /* muk =expm2(O*dt)*[mux;muy;muz]; not working because use double */
  }

  /*  x prediction Ê±¼ä´«µÝ */
  x_apr[0] = X_apo_k[0] + dt * wak[0];
  x_apr[1] = X_apo_k[1] + dt * wak[1];
  x_apr[2] = X_apo_k[2] + dt * wak[2];
  for (i = 0; i < 3; i++) {
    x_apr[i + 3] = wak[i];
  }

  for (i = 0; i < 3; i++) {
    x_apr[i + 6] = zek[i];
  }

  for (i = 0; i < 3; i++) {
    x_apr[i + 9] = muk[i];
  }

  /*  compute the apriori error covariance estimate from the previous */
  /* aposteriori estimate */
  eye(dv0);
  fv1[0] = 0.0F;
  fv1[1] = -X_apo_k[8];
  fv1[2] = X_apo_k[7];
  fv1[3] = X_apo_k[8];
  fv1[4] = 0.0F;
  fv1[5] = -X_apo_k[6];
  fv1[6] = -X_apo_k[7];
  fv1[7] = X_apo_k[6];
  fv1[8] = 0.0F;
  fv2[0] = 0.0F;
  fv2[1] = -X_apo_k[11];
  fv2[2] = X_apo_k[10];
  fv2[3] = X_apo_k[11];
  fv2[4] = 0.0F;
  fv2[5] = -X_apo_k[9];
  fv2[6] = -X_apo_k[10];
  fv2[7] = X_apo_k[9];
  fv2[8] = 0.0F;
  for (i = 0; i < 12; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      A_lin[i0 + 12 * i] = iv1[i0 + 3 * i];
    }

    for (i0 = 0; i0 < 3; i0++) {
      A_lin[(i0 + 12 * i) + 3] = 0.0F;
    }
  }

  for (i = 0; i < 3; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      A_lin[(i0 + 12 * i) + 6] = -fv1[i0 + 3 * i];
    }
  }

  for (i = 0; i < 3; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      A_lin[(i0 + 12 * (i + 3)) + 6] = 0.0F;
    }
  }

  for (i = 0; i < 3; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      A_lin[(i0 + 12 * (i + 6)) + 6] = Ji[i0 + 3 * i];
    }
  }

  for (i = 0; i < 3; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      A_lin[(i0 + 12 * (i + 9)) + 6] = 0.0F;
    }
  }

  for (i = 0; i < 3; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      A_lin[(i0 + 12 * i) + 9] = -fv2[i0 + 3 * i];
    }
  }

  for (i = 0; i < 3; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      A_lin[(i0 + 12 * (i + 3)) + 9] = 0.0F;
    }
  }

  for (i = 0; i < 3; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      A_lin[(i0 + 12 * (i + 6)) + 9] = 0.0F;
    }
  }

  for (i = 0; i < 3; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      A_lin[(i0 + 12 * (i + 9)) + 9] = Ji[i0 + 3 * i];
    }
  }

  for (i = 0; i < 12; i++) {
    for (i0 = 0; i0 < 12; i0++) {
      b_A_lin[i0 + 12 * i] = (float)dv0[i0 + 12 * i] + A_lin[i0 + 12 * i] * dt;
    }
  }

  for (i = 0; i < 12; i++) {
    for (i0 = 0; i0 < 12; i0++) {
      A_lin[i + 12 * i0] = 0.0F;
      for (i1 = 0; i1 < 12; i1++) {
        A_lin[i + 12 * i0] += b_A_lin[i + 12 * i1] * P_apo_k[i1 + 12 * i0];
      }
    }
  }

  for (i = 0; i < 12; i++) {
    for (i0 = 0; i0 < 12; i0++) {
      y = 0.0F;
      for (i1 = 0; i1 < 12; i1++) {
        y += A_lin[i + 12 * i1] * b_A_lin[i0 + 12 * i1];
      }

      P_apr[i + 12 * i0] = y + Q[i + 12 * i0];
    }
  }

  /* % update */
  if ((zFlag[0] == 1) && (zFlag[1] == 1) && (zFlag[2] == 1)) {
    /* observation matrix */
    /* [zw;ze;zmk]; */
    /* S_k=H_k*P_apr*H_k'+R; */
    for (i = 0; i < 9; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        a[i + 9 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          a[i + 9 * i0] += (float)b_a[i + 9 * i1] * P_apr[i1 + 12 * i0];
        }
      }

      for (i0 = 0; i0 < 9; i0++) {
        S_k[i + 9 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          S_k[i + 9 * i0] += a[i + 9 * i1] * (float)b[i1 + 12 * i0];
        }
      }
    }

    b_r_vector[0] = r_vector[0];
    b_r_vector[1] = r_vector[0];
    b_r_vector[2] = r_vector[0];
    b_r_vector[3] = r_vector[1];
    b_r_vector[4] = r_vector[1];
    b_r_vector[5] = r_vector[1];
    b_r_vector[6] = r_vector[2];
    b_r_vector[7] = r_vector[2];
    b_r_vector[8] = r_vector[2];
    for (i = 0; i < 9; i++) {
      b_Ji[i] = S_k[10 * i] + b_r_vector[i];
    }

    for (i = 0; i < 9; i++) {
      S_k[10 * i] = b_Ji[i];
    }

    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 9; i0++) {
        K_k[i + 12 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          K_k[i + 12 * i0] += P_apr[i + 12 * i1] * (float)b[i1 + 12 * i0];
        }
      }
    }

    b_mrdivide(K_k, S_k);
    for (i = 0; i < 9; i++) {
      y = 0.0F;
      for (i0 = 0; i0 < 12; i0++) {
        y += (float)b_a[i + 9 * i0] * x_apr[i0];
      }

      b_Ji[i] = z[i] - y;
    }

    for (i = 0; i < 12; i++) {
      y = 0.0F;
      for (i0 = 0; i0 < 9; i0++) {
        y += K_k[i + 12 * i0] * b_Ji[i0];
      }

      xa_apo[i] = x_apr[i] + y;
    }

    eye(dv0);
    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        y = 0.0F;
        for (i1 = 0; i1 < 9; i1++) {
          y += K_k[i + 12 * i1] * (float)b_a[i1 + 9 * i0];
        }

        Q[i + 12 * i0] = (float)dv0[i + 12 * i0] - y;
      }
    }

    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        Pa_apo[i + 12 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          Pa_apo[i + 12 * i0] += Q[i + 12 * i1] * P_apr[i1 + 12 * i0];
        }
      }
    }
  } else if ((zFlag[0] == 1) && (zFlag[1] == 0) && (zFlag[2] == 0)) {
    /* observation matrix */
    /*  S_k=H_k(1:3,1:12)*P_apr*H_k(1:3,1:12)'+R(1:3,1:3); */
    for (i = 0; i < 3; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        c_a[i + 3 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          c_a[i + 3 * i0] += (float)d_a[i + 3 * i1] * P_apr[i1 + 12 * i0];
        }
      }

      for (i0 = 0; i0 < 3; i0++) {
        Ji[i + 3 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          Ji[i + 3 * i0] += c_a[i + 3 * i1] * (float)b_b[i1 + 12 * i0];
        }
      }
    }

    c_r_vector[0] = r_vector[0];
    c_r_vector[1] = r_vector[0];
    c_r_vector[2] = r_vector[0];
    for (i = 0; i < 3; i++) {
      b_X_apo_k[i] = Ji[i << 2] + c_r_vector[i];
    }

    for (i = 0; i < 3; i++) {
      Ji[i << 2] = b_X_apo_k[i];
    }

    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 3; i0++) {
        b_P_apr[i + 12 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          b_P_apr[i + 12 * i0] += P_apr[i + 12 * i1] * (float)b_b[i1 + 12 * i0];
        }
      }
    }

    mrdivide(b_P_apr, Ji, b_K_k);
    for (i = 0; i < 3; i++) {
      y = 0.0F;
      for (i0 = 0; i0 < 12; i0++) {
        y += (float)d_a[i + 3 * i0] * x_apr[i0];
      }

      c_X_apo_k[i] = z[i] - y;
    }

    for (i = 0; i < 12; i++) {
      y = 0.0F;
      for (i0 = 0; i0 < 3; i0++) {
        y += b_K_k[i + 12 * i0] * c_X_apo_k[i0];
      }

      xa_apo[i] = x_apr[i] + y;
    }

    eye(dv0);
    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        y = 0.0F;
        for (i1 = 0; i1 < 3; i1++) {
          y += b_K_k[i + 12 * i1] * (float)d_a[i1 + 3 * i0];
        }

        Q[i + 12 * i0] = (float)dv0[i + 12 * i0] - y;
      }
    }

    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        Pa_apo[i + 12 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          Pa_apo[i + 12 * i0] += Q[i + 12 * i1] * P_apr[i1 + 12 * i0];
        }
      }
    }
  } else if ((zFlag[0] == 1) && (zFlag[1] == 1) && (zFlag[2] == 0)) {
    /* observation matrix */
    /*  S_k=H_k(1:6,1:12)*P_apr*H_k(1:6,1:12)'+R(1:6,1:6); */
    for (i = 0; i < 6; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        e_a[i + 6 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          e_a[i + 6 * i0] += (float)f_a[i + 6 * i1] * P_apr[i1 + 12 * i0];
        }
      }

      for (i0 = 0; i0 < 6; i0++) {
        b_S_k[i + 6 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          b_S_k[i + 6 * i0] += e_a[i + 6 * i1] * (float)c_b[i1 + 12 * i0];
        }
      }
    }

    d_r_vector[0] = r_vector[0];
    d_r_vector[1] = r_vector[0];
    d_r_vector[2] = r_vector[0];
    d_r_vector[3] = r_vector[1];
    d_r_vector[4] = r_vector[1];
    d_r_vector[5] = r_vector[1];
    for (i = 0; i < 6; i++) {
      c_S_k[i] = b_S_k[7 * i] + d_r_vector[i];
    }

    for (i = 0; i < 6; i++) {
      b_S_k[7 * i] = c_S_k[i];
    }

    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 6; i0++) {
        c_K_k[i + 12 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          c_K_k[i + 12 * i0] += P_apr[i + 12 * i1] * (float)c_b[i1 + 12 * i0];
        }
      }
    }

    c_mrdivide(c_K_k, b_S_k);
    for (i = 0; i < 6; i++) {
      y = 0.0F;
      for (i0 = 0; i0 < 12; i0++) {
        y += (float)f_a[i + 6 * i0] * x_apr[i0];
      }

      d_r_vector[i] = z[i] - y;
    }

    for (i = 0; i < 12; i++) {
      y = 0.0F;
      for (i0 = 0; i0 < 6; i0++) {
        y += c_K_k[i + 12 * i0] * d_r_vector[i0];
      }

      xa_apo[i] = x_apr[i] + y;
    }

    eye(dv0);
    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        y = 0.0F;
        for (i1 = 0; i1 < 6; i1++) {
          y += c_K_k[i + 12 * i1] * (float)f_a[i1 + 6 * i0];
        }

        Q[i + 12 * i0] = (float)dv0[i + 12 * i0] - y;
      }
    }

    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        Pa_apo[i + 12 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          Pa_apo[i + 12 * i0] += Q[i + 12 * i1] * P_apr[i1 + 12 * i0];
        }
      }
    }
  } else if ((zFlag[0] == 1) && (zFlag[1] == 0) && (zFlag[2] == 1)) {
    /* observation matrix */
    /* S_k=H_k(1:6,1:12)*P_apr*H_k(1:6,1:12)'+R(1:6,1:6); */
    for (i = 0; i < 6; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        e_a[i + 6 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          e_a[i + 6 * i0] += (float)g_a[i + 6 * i1] * P_apr[i1 + 12 * i0];
        }
      }

      for (i0 = 0; i0 < 6; i0++) {
        b_S_k[i + 6 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          b_S_k[i + 6 * i0] += e_a[i + 6 * i1] * (float)d_b[i1 + 12 * i0];
        }
      }
    }

    d_r_vector[0] = r_vector[0];
    d_r_vector[1] = r_vector[0];
    d_r_vector[2] = r_vector[0];
    d_r_vector[3] = r_vector[2];
    d_r_vector[4] = r_vector[2];
    d_r_vector[5] = r_vector[2];
    for (i = 0; i < 6; i++) {
      c_S_k[i] = b_S_k[7 * i] + d_r_vector[i];
    }

    for (i = 0; i < 6; i++) {
      b_S_k[7 * i] = c_S_k[i];
    }

    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 6; i0++) {
        c_K_k[i + 12 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          c_K_k[i + 12 * i0] += P_apr[i + 12 * i1] * (float)d_b[i1 + 12 * i0];
        }
      }
    }

    c_mrdivide(c_K_k, b_S_k);
    for (i = 0; i < 3; i++) {
      d_r_vector[i] = z[i];
    }

    for (i = 0; i < 3; i++) {
      d_r_vector[i + 3] = z[6 + i];
    }

    for (i = 0; i < 6; i++) {
      c_S_k[i] = 0.0F;
      for (i0 = 0; i0 < 12; i0++) {
        c_S_k[i] += (float)g_a[i + 6 * i0] * x_apr[i0];
      }

      b_z[i] = d_r_vector[i] - c_S_k[i];
    }

    for (i = 0; i < 12; i++) {
      y = 0.0F;
      for (i0 = 0; i0 < 6; i0++) {
        y += c_K_k[i + 12 * i0] * b_z[i0];
      }

      xa_apo[i] = x_apr[i] + y;
    }

    eye(dv0);
    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        y = 0.0F;
        for (i1 = 0; i1 < 6; i1++) {
          y += c_K_k[i + 12 * i1] * (float)g_a[i1 + 6 * i0];
        }

        Q[i + 12 * i0] = (float)dv0[i + 12 * i0] - y;
      }
    }

    for (i = 0; i < 12; i++) {
      for (i0 = 0; i0 < 12; i0++) {
        Pa_apo[i + 12 * i0] = 0.0F;
        for (i1 = 0; i1 < 12; i1++) {
          Pa_apo[i + 12 * i0] += Q[i + 12 * i1] * P_apr[i1 + 12 * i0];
        }
      }
    }
  } else {
    for (i = 0; i < 12; i++) {
      xa_apo[i] = x_apr[i];
    }

    memcpy(&Pa_apo[0], &P_apr[0], 144U * sizeof(float));
  }

  /* % euler anglels extraction */
  for (i = 0; i < 3; i++) {
    c_X_apo_k[i] = -xa_apo[i + 6];
  }

  rdivide(c_X_apo_k, norm(*(float (*)[3])&xa_apo[6]), muk);
  rdivide(*(float (*)[3])&xa_apo[9], norm(*(float (*)[3])&xa_apo[9]), zek);
  cross(muk, zek, wak);
  for (i = 0; i < 3; i++) {
    c_X_apo_k[i] = wak[i];
  }

  rdivide(c_X_apo_k, norm(wak), wak);
  cross(wak, muk, zek);
  for (i = 0; i < 3; i++) {
    c_X_apo_k[i] = zek[i];
  }

  rdivide(c_X_apo_k, norm(zek), zek);

  /*  rotation matrix from earth to body system */
  for (i = 0; i < 3; i++) {
    Rot_matrix[i] = zek[i];
    Rot_matrix[3 + i] = wak[i];
    Rot_matrix[6 + i] = muk[i];
  }

  eulerAngles[0] = rt_atan2f_snf(Rot_matrix[7], Rot_matrix[8]);
  eulerAngles[1] = -(real32_T)asin(Rot_matrix[6]);
  eulerAngles[2] = rt_atan2f_snf(Rot_matrix[3], Rot_matrix[0]);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void AttitudeEKF_0916_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void AttitudeEKF_0916_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for AttitudeEKF_0916.c
 *
 * [EOF]
 */
