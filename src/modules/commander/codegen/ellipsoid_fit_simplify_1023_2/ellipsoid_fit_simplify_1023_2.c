/*
 * File: ellipsoid_fit_simplify_1023_2.c
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 24-Oct-2015 14:16:57
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "ellipsoid_fit_simplify_1023_2.h"

/* Function Declarations */
static void b_eml_matlab_zlartg(const creal32_T f, const creal32_T g, float *cs,
  creal32_T *sn);
static void b_sqrt(creal32_T *x);
static int div_s32(int numerator, int denominator);
static void eml_matlab_zggev(creal32_T A[9], int *info, creal32_T alpha1[3],
  creal32_T beta1[3], creal32_T V[9]);
static void eml_matlab_zhgeqz(creal32_T A[9], int ilo, int ihi, creal32_T Z[9],
  int *info, creal32_T alpha1[3], creal32_T beta1[3]);
static void eml_matlab_zlartg(const creal32_T f, const creal32_T g, float *cs,
  creal32_T *sn, creal32_T *r);
static void eml_matlab_ztgevc(const creal32_T A[9], creal32_T V[9]);
static void emxFree_real32_T(emxArray_real32_T **pEmxArray);
static void emxInit_real32_T(emxArray_real32_T **pEmxArray, int numDimensions);
static void mldivide(const float A[81], float B[9]);
static float rt_hypotf_snf(float u0, float u1);

/* Function Definitions */

/*
 * Arguments    : const creal32_T f
 *                const creal32_T g
 *                float *cs
 *                creal32_T *sn
 * Return Type  : void
 */
static void b_eml_matlab_zlartg(const creal32_T f, const creal32_T g, float *cs,
  creal32_T *sn)
{
  float scale;
  float f2s;
  float g2;
  float fs_re;
  float fs_im;
  float gs_re;
  float gs_im;
  boolean_T guard1 = false;
  float g2s;
  scale = (real32_T)fabs(f.re);
  f2s = (real32_T)fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  f2s = (real32_T)fabs(g.re);
  g2 = (real32_T)fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = false;
  if (scale >= 5.49755814E+11F) {
    do {
      fs_re *= 1.8189894E-12F;
      fs_im *= 1.8189894E-12F;
      gs_re *= 1.8189894E-12F;
      gs_im *= 1.8189894E-12F;
      scale *= 1.8189894E-12F;
    } while (!(scale < 5.49755814E+11F));

    guard1 = true;
  } else if (scale <= 1.8189894E-12F) {
    if ((g.re == 0.0F) && (g.im == 0.0F)) {
      *cs = 1.0F;
      sn->re = 0.0F;
      sn->im = 0.0F;
    } else {
      do {
        fs_re *= 5.49755814E+11F;
        fs_im *= 5.49755814E+11F;
        gs_re *= 5.49755814E+11F;
        gs_im *= 5.49755814E+11F;
        scale *= 5.49755814E+11F;
      } while (!(scale > 1.8189894E-12F));

      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0F > g2) {
      f2s = 1.0F;
    }

    if (scale <= f2s * 1.97215226E-31F) {
      if ((f.re == 0.0F) && (f.im == 0.0F)) {
        *cs = 0.0F;
        scale = rt_hypotf_snf(gs_re, gs_im);
        sn->re = gs_re / scale;
        sn->im = -gs_im / scale;
      } else {
        g2s = (real32_T)sqrt(g2);
        *cs = rt_hypotf_snf(fs_re, fs_im) / g2s;
        f2s = (real32_T)fabs(f.re);
        g2 = (real32_T)fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0F) {
          scale = rt_hypotf_snf(f.re, f.im);
          fs_re = f.re / scale;
          fs_im = f.im / scale;
        } else {
          f2s = 5.49755814E+11F * f.re;
          g2 = 5.49755814E+11F * f.im;
          scale = rt_hypotf_snf(f2s, g2);
          fs_re = f2s / scale;
          fs_im = g2 / scale;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      f2s = (real32_T)sqrt(1.0F + g2 / scale);
      *cs = 1.0F / f2s;
      scale += g2;
      fs_re = f2s * fs_re / scale;
      fs_im = f2s * fs_im / scale;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

/*
 * Arguments    : creal32_T *x
 * Return Type  : void
 */
static void b_sqrt(creal32_T *x)
{
  float absxi;
  float absxr;
  if (x->im == 0.0F) {
    if (x->re < 0.0F) {
      absxi = 0.0F;
      absxr = (real32_T)sqrt((real32_T)fabs(x->re));
    } else {
      absxi = (real32_T)sqrt(x->re);
      absxr = 0.0F;
    }
  } else if (x->re == 0.0F) {
    if (x->im < 0.0F) {
      absxi = (real32_T)sqrt(-x->im / 2.0F);
      absxr = -absxi;
    } else {
      absxi = (real32_T)sqrt(x->im / 2.0F);
      absxr = absxi;
    }
  } else if (rtIsNaNF(x->re) || rtIsNaNF(x->im)) {
    absxi = ((real32_T)rtNaN);
    absxr = ((real32_T)rtNaN);
  } else if (rtIsInfF(x->im)) {
    absxi = ((real32_T)rtInf);
    absxr = x->im;
  } else if (rtIsInfF(x->re)) {
    if (x->re < 0.0F) {
      absxi = 0.0F;
      absxr = ((real32_T)rtInf);
    } else {
      absxi = ((real32_T)rtInf);
      absxr = 0.0F;
    }
  } else {
    absxr = (real32_T)fabs(x->re);
    absxi = (real32_T)fabs(x->im);
    if ((absxr > 8.50705867E+37F) || (absxi > 8.50705867E+37F)) {
      absxr *= 0.5F;
      absxi *= 0.5F;
      absxi = rt_hypotf_snf(absxr, absxi);
      if (absxi > absxr) {
        absxi = (real32_T)sqrt(absxi) * (real32_T)sqrt(1.0F + absxr / absxi);
      } else {
        absxi = (real32_T)sqrt(absxi) * 1.41421354F;
      }
    } else {
      absxi = (real32_T)sqrt((rt_hypotf_snf(absxr, absxi) + absxr) * 0.5F);
    }

    if (x->re > 0.0F) {
      absxr = 0.5F * (x->im / absxi);
    } else {
      if (x->im < 0.0F) {
        absxr = -absxi;
      } else {
        absxr = absxi;
      }

      absxi = 0.5F * (x->im / absxr);
    }
  }

  x->re = absxi;
  x->im = absxr;
}

/*
 * Arguments    : int numerator
 *                int denominator
 * Return Type  : int
 */
static int div_s32(int numerator, int denominator)
{
  int quotient;
  unsigned int absNumerator;
  unsigned int absDenominator;
  boolean_T quotientNeedsNegation;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    if (numerator >= 0) {
      absNumerator = (unsigned int)numerator;
    } else {
      absNumerator = (unsigned int)-numerator;
    }

    if (denominator >= 0) {
      absDenominator = (unsigned int)denominator;
    } else {
      absDenominator = (unsigned int)-denominator;
    }

    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    absNumerator /= absDenominator;
    if (quotientNeedsNegation) {
      quotient = -(int)absNumerator;
    } else {
      quotient = (int)absNumerator;
    }
  }

  return quotient;
}

/*
 * Arguments    : creal32_T A[9]
 *                int *info
 *                creal32_T alpha1[3]
 *                creal32_T beta1[3]
 *                creal32_T V[9]
 * Return Type  : void
 */
static void eml_matlab_zggev(creal32_T A[9], int *info, creal32_T alpha1[3],
  creal32_T beta1[3], creal32_T V[9])
{
  float anrm;
  int ii;
  boolean_T exitg7;
  float absxk;
  int i;
  boolean_T ilascl;
  float anrmto;
  float ctoc;
  boolean_T notdone;
  float cfrom1;
  float stemp_re;
  float stemp_im;
  creal32_T b_A[9];
  int rscale[3];
  int ilo;
  int ihi;
  int32_T exitg2;
  int j;
  boolean_T exitg5;
  int nzcount;
  int jj;
  boolean_T exitg6;
  boolean_T guard2 = false;
  creal32_T atmp;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T guard1 = false;
  static const signed char iv0[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  *info = 0;
  anrm = 0.0F;
  ii = 0;
  exitg7 = false;
  while ((!exitg7) && (ii < 9)) {
    absxk = rt_hypotf_snf(A[ii].re, A[ii].im);
    if (rtIsNaNF(absxk)) {
      anrm = ((real32_T)rtNaN);
      exitg7 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      ii++;
    }
  }

  if (!((!rtIsInfF(anrm)) && (!rtIsNaNF(anrm)))) {
    for (i = 0; i < 3; i++) {
      alpha1[i].re = ((real32_T)rtNaN);
      alpha1[i].im = 0.0F;
      beta1[i].re = ((real32_T)rtNaN);
      beta1[i].im = 0.0F;
    }

    for (ii = 0; ii < 9; ii++) {
      V[ii].re = ((real32_T)rtNaN);
      V[ii].im = 0.0F;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    if ((anrm > 0.0F) && (anrm < 9.09494702E-13F)) {
      anrmto = 9.09494702E-13F;
      ilascl = true;
    } else {
      if (anrm > 1.09951163E+12F) {
        anrmto = 1.09951163E+12F;
        ilascl = true;
      }
    }

    if (ilascl) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = true;
      while (notdone) {
        cfrom1 = absxk * 1.97215226E-31F;
        stemp_re = ctoc / 5.0706024E+30F;
        if ((cfrom1 > ctoc) && (ctoc != 0.0F)) {
          stemp_im = 1.97215226E-31F;
          absxk = cfrom1;
        } else if (stemp_re > absxk) {
          stemp_im = 5.0706024E+30F;
          ctoc = stemp_re;
        } else {
          stemp_im = ctoc / absxk;
          notdone = false;
        }

        for (ii = 0; ii < 9; ii++) {
          A[ii].re *= stemp_im;
          A[ii].im *= stemp_im;
        }
      }
    }

    memcpy(&b_A[0], &A[0], 9U * sizeof(creal32_T));
    for (i = 0; i < 3; i++) {
      rscale[i] = 0;
    }

    ilo = 1;
    ihi = 3;
    do {
      exitg2 = 0;
      i = 0;
      j = 0;
      notdone = false;
      ii = ihi;
      exitg5 = false;
      while ((!exitg5) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = ihi;
        jj = 1;
        exitg6 = false;
        while ((!exitg6) && (jj <= ihi)) {
          guard2 = false;
          if ((b_A[(ii + 3 * (jj - 1)) - 1].re != 0.0F) || (b_A[(ii + 3 * (jj -
                 1)) - 1].im != 0.0F) || (ii == jj)) {
            if (nzcount == 0) {
              j = jj;
              nzcount = 1;
              guard2 = true;
            } else {
              nzcount = 2;
              exitg6 = true;
            }
          } else {
            guard2 = true;
          }

          if (guard2) {
            jj++;
          }
        }

        if (nzcount < 2) {
          notdone = true;
          exitg5 = true;
        } else {
          ii--;
        }
      }

      if (!notdone) {
        exitg2 = 2;
      } else {
        if (i != ihi) {
          for (ii = 0; ii < 3; ii++) {
            atmp = b_A[(i + 3 * ii) - 1];
            b_A[(i + 3 * ii) - 1] = b_A[(ihi + 3 * ii) - 1];
            b_A[(ihi + 3 * ii) - 1] = atmp;
          }
        }

        if (j != ihi) {
          for (ii = 0; ii + 1 <= ihi; ii++) {
            atmp = b_A[ii + 3 * (j - 1)];
            b_A[ii + 3 * (j - 1)] = b_A[ii + 3 * (ihi - 1)];
            b_A[ii + 3 * (ihi - 1)] = atmp;
          }
        }

        rscale[ihi - 1] = j;
        ihi--;
        if (ihi == 1) {
          rscale[0] = 1;
          exitg2 = 1;
        }
      }
    } while (exitg2 == 0);

    if (exitg2 == 1) {
    } else {
      do {
        exitg1 = 0;
        i = 0;
        j = 0;
        notdone = false;
        jj = ilo;
        exitg3 = false;
        while ((!exitg3) && (jj <= ihi)) {
          nzcount = 0;
          i = ihi;
          j = jj;
          ii = ilo;
          exitg4 = false;
          while ((!exitg4) && (ii <= ihi)) {
            guard1 = false;
            if ((b_A[(ii + 3 * (jj - 1)) - 1].re != 0.0F) || (b_A[(ii + 3 * (jj
                   - 1)) - 1].im != 0.0F) || (ii == jj)) {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                guard1 = true;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              guard1 = true;
            }

            if (guard1) {
              ii++;
            }
          }

          if (nzcount < 2) {
            notdone = true;
            exitg3 = true;
          } else {
            jj++;
          }
        }

        if (!notdone) {
          exitg1 = 1;
        } else {
          if (i != ilo) {
            for (ii = ilo - 1; ii + 1 < 4; ii++) {
              atmp = b_A[(i + 3 * ii) - 1];
              b_A[(i + 3 * ii) - 1] = b_A[(ilo + 3 * ii) - 1];
              b_A[(ilo + 3 * ii) - 1] = atmp;
            }
          }

          if (j != ilo) {
            for (ii = 0; ii + 1 <= ihi; ii++) {
              atmp = b_A[ii + 3 * (j - 1)];
              b_A[ii + 3 * (j - 1)] = b_A[ii + 3 * (ilo - 1)];
              b_A[ii + 3 * (ilo - 1)] = atmp;
            }
          }

          rscale[ilo - 1] = j;
          ilo++;
          if (ilo == ihi) {
            rscale[ilo - 1] = ilo;
            exitg1 = 1;
          }
        }
      } while (exitg1 == 0);
    }

    for (ii = 0; ii < 9; ii++) {
      V[ii].re = iv0[ii];
      V[ii].im = 0.0F;
    }

    if (ihi < ilo + 2) {
    } else {
      ii = ilo;
      while (ii < 2) {
        eml_matlab_zlartg(b_A[1], b_A[2], &cfrom1, &atmp, &b_A[1]);
        b_A[2].re = 0.0F;
        b_A[2].im = 0.0F;
        for (j = 1; j + 1 < 4; j++) {
          stemp_re = cfrom1 * b_A[1 + 3 * j].re + (atmp.re * b_A[2 + 3 * j].re -
            atmp.im * b_A[2 + 3 * j].im);
          stemp_im = cfrom1 * b_A[1 + 3 * j].im + (atmp.re * b_A[2 + 3 * j].im +
            atmp.im * b_A[2 + 3 * j].re);
          absxk = b_A[1 + 3 * j].im;
          ctoc = b_A[1 + 3 * j].re;
          b_A[2 + 3 * j].re = cfrom1 * b_A[2 + 3 * j].re - (atmp.re * b_A[1 + 3 *
            j].re + atmp.im * b_A[1 + 3 * j].im);
          b_A[2 + 3 * j].im = cfrom1 * b_A[2 + 3 * j].im - (atmp.re * absxk -
            atmp.im * ctoc);
          b_A[1 + 3 * j].re = stemp_re;
          b_A[1 + 3 * j].im = stemp_im;
        }

        atmp.re = -atmp.re;
        atmp.im = -atmp.im;
        for (i = ilo - 1; i + 1 < 4; i++) {
          stemp_re = cfrom1 * b_A[6 + i].re + (atmp.re * b_A[3 + i].re - atmp.im
            * b_A[3 + i].im);
          stemp_im = cfrom1 * b_A[6 + i].im + (atmp.re * b_A[3 + i].im + atmp.im
            * b_A[3 + i].re);
          absxk = b_A[6 + i].im;
          ctoc = b_A[6 + i].re;
          b_A[3 + i].re = cfrom1 * b_A[3 + i].re - (atmp.re * b_A[6 + i].re +
            atmp.im * b_A[6 + i].im);
          b_A[3 + i].im = cfrom1 * b_A[3 + i].im - (atmp.re * absxk - atmp.im *
            ctoc);
          b_A[6 + i].re = stemp_re;
          b_A[6 + i].im = stemp_im;
        }

        for (i = 0; i < 3; i++) {
          stemp_re = cfrom1 * V[6 + i].re + (atmp.re * V[3 + i].re - atmp.im *
            V[3 + i].im);
          stemp_im = cfrom1 * V[6 + i].im + (atmp.re * V[3 + i].im + atmp.im *
            V[3 + i].re);
          absxk = V[6 + i].im;
          ctoc = V[6 + i].re;
          V[3 + i].re = cfrom1 * V[3 + i].re - (atmp.re * V[6 + i].re + atmp.im *
            V[6 + i].im);
          V[3 + i].im = cfrom1 * V[3 + i].im - (atmp.re * absxk - atmp.im * ctoc);
          V[6 + i].re = stemp_re;
          V[6 + i].im = stemp_im;
        }

        ii = 2;
      }
    }

    eml_matlab_zhgeqz(b_A, ilo, ihi, V, info, alpha1, beta1);
    if (*info != 0) {
    } else {
      eml_matlab_ztgevc(b_A, V);
      if (ilo > 1) {
        for (i = ilo - 2; i + 1 >= 1; i--) {
          ii = rscale[i] - 1;
          if (rscale[i] != i + 1) {
            for (j = 0; j < 3; j++) {
              atmp = V[i + 3 * j];
              V[i + 3 * j] = V[ii + 3 * j];
              V[ii + 3 * j] = atmp;
            }
          }
        }
      }

      if (ihi < 3) {
        while (ihi + 1 < 4) {
          ii = rscale[ihi] - 1;
          if (rscale[ihi] != ihi + 1) {
            for (j = 0; j < 3; j++) {
              atmp = V[ihi + 3 * j];
              V[ihi + 3 * j] = V[ii + 3 * j];
              V[ii + 3 * j] = atmp;
            }
          }

          ihi++;
        }
      }

      for (ii = 0; ii < 3; ii++) {
        absxk = (real32_T)fabs(V[3 * ii].re) + (real32_T)fabs(V[3 * ii].im);
        for (nzcount = 0; nzcount < 2; nzcount++) {
          ctoc = (real32_T)fabs(V[(nzcount + 3 * ii) + 1].re) + (real32_T)fabs
            (V[(nzcount + 3 * ii) + 1].im);
          if (ctoc > absxk) {
            absxk = ctoc;
          }
        }

        if (absxk >= 9.09494702E-13F) {
          absxk = 1.0F / absxk;
          for (nzcount = 0; nzcount < 3; nzcount++) {
            V[nzcount + 3 * ii].re *= absxk;
            V[nzcount + 3 * ii].im *= absxk;
          }
        }
      }

      if (ilascl) {
        notdone = true;
        while (notdone) {
          cfrom1 = anrmto * 1.97215226E-31F;
          stemp_re = anrm / 5.0706024E+30F;
          if ((cfrom1 > anrm) && (anrm != 0.0F)) {
            stemp_im = 1.97215226E-31F;
            anrmto = cfrom1;
          } else if (stemp_re > anrmto) {
            stemp_im = 5.0706024E+30F;
            anrm = stemp_re;
          } else {
            stemp_im = anrm / anrmto;
            notdone = false;
          }

          for (ii = 0; ii < 3; ii++) {
            alpha1[ii].re *= stemp_im;
            alpha1[ii].im *= stemp_im;
          }
        }
      }
    }
  }
}

/*
 * Arguments    : creal32_T A[9]
 *                int ilo
 *                int ihi
 *                creal32_T Z[9]
 *                int *info
 *                creal32_T alpha1[3]
 *                creal32_T beta1[3]
 * Return Type  : void
 */
static void eml_matlab_zhgeqz(creal32_T A[9], int ilo, int ihi, creal32_T Z[9],
  int *info, creal32_T alpha1[3], creal32_T beta1[3])
{
  int i;
  float eshift_re;
  float eshift_im;
  creal32_T ctemp;
  float rho_re;
  float rho_im;
  float anorm;
  float sumsq;
  boolean_T firstNonZero;
  int j;
  int jp1;
  float reAij;
  float imAij;
  float temp1;
  float temp2;
  float b_atol;
  boolean_T failed;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  int ifirst;
  int istart;
  int ilast;
  int ilastm1;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T b_guard1 = false;
  creal32_T t1;
  creal32_T d;
  boolean_T exitg2;
  *info = 0;
  for (i = 0; i < 3; i++) {
    alpha1[i].re = 0.0F;
    alpha1[i].im = 0.0F;
    beta1[i].re = 1.0F;
    beta1[i].im = 0.0F;
  }

  eshift_re = 0.0F;
  eshift_im = 0.0F;
  ctemp.re = 0.0F;
  ctemp.im = 0.0F;
  rho_re = 0.0F;
  rho_im = 0.0F;
  anorm = 0.0F;
  if (ilo > ihi) {
  } else {
    anorm = 0.0F;
    sumsq = 0.0F;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      jp1 = j + 1;
      if (ihi < j + 1) {
        jp1 = ihi;
      }

      for (i = ilo; i <= jp1; i++) {
        reAij = A[(i + 3 * (j - 1)) - 1].re;
        imAij = A[(i + 3 * (j - 1)) - 1].im;
        if (reAij != 0.0F) {
          temp1 = (real32_T)fabs(reAij);
          if (firstNonZero) {
            sumsq = 1.0F;
            anorm = temp1;
            firstNonZero = false;
          } else if (anorm < temp1) {
            temp2 = anorm / temp1;
            sumsq = 1.0F + sumsq * temp2 * temp2;
            anorm = temp1;
          } else {
            temp2 = temp1 / anorm;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0F) {
          temp1 = (real32_T)fabs(imAij);
          if (firstNonZero) {
            sumsq = 1.0F;
            anorm = temp1;
            firstNonZero = false;
          } else if (anorm < temp1) {
            temp2 = anorm / temp1;
            sumsq = 1.0F + sumsq * temp2 * temp2;
            anorm = temp1;
          } else {
            temp2 = temp1 / anorm;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    anorm *= (real32_T)sqrt(sumsq);
  }

  reAij = 1.1920929E-7F * anorm;
  b_atol = 1.17549435E-38F;
  if (reAij > 1.17549435E-38F) {
    b_atol = reAij;
  }

  temp1 = 1.17549435E-38F;
  if (anorm > 1.17549435E-38F) {
    temp1 = anorm;
  }

  imAij = 1.0F / temp1;
  failed = true;
  for (j = ihi; j + 1 < 4; j++) {
    alpha1[j] = A[j + 3 * j];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1)) {
        if (ilast + 1 == ilo) {
          goto60 = true;
        } else if ((real32_T)fabs(A[ilast + 3 * ilastm1].re) + (real32_T)fabs
                   (A[ilast + 3 * ilastm1].im) <= b_atol) {
          A[ilast + 3 * ilastm1].re = 0.0F;
          A[ilast + 3 * ilastm1].im = 0.0F;
          goto60 = true;
        } else {
          j = ilastm1;
          exitg3 = false;
          while ((!exitg3) && (j + 1 >= ilo)) {
            if (j + 1 == ilo) {
              firstNonZero = true;
            } else if ((real32_T)fabs(A[j + 3 * (j - 1)].re) + (real32_T)fabs
                       (A[j + 3 * (j - 1)].im) <= b_atol) {
              A[j + 3 * (j - 1)].re = 0.0F;
              A[j + 3 * (j - 1)].im = 0.0F;
              firstNonZero = true;
            } else {
              firstNonZero = false;
            }

            if (firstNonZero) {
              ifirst = j + 1;
              goto70 = true;
              exitg3 = true;
            } else {
              j--;
            }
          }
        }

        if (goto60 || goto70) {
          firstNonZero = true;
        } else {
          firstNonZero = false;
        }

        if (!firstNonZero) {
          for (i = 0; i < 3; i++) {
            alpha1[i].re = ((real32_T)rtNaN);
            alpha1[i].im = 0.0F;
            beta1[i].re = ((real32_T)rtNaN);
            beta1[i].im = 0.0F;
          }

          for (jp1 = 0; jp1 < 9; jp1++) {
            Z[jp1].re = ((real32_T)rtNaN);
            Z[jp1].im = 0.0F;
          }

          *info = 1;
          exitg1 = 1;
        } else {
          b_guard1 = false;
          if (goto60) {
            goto60 = false;
            alpha1[ilast] = A[ilast + 3 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0F;
              eshift_im = 0.0F;
              b_guard1 = true;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              if (iiter - div_s32(iiter, 10) * 10 != 0) {
                reAij = -(A[ilast + 3 * ilast].re - A[ilastm1 + 3 * ilastm1].re);
                temp1 = -(A[ilast + 3 * ilast].im - A[ilastm1 + 3 * ilastm1].im);
                if (temp1 == 0.0F) {
                  t1.re = reAij / 2.0F;
                  t1.im = 0.0F;
                } else if (reAij == 0.0F) {
                  t1.re = 0.0F;
                  t1.im = temp1 / 2.0F;
                } else {
                  t1.re = reAij / 2.0F;
                  t1.im = temp1 / 2.0F;
                }

                d.re = (t1.re * t1.re - t1.im * t1.im) + (A[ilastm1 + 3 * ilast]
                  .re * A[ilast + 3 * ilastm1].re - A[ilastm1 + 3 * ilast].im *
                  A[ilast + 3 * ilastm1].im);
                d.im = (t1.re * t1.im + t1.im * t1.re) + (A[ilastm1 + 3 * ilast]
                  .re * A[ilast + 3 * ilastm1].im + A[ilastm1 + 3 * ilast].im *
                  A[ilast + 3 * ilastm1].re);
                b_sqrt(&d);
                rho_re = A[ilastm1 + 3 * ilastm1].re - (t1.re - d.re);
                rho_im = A[ilastm1 + 3 * ilastm1].im - (t1.im - d.im);
                anorm = A[ilastm1 + 3 * ilastm1].re - (t1.re + d.re);
                reAij = A[ilastm1 + 3 * ilastm1].im - (t1.im + d.im);
                if (rt_hypotf_snf(rho_re - A[ilast + 3 * ilast].re, rho_im -
                                  A[ilast + 3 * ilast].im) <= rt_hypotf_snf
                    (anorm - A[ilast + 3 * ilast].re, reAij - A[ilast + 3 *
                     ilast].im)) {
                  anorm = rho_re;
                  reAij = rho_im;
                  rho_re = t1.re - d.re;
                  rho_im = t1.im - d.im;
                } else {
                  rho_re = t1.re + d.re;
                  rho_im = t1.im + d.im;
                }
              } else {
                eshift_re += A[ilast + 3 * ilastm1].re;
                eshift_im += A[ilast + 3 * ilastm1].im;
                anorm = eshift_re;
                reAij = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                ctemp.re = A[j + 3 * j].re - anorm;
                ctemp.im = A[j + 3 * j].im - reAij;
                temp1 = imAij * ((real32_T)fabs(ctemp.re) + (real32_T)fabs
                                 (ctemp.im));
                temp2 = imAij * ((real32_T)fabs(A[jp1 + 3 * j].re) + (real32_T)
                                 fabs(A[jp1 + 3 * j].im));
                sumsq = temp1;
                if (temp2 > temp1) {
                  sumsq = temp2;
                }

                if ((sumsq < 1.0F) && (sumsq != 0.0F)) {
                  temp1 /= sumsq;
                  temp2 /= sumsq;
                }

                if (((real32_T)fabs(A[j + 3 * (j - 1)].re) + (real32_T)fabs(A[j
                      + 3 * (j - 1)].im)) * temp2 <= temp1 * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  jp1 = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                if (ifirst == ilastm1 + 1) {
                  ctemp.re = rho_re;
                  ctemp.im = rho_im;
                } else {
                  ctemp.re = A[(ifirst + 3 * (ifirst - 1)) - 1].re - anorm;
                  ctemp.im = A[(ifirst + 3 * (ifirst - 1)) - 1].im - reAij;
                }

                goto90 = true;
              }
            }

            if (goto90) {
              goto90 = false;
              b_eml_matlab_zlartg(ctemp, A[istart + 3 * (istart - 1)], &temp1,
                                  &t1);
              j = istart;
              jp1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  eml_matlab_zlartg(A[(j + 3 * jp1) - 1], A[2 + 3 * jp1], &temp1,
                                    &t1, &A[(j + 3 * jp1) - 1]);
                  A[j + 3 * jp1].re = 0.0F;
                  A[j + 3 * jp1].im = 0.0F;
                }

                for (jp1 = j - 1; jp1 + 1 < 4; jp1++) {
                  d.re = temp1 * A[(j + 3 * jp1) - 1].re + (t1.re * A[j + 3 *
                    jp1].re - t1.im * A[j + 3 * jp1].im);
                  d.im = temp1 * A[(j + 3 * jp1) - 1].im + (t1.re * A[j + 3 *
                    jp1].im + t1.im * A[j + 3 * jp1].re);
                  anorm = A[(j + 3 * jp1) - 1].im;
                  reAij = A[(j + 3 * jp1) - 1].re;
                  A[j + 3 * jp1].re = temp1 * A[j + 3 * jp1].re - (t1.re * A[(j
                    + 3 * jp1) - 1].re + t1.im * A[(j + 3 * jp1) - 1].im);
                  A[j + 3 * jp1].im = temp1 * A[j + 3 * jp1].im - (t1.re * anorm
                    - t1.im * reAij);
                  A[(j + 3 * jp1) - 1] = d;
                }

                t1.re = -t1.re;
                t1.im = -t1.im;
                jp1 = j;
                if (ilast + 1 < j + 2) {
                  jp1 = ilast - 1;
                }

                for (i = 0; i + 1 <= jp1 + 2; i++) {
                  d.re = temp1 * A[i + 3 * j].re + (t1.re * A[i + 3 * (j - 1)].
                    re - t1.im * A[i + 3 * (j - 1)].im);
                  d.im = temp1 * A[i + 3 * j].im + (t1.re * A[i + 3 * (j - 1)].
                    im + t1.im * A[i + 3 * (j - 1)].re);
                  anorm = A[i + 3 * j].im;
                  reAij = A[i + 3 * j].re;
                  A[i + 3 * (j - 1)].re = temp1 * A[i + 3 * (j - 1)].re - (t1.re
                    * A[i + 3 * j].re + t1.im * A[i + 3 * j].im);
                  A[i + 3 * (j - 1)].im = temp1 * A[i + 3 * (j - 1)].im - (t1.re
                    * anorm - t1.im * reAij);
                  A[i + 3 * j] = d;
                }

                for (i = 0; i < 3; i++) {
                  d.re = temp1 * Z[i + 3 * j].re + (t1.re * Z[i + 3 * (j - 1)].
                    re - t1.im * Z[i + 3 * (j - 1)].im);
                  d.im = temp1 * Z[i + 3 * j].im + (t1.re * Z[i + 3 * (j - 1)].
                    im + t1.im * Z[i + 3 * (j - 1)].re);
                  anorm = Z[i + 3 * j].im;
                  reAij = Z[i + 3 * j].re;
                  Z[i + 3 * (j - 1)].re = temp1 * Z[i + 3 * (j - 1)].re - (t1.re
                    * Z[i + 3 * j].re + t1.im * Z[i + 3 * j].im);
                  Z[i + 3 * (j - 1)].im = temp1 * Z[i + 3 * (j - 1)].im - (t1.re
                    * anorm - t1.im * reAij);
                  Z[i + 3 * j] = d;
                }

                jp1 = j - 1;
                j++;
              }
            }

            b_guard1 = true;
          }

          if (b_guard1) {
            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (failed) {
      *info = ilast + 1;
      for (jp1 = 0; jp1 + 1 <= ilast + 1; jp1++) {
        alpha1[jp1].re = ((real32_T)rtNaN);
        alpha1[jp1].im = 0.0F;
        beta1[jp1].re = ((real32_T)rtNaN);
        beta1[jp1].im = 0.0F;
      }

      for (jp1 = 0; jp1 < 9; jp1++) {
        Z[jp1].re = ((real32_T)rtNaN);
        Z[jp1].im = 0.0F;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j + 1 < ilo; j++) {
      alpha1[j] = A[j + 3 * j];
    }
  }
}

/*
 * Arguments    : const creal32_T f
 *                const creal32_T g
 *                float *cs
 *                creal32_T *sn
 *                creal32_T *r
 * Return Type  : void
 */
static void eml_matlab_zlartg(const creal32_T f, const creal32_T g, float *cs,
  creal32_T *sn, creal32_T *r)
{
  float scale;
  float f2s;
  float g2;
  float fs_re;
  float fs_im;
  float gs_re;
  float gs_im;
  int count;
  int rescaledir;
  boolean_T guard1 = false;
  float g2s;
  scale = (real32_T)fabs(f.re);
  f2s = (real32_T)fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  f2s = (real32_T)fabs(g.re);
  g2 = (real32_T)fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = 0;
  rescaledir = 0;
  guard1 = false;
  if (scale >= 5.49755814E+11F) {
    do {
      count++;
      fs_re *= 1.8189894E-12F;
      fs_im *= 1.8189894E-12F;
      gs_re *= 1.8189894E-12F;
      gs_im *= 1.8189894E-12F;
      scale *= 1.8189894E-12F;
    } while (!(scale < 5.49755814E+11F));

    rescaledir = 1;
    guard1 = true;
  } else if (scale <= 1.8189894E-12F) {
    if ((g.re == 0.0F) && (g.im == 0.0F)) {
      *cs = 1.0F;
      sn->re = 0.0F;
      sn->im = 0.0F;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 5.49755814E+11F;
        fs_im *= 5.49755814E+11F;
        gs_re *= 5.49755814E+11F;
        gs_im *= 5.49755814E+11F;
        scale *= 5.49755814E+11F;
      } while (!(scale > 1.8189894E-12F));

      rescaledir = -1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0F > g2) {
      f2s = 1.0F;
    }

    if (scale <= f2s * 1.97215226E-31F) {
      if ((f.re == 0.0F) && (f.im == 0.0F)) {
        *cs = 0.0F;
        r->re = rt_hypotf_snf(g.re, g.im);
        r->im = 0.0F;
        f2s = rt_hypotf_snf(gs_re, gs_im);
        sn->re = gs_re / f2s;
        sn->im = -gs_im / f2s;
      } else {
        g2s = (real32_T)sqrt(g2);
        *cs = rt_hypotf_snf(fs_re, fs_im) / g2s;
        f2s = (real32_T)fabs(f.re);
        g2 = (real32_T)fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0F) {
          f2s = rt_hypotf_snf(f.re, f.im);
          fs_re = f.re / f2s;
          fs_im = f.im / f2s;
        } else {
          g2 = 5.49755814E+11F * f.re;
          scale = 5.49755814E+11F * f.im;
          f2s = rt_hypotf_snf(g2, scale);
          fs_re = g2 / f2s;
          fs_im = scale / f2s;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = (real32_T)sqrt(1.0F + g2 / scale);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0F / f2s;
      f2s = scale + g2;
      g2 = r->re / f2s;
      f2s = r->im / f2s;
      sn->re = g2 * gs_re - f2s * -gs_im;
      sn->im = g2 * -gs_im + f2s * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 1; rescaledir <= count; rescaledir++) {
          r->re *= 5.49755814E+11F;
          r->im *= 5.49755814E+11F;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 1; rescaledir <= count; rescaledir++) {
            r->re *= 1.8189894E-12F;
            r->im *= 1.8189894E-12F;
          }
        }
      }
    }
  }
}

/*
 * Arguments    : const creal32_T A[9]
 *                creal32_T V[9]
 * Return Type  : void
 */
static void eml_matlab_ztgevc(const creal32_T A[9], creal32_T V[9])
{
  float rworka[3];
  int j;
  float anorm;
  int b_j;
  float d_re;
  float xmx;
  float ascale;
  int je;
  float temp;
  float salpha_re;
  float salpha_im;
  float acoeff;
  boolean_T b0;
  boolean_T b1;
  float scale;
  float acoefa;
  creal32_T work1[3];
  int jr;
  float dmin;
  float d_im;
  float work1_re;
  float work1_im;
  creal32_T work2[3];
  for (j = 0; j < 3; j++) {
    rworka[j] = 0.0F;
  }

  anorm = (real32_T)fabs(A[0].re) + (real32_T)fabs(A[0].im);
  for (b_j = 0; b_j < 2; b_j++) {
    for (j = 0; j <= b_j; j++) {
      rworka[b_j + 1] += (real32_T)fabs(A[j + 3 * (b_j + 1)].re) + (real32_T)
        fabs(A[j + 3 * (b_j + 1)].im);
    }

    d_re = rworka[b_j + 1] + ((real32_T)fabs(A[(b_j + 3 * (b_j + 1)) + 1].re) +
      (real32_T)fabs(A[(b_j + 3 * (b_j + 1)) + 1].im));
    if (d_re > anorm) {
      anorm = d_re;
    }
  }

  xmx = anorm;
  if (1.17549435E-38F > anorm) {
    xmx = 1.17549435E-38F;
  }

  ascale = 1.0F / xmx;
  for (je = 0; je < 3; je++) {
    xmx = ((real32_T)fabs(A[(3 * (2 - je) - je) + 2].re) + (real32_T)fabs(A[(3 *
             (2 - je) - je) + 2].im)) * ascale;
    if (1.0F > xmx) {
      xmx = 1.0F;
    }

    temp = 1.0F / xmx;
    salpha_re = ascale * (temp * A[(3 * (2 - je) - je) + 2].re);
    salpha_im = ascale * (temp * A[(3 * (2 - je) - je) + 2].im);
    acoeff = temp * ascale;
    if (((real32_T)fabs(temp) >= 1.17549435E-38F) && ((real32_T)fabs(acoeff) <
         2.95822839E-31F)) {
      b0 = true;
    } else {
      b0 = false;
    }

    if (((real32_T)fabs(salpha_re) + (real32_T)fabs(salpha_im) >=
         1.17549435E-38F) && ((real32_T)fabs(salpha_re) + (real32_T)fabs
         (salpha_im) < 2.95822839E-31F)) {
      b1 = true;
    } else {
      b1 = false;
    }

    scale = 1.0F;
    if (b0) {
      xmx = anorm;
      if (3.3804017E+30F < anorm) {
        xmx = 3.3804017E+30F;
      }

      scale = 2.95822839E-31F / (real32_T)fabs(temp) * xmx;
    }

    if (b1) {
      xmx = 2.95822839E-31F / ((real32_T)fabs(salpha_re) + (real32_T)fabs
        (salpha_im));
      if (xmx > scale) {
        scale = xmx;
      }
    }

    if (b0 || b1) {
      d_re = (real32_T)fabs(acoeff);
      xmx = (real32_T)fabs(salpha_re) + (real32_T)fabs(salpha_im);
      if (1.0F > d_re) {
        d_re = 1.0F;
      }

      if (xmx > d_re) {
        d_re = xmx;
      }

      d_re = 1.0F / (1.17549435E-38F * d_re);
      if (d_re < scale) {
        scale = d_re;
      }

      if (b0) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }

      if (b1) {
        salpha_re *= scale;
        salpha_im *= scale;
      } else {
        salpha_re *= scale;
        salpha_im *= scale;
      }
    }

    acoefa = (real32_T)fabs(acoeff);
    for (jr = 0; jr < 3; jr++) {
      work1[jr].re = 0.0F;
      work1[jr].im = 0.0F;
    }

    work1[2 - je].re = 1.0F;
    work1[2 - je].im = 0.0F;
    xmx = 1.1920929E-7F * ((real32_T)fabs(salpha_re) + (real32_T)fabs(salpha_im));
    dmin = 1.1920929E-7F * acoefa * anorm;
    if (xmx > dmin) {
      dmin = xmx;
    }

    if (1.17549435E-38F > dmin) {
      dmin = 1.17549435E-38F;
    }

    for (jr = 0; jr <= 1 - je; jr++) {
      work1[jr].re = acoeff * A[jr + 3 * (2 - je)].re;
      work1[jr].im = acoeff * A[jr + 3 * (2 - je)].im;
    }

    work1[2 - je].re = 1.0F;
    work1[2 - je].im = 0.0F;
    for (b_j = -1; b_j + 1 <= 1 - je; b_j++) {
      j = -(je + b_j);
      d_re = acoeff * A[j + 3 * j].re - salpha_re;
      d_im = acoeff * A[j + 3 * j].im - salpha_im;
      if ((real32_T)fabs(d_re) + (real32_T)fabs(d_im) <= dmin) {
        d_re = dmin;
        d_im = 0.0F;
      }

      if (((real32_T)fabs(d_re) + (real32_T)fabs(d_im) < 1.0F) && ((real32_T)
           fabs(work1[j].re) + (real32_T)fabs(work1[j].im) >= 2.83568648E+37F *
           ((real32_T)fabs(d_re) + (real32_T)fabs(d_im)))) {
        temp = 1.0F / ((real32_T)fabs(work1[j].re) + (real32_T)fabs(work1[j].im));
        for (jr = 0; jr <= 2 - je; jr++) {
          work1[jr].re *= temp;
          work1[jr].im *= temp;
        }
      }

      work1_re = -work1[j].re;
      work1_im = -work1[j].im;
      if (d_im == 0.0F) {
        if (-work1[j].im == 0.0F) {
          work1[j].re = -work1[j].re / d_re;
          work1[j].im = 0.0F;
        } else if (-work1[j].re == 0.0F) {
          work1[j].re = 0.0F;
          work1[j].im = work1_im / d_re;
        } else {
          work1[j].re = -work1[j].re / d_re;
          work1[j].im = work1_im / d_re;
        }
      } else if (d_re == 0.0F) {
        if (-work1[j].re == 0.0F) {
          work1[j].re = -work1[j].im / d_im;
          work1[j].im = 0.0F;
        } else if (-work1[j].im == 0.0F) {
          work1[j].re = 0.0F;
          work1[j].im = -(work1_re / d_im);
        } else {
          work1[j].re = -work1[j].im / d_im;
          work1[j].im = -(work1_re / d_im);
        }
      } else {
        temp = (real32_T)fabs(d_re);
        xmx = (real32_T)fabs(d_im);
        if (temp > xmx) {
          scale = d_im / d_re;
          xmx = d_re + scale * d_im;
          work1[j].re = (-work1[j].re + scale * -work1[j].im) / xmx;
          work1[j].im = (work1_im - scale * work1_re) / xmx;
        } else if (xmx == temp) {
          if (d_re > 0.0F) {
            scale = 0.5F;
          } else {
            scale = -0.5F;
          }

          if (d_im > 0.0F) {
            xmx = 0.5F;
          } else {
            xmx = -0.5F;
          }

          work1[j].re = (-work1[j].re * scale + -work1[j].im * xmx) / temp;
          work1[j].im = (work1_im * scale - work1_re * xmx) / temp;
        } else {
          scale = d_re / d_im;
          xmx = d_im + scale * d_re;
          work1[j].re = (scale * -work1[j].re + -work1[j].im) / xmx;
          work1[j].im = (scale * work1_im - work1_re) / xmx;
        }
      }

      if (j + 1 > 1) {
        if ((real32_T)fabs(work1[1].re) + (real32_T)fabs(work1[1].im) > 1.0F) {
          temp = 1.0F / ((real32_T)fabs(work1[1].re) + (real32_T)fabs(work1[1].
            im));
          if (acoefa * rworka[1] >= 2.83568648E+37F * temp) {
            for (jr = 0; jr <= 2 - je; jr++) {
              work1[jr].re *= temp;
              work1[jr].im *= temp;
            }
          }
        }

        xmx = acoeff * work1[1].re;
        scale = acoeff * work1[1].im;
        work1[0].re += xmx * A[3].re - scale * A[3].im;
        work1[0].im += xmx * A[3].im + scale * A[3].re;
      }
    }

    for (jr = 0; jr < 3; jr++) {
      work2[jr].re = 0.0F;
      work2[jr].im = 0.0F;
    }

    for (j = 0; j <= 2 - je; j++) {
      for (jr = 0; jr < 3; jr++) {
        work2[jr].re += V[jr + 3 * j].re * work1[j].re - V[jr + 3 * j].im *
          work1[j].im;
        work2[jr].im += V[jr + 3 * j].re * work1[j].im + V[jr + 3 * j].im *
          work1[j].re;
      }
    }

    xmx = (real32_T)fabs(work2[0].re) + (real32_T)fabs(work2[0].im);
    for (jr = 0; jr < 2; jr++) {
      d_re = (real32_T)fabs(work2[jr + 1].re) + (real32_T)fabs(work2[jr + 1].im);
      if (d_re > xmx) {
        xmx = d_re;
      }
    }

    if (xmx > 1.17549435E-38F) {
      temp = 1.0F / xmx;
      for (jr = 0; jr < 3; jr++) {
        V[jr + 3 * (2 - je)].re = temp * work2[jr].re;
        V[jr + 3 * (2 - je)].im = temp * work2[jr].im;
      }
    } else {
      for (jr = 0; jr < 3; jr++) {
        V[jr + 3 * (2 - je)].re = 0.0F;
        V[jr + 3 * (2 - je)].im = 0.0F;
      }
    }
  }
}

/*
 * Arguments    : emxArray_real32_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_real32_T(emxArray_real32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real32_T *)NULL) {
    if (((*pEmxArray)->data != (float *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real32_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_real32_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_real32_T(emxArray_real32_T **pEmxArray, int numDimensions)
{
  emxArray_real32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real32_T *)malloc(sizeof(emxArray_real32_T));
  emxArray = *pEmxArray;
  emxArray->data = (float *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : const float A[81]
 *                float B[9]
 * Return Type  : void
 */
static void mldivide(const float A[81], float B[9])
{
  float b_A[81];
  signed char ipiv[9];
  int i0;
  int j;
  int c;
  int kAcol;
  int ix;
  float temp;
  int k;
  float s;
  int jy;
  int ijA;
  memcpy(&b_A[0], &A[0], 81U * sizeof(float));
  for (i0 = 0; i0 < 9; i0++) {
    ipiv[i0] = (signed char)(1 + i0);
  }

  for (j = 0; j < 8; j++) {
    c = j * 10;
    kAcol = 0;
    ix = c;
    temp = (real32_T)fabs(b_A[c]);
    for (k = 2; k <= 9 - j; k++) {
      ix++;
      s = (real32_T)fabs(b_A[ix]);
      if (s > temp) {
        kAcol = k - 1;
        temp = s;
      }
    }

    if (b_A[c + kAcol] != 0.0F) {
      if (kAcol != 0) {
        ipiv[j] = (signed char)((j + kAcol) + 1);
        ix = j;
        kAcol += j;
        for (k = 0; k < 9; k++) {
          temp = b_A[ix];
          b_A[ix] = b_A[kAcol];
          b_A[kAcol] = temp;
          ix += 9;
          kAcol += 9;
        }
      }

      i0 = (c - j) + 9;
      for (jy = c + 1; jy + 1 <= i0; jy++) {
        b_A[jy] /= b_A[c];
      }
    }

    kAcol = c;
    jy = c + 9;
    for (k = 1; k <= 8 - j; k++) {
      temp = b_A[jy];
      if (b_A[jy] != 0.0F) {
        ix = c + 1;
        i0 = (kAcol - j) + 18;
        for (ijA = 10 + kAcol; ijA + 1 <= i0; ijA++) {
          b_A[ijA] += b_A[ix] * -temp;
          ix++;
        }
      }

      jy += 9;
      kAcol += 9;
    }
  }

  for (kAcol = 0; kAcol < 8; kAcol++) {
    if (ipiv[kAcol] != kAcol + 1) {
      temp = B[kAcol];
      B[kAcol] = B[ipiv[kAcol] - 1];
      B[ipiv[kAcol] - 1] = temp;
    }
  }

  for (k = 0; k < 9; k++) {
    kAcol = 9 * k;
    if (B[k] != 0.0F) {
      for (jy = k + 1; jy + 1 < 10; jy++) {
        B[jy] -= B[k] * b_A[jy + kAcol];
      }
    }
  }

  for (k = 8; k > -1; k += -1) {
    kAcol = 9 * k;
    if (B[k] != 0.0F) {
      B[k] /= b_A[k + kAcol];
      for (jy = 0; jy + 1 <= k; jy++) {
        B[jy] -= B[k] * b_A[jy + kAcol];
      }
    }
  }
}

/*
 * Arguments    : float u0
 *                float u1
 * Return Type  : float
 */
static float rt_hypotf_snf(float u0, float u1)
{
  float y;
  float a;
  float b;
  a = (real32_T)fabs(u0);
  b = (real32_T)fabs(u1);
  if (a < b) {
    a /= b;
    y = b * (real32_T)sqrt(a * a + 1.0F);
  } else if (a > b) {
    b /= a;
    y = a * (real32_T)sqrt(b * b + 1.0F);
  } else if (rtIsNaNF(b)) {
    y = b;
  } else {
    y = a * 1.41421354F;
  }

  return y;
}

/*
 * % calculate DTD in algebraic method
 * Arguments    : const emxArray_real32_T *x
 *                const emxArray_real32_T *y
 *                const emxArray_real32_T *z
 *                float N
 *                float center[3]
 *                float radii[3]
 * Return Type  : void
 */
void ellipsoid_fit_simplify_1023_2(const emxArray_real32_T *x, const
  emxArray_real32_T *y, const emxArray_real32_T *z, float N, float center[3],
  float radii[3])
{
  float DTD[81];
  float DTb[9];
  int r1;
  int r3;
  float b_x[9];
  float X[9];
  int rtemp;
  int r2;
  float A[16];
  float evals_t[9];
  float B[3];
  float maxval;
  float a21;
  double T[16];
  static const double dv0[16] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  float b_T[16];
  float R[16];
  creal32_T b_R[9];
  creal32_T V[9];
  creal32_T beta1[3];
  creal32_T alpha1[3];
  float brm;
  float fv0[3];
  memset(&DTD[0], 0, 81U * sizeof(float));
  for (r1 = 0; r1 < 9; r1++) {
    DTb[r1] = 0.0F;
  }

  for (r3 = 0; r3 < (int)N; r3++) {
    b_x[0] = x->data[(int)(1.0F + (float)r3) - 1] * x->data[(int)(1.0F + (float)
      r3) - 1];
    b_x[1] = y->data[(int)(1.0F + (float)r3) - 1] * y->data[(int)(1.0F + (float)
      r3) - 1];
    b_x[2] = z->data[(int)(1.0F + (float)r3) - 1] * z->data[(int)(1.0F + (float)
      r3) - 1];
    b_x[3] = 2.0F * x->data[(int)(1.0F + (float)r3) - 1] * y->data[(int)(1.0F +
      (float)r3) - 1];
    b_x[4] = 2.0F * x->data[(int)(1.0F + (float)r3) - 1] * z->data[(int)(1.0F +
      (float)r3) - 1];
    b_x[5] = 2.0F * y->data[(int)(1.0F + (float)r3) - 1] * z->data[(int)(1.0F +
      (float)r3) - 1];
    b_x[6] = 2.0F * x->data[(int)(1.0F + (float)r3) - 1];
    b_x[7] = 2.0F * y->data[(int)(1.0F + (float)r3) - 1];
    b_x[8] = 2.0F * z->data[(int)(1.0F + (float)r3) - 1];
    for (rtemp = 0; rtemp < 9; rtemp++) {
      X[rtemp] = b_x[rtemp];
    }

    for (r1 = 0; r1 < 9; r1++) {
      for (r2 = 0; r2 < 9; r2++) {
        DTD[r1 + 9 * r2] += X[r1] * X[r2] / N;
      }

      DTb[r1] += X[r1] / N;
    }
  }

  /*  solve the normal system of equations */
  mldivide(DTD, DTb);
  A[0] = DTb[0];
  A[4] = DTb[3];
  A[8] = DTb[4];
  A[12] = DTb[6];
  A[1] = DTb[3];
  A[5] = DTb[1];
  A[9] = DTb[5];
  A[13] = DTb[7];
  A[2] = DTb[4];
  A[6] = DTb[5];
  A[10] = DTb[2];
  A[14] = DTb[8];
  A[3] = DTb[6];
  A[7] = DTb[7];
  A[11] = DTb[8];
  A[15] = -1.0F;
  for (rtemp = 0; rtemp < 3; rtemp++) {
    for (r1 = 0; r1 < 3; r1++) {
      evals_t[r1 + 3 * rtemp] = -A[r1 + (rtemp << 2)];
    }
  }

  B[0] = DTb[6];
  B[1] = DTb[7];
  B[2] = DTb[8];
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = (real32_T)fabs(evals_t[0]);
  a21 = (real32_T)fabs(evals_t[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if ((real32_T)fabs(evals_t[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  evals_t[r2] /= evals_t[r1];
  evals_t[r3] /= evals_t[r1];
  evals_t[3 + r2] -= evals_t[r2] * evals_t[3 + r1];
  evals_t[3 + r3] -= evals_t[r3] * evals_t[3 + r1];
  evals_t[6 + r2] -= evals_t[r2] * evals_t[6 + r1];
  evals_t[6 + r3] -= evals_t[r3] * evals_t[6 + r1];
  if ((real32_T)fabs(evals_t[3 + r3]) > (real32_T)fabs(evals_t[3 + r2])) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }

  evals_t[3 + r3] /= evals_t[3 + r2];
  evals_t[6 + r3] -= evals_t[3 + r3] * evals_t[6 + r2];
  center[1] = B[r2] - B[r1] * evals_t[r2];
  center[2] = (B[r3] - B[r1] * evals_t[r3]) - center[1] * evals_t[3 + r3];
  center[2] /= evals_t[6 + r3];
  center[0] = B[r1] - center[2] * evals_t[6 + r1];
  center[1] -= center[2] * evals_t[6 + r2];
  center[1] /= evals_t[3 + r2];
  center[0] -= center[1] * evals_t[3 + r1];
  center[0] /= evals_t[r1];
  memcpy(&T[0], &dv0[0], sizeof(double) << 4);
  for (rtemp = 0; rtemp < 3; rtemp++) {
    T[3 + (rtemp << 2)] = center[rtemp];
  }

  for (rtemp = 0; rtemp < 4; rtemp++) {
    for (r1 = 0; r1 < 4; r1++) {
      b_T[rtemp + (r1 << 2)] = 0.0F;
      for (r2 = 0; r2 < 4; r2++) {
        b_T[rtemp + (r1 << 2)] += (float)T[rtemp + (r2 << 2)] * A[r2 + (r1 << 2)];
      }
    }

    for (r1 = 0; r1 < 4; r1++) {
      R[rtemp + (r1 << 2)] = 0.0F;
      for (r2 = 0; r2 < 4; r2++) {
        R[rtemp + (r1 << 2)] += b_T[rtemp + (r2 << 2)] * (float)T[r1 + (r2 << 2)];
      }
    }
  }

  for (rtemp = 0; rtemp < 3; rtemp++) {
    for (r1 = 0; r1 < 3; r1++) {
      b_R[r1 + 3 * rtemp].re = R[r1 + (rtemp << 2)] / -R[15];
      b_R[r1 + 3 * rtemp].im = 0.0F;
    }
  }

  eml_matlab_zggev(b_R, &r1, alpha1, beta1, V);
  for (rtemp = 0; rtemp < 9; rtemp++) {
    V[rtemp].re = 0.0F;
    V[rtemp].im = 0.0F;
  }

  for (r1 = 0; r1 < 3; r1++) {
    if (beta1[r1].im == 0.0F) {
      if (alpha1[r1].im == 0.0F) {
        V[r1 + 3 * r1].re = alpha1[r1].re / beta1[r1].re;
        V[r1 + 3 * r1].im = 0.0F;
      } else if (alpha1[r1].re == 0.0F) {
        V[r1 + 3 * r1].re = 0.0F;
        V[r1 + 3 * r1].im = alpha1[r1].im / beta1[r1].re;
      } else {
        V[r1 + 3 * r1].re = alpha1[r1].re / beta1[r1].re;
        V[r1 + 3 * r1].im = alpha1[r1].im / beta1[r1].re;
      }
    } else if (beta1[r1].re == 0.0F) {
      if (alpha1[r1].re == 0.0F) {
        V[r1 + 3 * r1].re = alpha1[r1].im / beta1[r1].im;
        V[r1 + 3 * r1].im = 0.0F;
      } else if (alpha1[r1].im == 0.0F) {
        V[r1 + 3 * r1].re = 0.0F;
        V[r1 + 3 * r1].im = -(alpha1[r1].re / beta1[r1].im);
      } else {
        V[r1 + 3 * r1].re = alpha1[r1].im / beta1[r1].im;
        V[r1 + 3 * r1].im = -(alpha1[r1].re / beta1[r1].im);
      }
    } else {
      brm = (real32_T)fabs(beta1[r1].re);
      maxval = (real32_T)fabs(beta1[r1].im);
      if (brm > maxval) {
        maxval = beta1[r1].im / beta1[r1].re;
        a21 = beta1[r1].re + maxval * beta1[r1].im;
        V[r1 + 3 * r1].re = (alpha1[r1].re + maxval * alpha1[r1].im) / a21;
        V[r1 + 3 * r1].im = (alpha1[r1].im - maxval * alpha1[r1].re) / a21;
      } else if (maxval == brm) {
        if (beta1[r1].re > 0.0F) {
          maxval = 0.5F;
        } else {
          maxval = -0.5F;
        }

        if (beta1[r1].im > 0.0F) {
          a21 = 0.5F;
        } else {
          a21 = -0.5F;
        }

        V[r1 + 3 * r1].re = (alpha1[r1].re * maxval + alpha1[r1].im * a21) / brm;
        V[r1 + 3 * r1].im = (alpha1[r1].im * maxval - alpha1[r1].re * a21) / brm;
      } else {
        maxval = beta1[r1].re / beta1[r1].im;
        a21 = beta1[r1].im + maxval * beta1[r1].re;
        V[r1 + 3 * r1].re = (maxval * alpha1[r1].re + alpha1[r1].im) / a21;
        V[r1 + 3 * r1].im = (maxval * alpha1[r1].im - alpha1[r1].re) / a21;
      }
    }
  }

  for (rtemp = 0; rtemp < 9; rtemp++) {
    evals_t[rtemp] = V[rtemp].re;
  }

  fv0[0] = 1.0F / evals_t[0];
  fv0[1] = 1.0F / evals_t[4];
  fv0[2] = 1.0F / evals_t[8];
  for (rtemp = 0; rtemp < 3; rtemp++) {
    radii[rtemp] = fv0[rtemp];
  }

  /*      radii = single([1,1,1]'); */
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void ellipsoid_fit_simplify_1023_2_initialize(void)
{
  // rt_InitInfAndNaN(8U);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void ellipsoid_fit_simplify_1023_2_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * Arguments    : int numDimensions
 *                int *size
 * Return Type  : emxArray_real32_T *
 */
emxArray_real32_T *emxCreateND_real32_T(int numDimensions, int *size)
{
  emxArray_real32_T *emx;
  int numEl;
  int i;
  emxInit_real32_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (float *)calloc((unsigned int)numEl, sizeof(float));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

/*
 * Arguments    : float *data
 *                int numDimensions
 *                int *size
 * Return Type  : emxArray_real32_T *
 */
emxArray_real32_T *emxCreateWrapperND_real32_T(float *data, int numDimensions,
  int *size)
{
  emxArray_real32_T *emx;
  int numEl;
  int i;
  emxInit_real32_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

/*
 * Arguments    : float *data
 *                int rows
 *                int cols
 * Return Type  : emxArray_real32_T *
 */
emxArray_real32_T *emxCreateWrapper_real32_T(float *data, int rows, int cols)
{
  emxArray_real32_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real32_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

/*
 * Arguments    : int rows
 *                int cols
 * Return Type  : emxArray_real32_T *
 */
emxArray_real32_T *emxCreate_real32_T(int rows, int cols)
{
  emxArray_real32_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real32_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (float *)calloc((unsigned int)numEl, sizeof(float));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

/*
 * Arguments    : emxArray_real32_T *emxArray
 * Return Type  : void
 */
void emxDestroyArray_real32_T(emxArray_real32_T *emxArray)
{
  emxFree_real32_T(&emxArray);
}

/*
 * Arguments    : emxArray_real32_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInitArray_real32_T(emxArray_real32_T **pEmxArray, int numDimensions)
{
  emxInit_real32_T(pEmxArray, numDimensions);
}

/*
 * File trailer for ellipsoid_fit_simplify_1023_2.c
 *
 * [EOF]
 */
