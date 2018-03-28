/* Document here, INSIDE A COMMENT, any other students you
   collaborated with (other than one you partnered with for this
   assignment), and anything else (besides FSV) that helped you. */
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (0L in ref)

#ifndef VCR_HAS_BEEN_INCLUDED
#define VCR_HAS_BEEN_INCLUDED

#include <teem/meet.h>

#ifdef __cplusplus
extern "C" {
#endif

/* "real" is a typedef for either float or double; same as last time.
   You can compile with:

   make CFLAGS="-DVCR_REAL_IS_DOUBLE=1"

   to make "real" mean "double"; otherwise (with regular "make") "real"
   means "float". Your code must compile and work with real either as
   float or double, and give more accurate results with real as double.
   Make sure to to "make clean" when changing the accuracy of real
   (the Makefile doesn't know which files have been compiled how). */
#ifndef VCR_REAL_IS_DOUBLE
#define VCR_REAL_IS_DOUBLE 0
#endif
#if VCR_REAL_IS_DOUBLE
typedef double real;
#else
typedef float real;
#endif

/* this is non-standard; but "uint" is easier to type than any of
   C99's uint8_t, uint16_t, uint32_t, and there is no need here for
   any specific size; just an "unsigned int" */
typedef unsigned int uint;

/*
** The different kinds of integration we support; the enum values
** here are packaged into airEnum vcrIntg_ae
*/
typedef enum {
  vcrIntgUnknown,   // 0: don't know */
  vcrIntgEuler,     // 1: Euler integration
  vcrIntgMidpoint,  // 2: Midpoint method aka RK2
  vcrIntgRK4        // 3: Runge-Kutta fourth-order
} vcrIntg;
#define VCR_INTG_MAX   3

/*
** different (scalar-valued) properties of a vector field, which can be
** measured at each pixel during LIC computation (for the seed point).
** Note that we're using "Jacobian" for the matrix that others may call
** the "Jacobian matrix".
*/
typedef enum {
  vcrPropUnknown,  // 0: don't know
  vcrPropOne,      // 1: always 1.0
  vcrPropLength,   // 2: length of vector
  vcrPropVort,     // 3: field vorticity (from Jacobian)
  vcrPropDiverg,   // 4: field divergence (from Jacobian)
  vcrPropJacFrob,  // 5: Frobenius norm of field Jacobian
  vcrPropJacDet    // 6: determinant of field Jacobian
} vcrProp;
#define VCR_PROP_MAX  6

/*
** vcrType: the scalar pixel types supported in vcrImage
*/
typedef enum {
  vcrTypeUnknown,  // (0) no type known
  vcrTypeUChar,    // (1) unsigned char
  vcrTypeFloat,    // (2) float
} vcrType;
#define VCR_TYPE_MAX   2

/*
** The vcrKernel stores everything about a reconstruction kernel. The kernel
** is non-zero only within [-support/2,support/2], for integer "support" which
** may be odd or even (but always positive). The kernels are set up at
** compile-time in such a way that each kernel knows its own derivative
** "deriv", so the derivative of kernel k is k->deriv. Implemented in kernel.c
** and kparse.c.
*/
typedef struct vcrKernel_t {
  const char *name;                // short identifying string
  const char *desc;                // short descriptive string
  uint support;                    // # samples needed for convolution
  real (*eval)(real xx);           // how to evaluate the kernel
  const struct vcrKernel_t *deriv; // the derivative of this kernel
} vcrKernel;

/*
** struct vcrImage: like project 2's mprImage. Implemented in image.c
*/
typedef struct vcrImage_t {
  uint channel;   // 1 for scalar, 2 for vector, 3 for color
  uint size[2];   /* # samples on faster (size[0]), slower (size[1]) spatial
                     axis; these are the second and third axes in the
                     linearization; the fastest axis is one holding the two
                     vector components, or "channel" image components */
  real ItoW[9];   /* homogeneous coordinate mapping from index-space (faster
                     coordinate first) to the "right-up" world-space */
  vcrType dtype;  /* type of the data; determines which of the union members
                     below to use */
  union {         /* union for the pointer to the image data; the pointer
                     values are all the same; this is just to avoid casting.
                     The right union member to use (data.uc vs data.fl)
                     determined at run-time by value of dtype */
    void *vd;
    unsigned char *uc;
    float *fl;
  } data;
} vcrImage;

/*
** struct vcrCtx: all the state associated with computing convolution,
** and streamlines, and LIC (line integral convolution). Mostly
** implemented in ctx.c
*/
typedef struct vcrCtx_t {
  const vcrImage *vfl;     /* vector field to process. Because of the checks
                              at the start of vcrCtxNew, the only type of data
                              in which you have to implement convolution is
                              2-vectors of floats. The data will be linearized
                              with the 2 vector coordinates on the fastest
                              axis, then the vfl->size[0] (faster) and
                              vfl->size[1] (slower) axes sampling the domain
                              of the vector field */
  const vcrKernel *kern;   // kernel to reconstruct with

  /* put here the same sort of extra information (that depends on the
     vector field and kernel given to vcrCtxNew) that you put
     inside the mprContext from Project 2 */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real* WtoI;
  int odd_support;
  int upper;
  int lower;
  real* deriv_ItoW;
  real* kvalues1;
  real* kvalues2;
  real* the_kvalues1;
  real* the_kvalues2;
  const vcrKernel* the_kern;
  real* dpos;
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (13L in ref)

  /* output to be set (by you) in vcrConvoEval:
     inside: the following convolution results could be computed
     vec: the vector reconstructed componentwise by convolution
     jac: Jacobian matrix of vector field, linearized as:

       jac[0]  jac[1]        d vec[0]/dx   d vec[0]/dy
                        ==
       jac[2]  jac[3]        d vec[1]/dx   d vec[1]/dy

     Note especially the meaning of the off-diagonal entries.
  */
  int inside;
  real vec[2], jac[4];
} vcrCtx;

/*
** vcrSline: A container for a single streamline. vcrSlineTrace uses this to
** store its results. "halfLen" determines for how many positions the "pos"
** array is allocated (as detailed below), and all the remaining fields are
** set by vcrSlineTrace.  Implemented in sline.c.
*/
typedef struct vcrSline_t {
  uint halfLen;     /* pos is allocated for 2*(1 + 2*halfLen) reals, i.e.
                       for (1 + 2*halfLen) coordinate 2-vectors.
                       pos+2*halfLen is the position of the starting (seed)
                       point of the streamline.  Even if halfLen is 0, pos
                       stores a position, in which case this struct is just a
                       way of storing a vector at some location */
  real *pos;        /* world-space coords of all points along streamline, as
                       a 2 (fast) by 2*(1+2*halfLen) (slow) array.  This is
                       NULL upon return from vcrSlineNew() */
  int seedInside;   /* the seed point was indeed inside the field, so the
                       streamline integration could proceed.  At least the
                       2-vector at pos+2*halfLen is set (to the location of
                       the seed point), and more positions are set according
                       to how long the integration could proceed forward and
                       backward, as determined by vcrSlineTrace(), and as
                       recorded in forwNum and backNum. */
  uint forwNum,     /* the forward (downstream, following the vectors) part
                       of streamline is from pos+2*halfLen (the seed point)
                       to pos+2*(halfLen+forwNum) (the head of the
                       streamline, considered as one arrow) */
    backNum;        /* the backward (upstream, following the negation of the
                       vectors) part of streamline is from
                       pos+2*(halfLen-backNum) (the tail of the streamline,
                       considered as one arrow) to pos+2*halfLen (the
                       seedpoint) */
  real vecSeed[2],  // vector (from vector field) at seed point
    jacSeed[4];     // field jacobian at seed point
} vcrSline;

// misc.c: miscellaneous little things (nothing for you to do)
#ifndef NDEBUG
extern int vcrVerbose;
#else
#define vcrVerbose 0
#endif
extern const int vcrRealIsDouble;
extern const char *vcrBiffKey;
#define VCR vcrBiffKey
extern real vcrNan(unsigned short payload);
extern unsigned short vcrNanPayload(real nval);
extern const airEnum *const vcrIntg_ae;
extern const airEnum *const vcrProp_ae;
extern const airEnum *const vcrType_ae;
extern uint vecSprintLenMax(uint oldmax, const real *vv, uint vvNum);
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (0L in ref)

/* kernel.c: gives compile-time definition of some reconstruction kernels
   and their derivatives (nothing for you to do) */
extern const vcrKernel *const vcrKernelZero;
extern const vcrKernel *const vcrKernelBox;
extern const vcrKernel *const vcrKernelCos;
extern const vcrKernel *const vcrKernelSin;
extern const vcrKernel *const vcrKernelTent;
extern const vcrKernel *const vcrKernelTooth;
extern const vcrKernel *const vcrKernelBspln2;
extern const vcrKernel *const vcrKernelBspln3;
extern const vcrKernel *const vcrKernelCtmr;
extern const vcrKernel *const vcrKernelC4hexic;
extern const vcrKernel *const vcrKernelAll[];

// kparse.c: for identifying kernels (nothing for you to do)
extern const vcrKernel *vcrKernelParse(const char *kstr);
extern hestCB *vcrKernelHest;

// image.c: for representing vector data (nothing for you to do)
extern vcrImage *vcrImageNew(void);
extern vcrImage *vcrImageNix(vcrImage *vol);
extern int vcrImageAlloc(vcrImage *img, uint channel,
                         uint size0, uint size1,
                         vcrType dtype);
extern int vcrImageLoad(vcrImage *img, const char *fname);
extern int vcrImageSave(const char *fname, const vcrImage *img);

/* ctx.c: for setting up and using the vcrCtx. You can augment
   vcrCtxNew() and vcrCtxNix(), as well as the vcrCtx struct
   itself (above) */
extern vcrCtx *vcrCtxNew(const vcrImage *vec, const vcrKernel *kern);
extern vcrCtx *vcrCtxNix(vcrCtx *ctx);

// convo.c: for convolution
extern void vcrConvoEval(vcrCtx *ctx, real xw, real yw, int getjac);

// sline.c: for storing and computing streamlines
extern vcrSline *vcrSlineNew(void);
extern int vcrSlineAlloc(vcrSline *sln, uint halfLen);
extern vcrSline *vcrSlineNix(vcrSline *sln);
extern int vcrSlineTrace(vcrSline *sln, real seedX, real seedY,
                         uint halfLen, real hh,
                         int normalize, int intg, vcrCtx *ctx);

// lic.c: for Line Integral Convolution, as well as "property" calculation
extern real vcrPropCalc(int prop, const real vec[2], const real jac[4]);
extern int vcrLICEval(real *result, vcrSline *sln, real wx, real wy,
                      const vcrKernel *wkern, const vcrImage *rnd, int rndLinterp,
                      uint halfLen, real hh, int normalize, int intg,
                      vcrCtx *ctx);
extern int vcrLIC(vcrImage *lmg, vcrImage *pmg, int prop,
                  const vcrKernel *wkern, const vcrImage *rnd, int rndLinterp,
                  uint halfLen, real hh, int normalize, int intg,
                  vcrCtx *ctx);

#ifdef __cplusplus
}
#endif
#endif // VCR_HAS_BEEN_INCLUDED
