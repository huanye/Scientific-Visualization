/*
** mpr and mapr: SciVis-2018 Project 2
** Copyright (C)  2018 University of Chicago. All rights reserved.
**
** This file is distributed to students in the 2018 CMSC 23710/33710
** ("SciVis") class, and no one else, for use in the SciVis class. This
** file is not open-source licensed, nor licensed for any other distribution.
** You should not allow this file to be copied or downloaded by someone who
** is not a 2018 SciVis student, as that would violate the copyright
** held by University of Chicago.
*/

/* Document here, INSIDE A COMMENT, any other students you
   collaborated with (other than one you partnered with for this
   assignment), and anything else (besides FSV) that helped you. */
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
// ( 0 lines in reference implementation )
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

#ifndef MPR_HAS_BEEN_INCLUDED
#define MPR_HAS_BEEN_INCLUDED

#include <teem/meet.h>

#ifdef __cplusplus
extern "C" {
#endif

/* This header uses "real" as a typedef for either float (single precision)
   or double (double precision): you choose which at compile time, so that,
   for example, you can more easily experiment with the performance/accuracy
   trade-off. This is enabled by the above inclusion in mprPrivate.h of
   <tgmath.h>: it is a huge hack (see http://goo.gl/dEQE4R), but it allows
   "cos(x)" to become cosf(x) for float arguments and regular cos(x) for
   double arguments. You can compile with:

   make CFLAGS="-DMPR_REAL_IS_DOUBLE=1"

   to make "real" mean "double"; otherwise (with regular "make") "real"
   means "float". */
#ifndef MPR_REAL_IS_DOUBLE
#define MPR_REAL_IS_DOUBLE 0
#endif
#if MPR_REAL_IS_DOUBLE
typedef double real;
#else
typedef float real;
#endif

/* this is non-standard; but "uint" is easier to type than any of
   C99's uint8_t, uint16_t, uint32_t, and there is no need here for
   any specific size; just an "unsigned int" */
typedef unsigned int uint;

/*
** mprType: the scalar pixel types supported in mprImage
** mprTypeUChar is only used for output images.
** All input images (and some outputs) are mprTypeFloat.
** The single-precision float is used even with mprRealIsDouble.
*/
typedef enum {
  mprTypeUnknown,   // (0) (no type known)
  mprTypeUChar,     // (1) 1-byte unsigned char
  mprTypeFloat,     // (2) 4-byte float
} mprType;

/*
** mprColor: the colorspaces supported as the coordinates of control
** points in colormaps. In all of them, all components vary in [0,1]
*/
typedef enum {
  mprColorUnknown,  // (0) (no space known)
  mprColorRGB,      // (1) red-green-blue
  mprColorHSV,      // (2) hue-saturation-value single hexcone
  mprColorLast,     // (3) (not a valid color)
} mprColor;

/*
** mprMode: the different ways that mprPictureSample() can operate.  Each mode
** is described here by what kind of values are set in the output image by
** mprPictureEval(), and noting if it needs the gradient to be reconstructed.
** Whether the gradient is needed is also returned by mprModeNeedsGradient().
** These enum values are included in airEnum mprMode_ae in misc.c.
*/
typedef enum {
  mprModeUnknown,         // (0) mode not known
  mprModeWorldPos,        /* (1) 2-vec (real): for debugging: the world-
                             space position being sampled
                             (does not need gradient) */
  mprModeIndexPos,        /* (2) 2-vec (real): for debugging: index-
                             space position in the input image being sampled
                             (does not need gradient) */
  mprModeOutside,         /* (3) scalar (real): for debugging: how many
                             samples on faster and slower axis were missing
                             in convolution sum
                             (does not need gradient) */
  mprModeValue,           /* (4) scalar (real): reconstructeded value
                             (does not need gradient) */
  mprModeGradient,        /* (5) 2-vecs (real): reconstructed gradient;
                             DOES NEED gradient */
  mprModeFuzzyIso,        /* (6) scalar (real): highest when value is
                             ctx->fuzzyIsoVal, with thickness roughly
                             ctx->fuzzyIsoThick;
                             DOES NEED gradient */
  mprModeValueCmap,       /* (7) rgb color (uchar): colormapped value
                             (does not need gradient) */
  mprModeValueCmapShaded, /* (8) rgb color (uchar): colormapped value,
                             shaded according to gradient direction
                             relative to some synthetic light source,
                             DOES NEED gradient */
} mprMode;
#define MPR_MODE_MAX 8

/*
** The mprKernel stores everything about a reconstruction kernel. The kernel
** is non-zero only within [-support/2,support/2], for integer "support"
** which may be odd or even (but always positive). The kernels are set up at
** compile-time in such a way that each kernel knows its own derivative
** "deriv", so the derivative of kernel k is k->deriv. Implemented in
** kernel.c and kparse.c.
*/
typedef struct mprKernel_t {
  const char *name;                // short identifying string
  const char *desc;                // short descriptive string
  uint support;                    // # samples needed for convolution
  real (*eval)(real xx);           // how to evaluate the kernel
  const struct mprKernel_t *deriv; /* derivative of this kernel; will point
                                      back to itself when kernel is zero */
} mprKernel;

/*
** The mprCmap stores a univariate colormap, consisting of "num" control
** points, and a single "outside" color, all specified in colorspace "space".
** The control points are in the "data" array, which is logically a 4-by-num
** (faster-by-slower) array of reals, or a length-"num" array of (X,A,B,C)
** 4-vectors, where X is location of control point, and (A,B,C) are color
** coordinates.  The X values have to strictly increase with control point
** index. Implemented in cmap.c.
*/
typedef struct {
  mprColor space;       /* the color space in which of colormap control
                           points (and outside[]) have their coordinates */
  real outside[3];      // color of not having a convolution result
  real *data;           // control point data
  uint num;             // number of color control points
  /* These are set via mprCmapLutGen to generate a faster? look-up table
     ("lut") version of colormap, only as RGB. Some of this is redundant:
     lutMin is data[0] and lutMax is data[0 + 3*(num-1)] */
  real lutMin, lutMax;  // range of the cell-centered lut samples
  real *lut;            // 3*lutLen array of R,G,B values
  uint lutLen;          // logical length of lut
} mprCmap;

/*
** mprImage is a container for the oriented image data that is processed by
** mapr. Image orientation is defined in terms of a particular world-space
** basis: an orthonormal {right, up}. Implemented in image.c
*/
typedef struct {
  uint channel,        /* how many values are at each pixel; this always
                          fastest axis */
    size[2];           /* # of samples along faster (size[0]) and
                          slower (size[1]) spatial axes */
  real ItoW[9];        /* homogeneous coordinate mapping from index-space
                          (faster coordinate first) to the "right-up"
                          world-space */
  mprType dtype;       /* type of the data; determines which of the union
                          members below to use */
  union {              /* union for the pointer to the image data; the
                          pointer values are all the same; this is just to
                          avoid casting.  The right union member to use
                          (data.uc vs data.fl) determined at run-time by
                          value of dtype */
    void *vd;
    unsigned char *uc;
    float *fl;         /* Floating-point image data doesn't care whether
                          we've typedef'd "real" as double or float; it
                          will always come as single-precision float */
  } data;
} mprImage;

/*
** mprArray is a container for a dynamically-resized array, to be used with
** mprIsocontour to manage the list of line segments comprising the
** isocontour. Implemented in misc.c.
*/
typedef struct {
  void *data;   // where is the data
  uint unit,    // size in bytes of one element in the array
    len,        /* nominal length of array == # units for which there is
                   considered to be data <= total # units allocated. */
    size,       // # units allocated
    initial;    // initial and minimal size
} mprArray;

/*
** mprCtx is a container for all the state associated with doing convolution,
** colormapping, and shading. Implemented in ctx.c, but used in other files.
** TODO: You add to this struct definition as noted in the comment below.
*/
typedef struct {
  // input
  int verbose;
  const mprImage *image;
  const mprKernel *kern;
  mprMode mode;
  real fuzzyIsoVal, fuzzyIsoThick; // for mprModeFuzzyIso
  const mprCmap *cmap; // for mprModeValueCmap, mprModeValueCmapShaded
  real imgMinMax[2];  // for mprModeValueCmap, mprModeValueCmapShaded
  real ldir[3], shading, zscl; // for mprModeValueCmapShaded

  /* The following section is for the declaration of any information or
     buffers that you want to use for mprConvoEval() and mprPictureEval().
     Arrays here should be allocated by mprCtxNew(), and later freed by
     mprCtxNix().  To figure out what information should be put here, answer
     the question: "Of the information or buffers needed to to complete
     mprConvoEval() and mprPictureEval(), what can be computed or allocated
     once in mprCtxNew(), ahead of time, based on knowing the "image",
     "kern", "mode", and "cmap?"  You want to minimize needlessly repeating
     the same computation or allocation with each call to mprConvoEval() or
     mprPictureEval(). */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real* WtoI;
  int odd_support;
  int upper;
  int lower;
  real* deriv_ItoW;
  // ( 14 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

  // output fields set by mprConvoEval
  real wpos[2],  // copy of world-space pos passed to mprConvoEval
    ipos[2];     /* wpos converted to index-space, via inverse of
                    ctx->image->ItoW; note that this is a 2-vector
                    not a 3-vector of homogeneous coords. */
  /* If the requested convolution location puts the kernel support
     outside the valid image index domain, "outside" records the
     number of indices missing on the fast axis, plus the number of
     indices missing on the slow axes. */
  uint outside;
  real value;       // if !outside, convolution result
  real gradient[2]; /* if mprModeNeedsGradient(mode) is non-zero
                        and !outside, the *world-space* gradient */

  /* output set by mprPictureSample: the rate, in kHz, at which output
     pixels could be computed (includes convolution and possibly
     colormapping) */
  real srate;
} mprCtx;

/* misc.c: miscellaneous things (nothing for you to implement)
   mprVerbose is like pltVerbose, read p1plotr/plt.h for more */
#ifndef NDEBUG
extern int mprVerbose;
#else
#define mprVerbose 0
#endif
extern const int mprRealIsDouble;
extern const char *mprBiffKey;
#define MPR mprBiffKey  // identifies library in error messages
extern real mprNan(unsigned short payload);
extern unsigned short mprNanPayload(real nval);
extern const airEnum *const mprType_ae;
extern const airEnum *const mprColor_ae;
extern const airEnum *const mprMode_ae;
extern int mprModeNeedsGradient(mprMode mode);
extern mprArray *mprArrayNew(uint unit, uint initial);
extern mprArray *mprArrayNix(mprArray *arr);
extern uint mprArrayLenIncr(mprArray *arr, uint incr);
extern uint mprSprintLenMax(uint oldmax, real *vv, uint vvNum);

/* util.c: utility stuff. TODO: finish mprQuantize */
extern uint mprQuantize(real min, real val, real max, uint num);

/* kernel.c: gives compile-time definition of some reconstruction
   kernels and their derivatives (nothing for you to implement) */
extern const mprKernel *const mprKernelZero;
extern const mprKernel *const mprKernelBox;
extern const mprKernel *const mprKernelTent;
extern const mprKernel *const mprKernelBspln2;
extern const mprKernel *const mprKernelBspln3;
extern const mprKernel *const mprKernelCtmr;
extern const mprKernel *const mprKernelSpark;
extern const mprKernel *const mprKernelLuna;
extern const mprKernel *const mprKernelCelie;
extern const mprKernel *const mprKernelAll[];

/* kparse.c: for identifying kernels by strings (nothing for you to
   implement) */
extern const mprKernel *mprKernelParse(const char *kstr);
extern hestCB *mprKernelHest;

/* image.c: for working with mprImage struct (nothing for you to
   implement) */
extern mprImage *mprImageNew(void);
extern mprImage *mprImageNix(mprImage *img);
extern int mprImageAlloc(mprImage *img, uint channel,
                         uint size0, uint size1,
                         mprType type);
extern int mprImageNrrdWrap(Nrrd *nout, const mprImage *img);
extern int mprImageSave(const char *fname, const mprImage *img, real srate);
extern int mprImageLoad(mprImage *img, const char *fname);
extern int mprImageMinMax(real minmax[2], const mprImage *img);

/* cmap.c: for working with mprCamp struct.
   TODO: implement mprCmapEval */
extern mprCmap *mprCmapNew(void);
extern mprCmap *mprCmapNix(mprCmap *cmap);
extern int mprCmapAlloc(mprCmap *cmap, mprColor space,
                        const real outside[3],
                        uint pointNum);
extern int mprCmapLoad(mprCmap *cmap, const char *fname);
extern void mprCmapEval(real *rgb, const mprCmap *cmap, real val);
extern int mprCmapLutGen(mprCmap *cmap, uint len);
extern int mprCmapLutDraw(mprImage *img, const mprCmap *cmap,
                          uint height);

/* ctx.c: for working with the mprCtx struct.  TODO: finish mprCtxNew
   and mprCtxNix to allocate and free (respectively) anything
   dynamically-allocated that you add to the mprCtx struct */
extern mprCtx *mprCtxNew(const mprImage *img, const mprKernel *kern,
                         mprMode mode, const mprCmap *cmap,
                         const real *imm);
extern mprCtx *mprCtxNix(mprCtx *ctx);
extern int mprCtxModeFuzzyIsoParmSet(mprCtx *ctx,
                                     real isoval, real thickness);
extern int mprCtxModeValueCmapShadedParmSet(mprCtx *ctx,
                                            const real ldir[3], real shading,
                                            real zscl);

/* convo.c: for doing convolution. TODO: finish mprConvoEval */
extern void mprConvoEval(mprCtx *ctx, real xw, real yw);

/* picture.c: for creating a picture of the repeated convolution
   results.  TODO: finish mprPictureRGB */
extern void mprHSVfromRGB(real *hsv, const real *rgb);
extern void mprRGBfromHSV(real *rgb, const real *hsv);
extern void mprToRGB(real *rgb, const real *col, mprColor space);
extern int mprPictureItoW(real ItoW[9],
                          uint size0, uint size1,
                          real centX, real centY, real fov,
                          real scl0, real scl1,
                          real angle, real shear0, real shear1);
extern void mprPictureRGB(unsigned char *rgb, const mprCtx *ctx);
extern int mprPictureSample(mprImage *omg, mprCtx *ctx,
                            uint size0, uint size1,
                            real centX, real centY, real fov,
                            real angle);

/* march.c: computing Marching Squares isosurfaces.
   TODO: finish mprIsocontour */
extern int mprIsocontour(mprArray *segArr, const mprImage *img,
                         real isoval, int skipOutOfView,
                         uint size0, uint size1,
                         real centX, real centY, real fov,
                         real angle);

/*
** The "mapr" executable is compiled from mapr.c, which relies on
** the various mpr_XCmd objects (each defined in mapr_X.c) for the
** mapr commands X={about, klist, ...}  None of those objects need
** declaring here, because mapr is the only thing using them, and
** mapr.c generates (via C pre-processor tricks) its own extern
** declarations of them.
**
** You should run "rmapr" to see all the commands, then run "rmapr
** about" to learn more about what work is required to finish this
** project. For mapr command X, you can look at mapr_X.c to see how
** it works, and to see the context in which the mpr library
** functions above are called.
*/

#ifdef __cplusplus
}
#endif
#endif // MPR_HAS_BEEN_INCLUDED
