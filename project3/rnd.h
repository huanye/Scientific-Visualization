/* Document here, INSIDE A COMMENT, any other students you
   collaborated with (other than one you partnered with for this
   assignment), and anything else (besides FSV) that helped you. */
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
// ( 0 lines in reference implementation )
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

#ifndef RND_HAS_BEEN_INCLUDED
#define RND_HAS_BEEN_INCLUDED

#include <teem/meet.h>
#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

/* This header uses "real" as a typedef for either float (single precision)
   or double (double precision): you choose which at compile time, so that,
   for example, you can more easily experiment with the performance/accuracy
   trade-off. This is enabled by the above inclusion in rndPrivate.h of
   <tgmath.h>: it is a huge hack (see http://goo.gl/dEQE4R), but it allows
   "cos(x)" to become cosf(x) for float arguments and regular cos(x) for
   double arguments. You can compile with:

   make CFLAGS="-DRND_REAL_IS_DOUBLE=1"

   to make "real" mean "double"; otherwise (with regular "make") "real"
   means "float". Your code must compile and work with real either as
   float or double, and give more accurate results with real as double.
   Make sure to to "make clean" when changing the accuracy of real
   (the Makefile doesn't know which files have been compiled how). */
#ifndef RND_REAL_IS_DOUBLE
#define RND_REAL_IS_DOUBLE 0
#endif
#if RND_REAL_IS_DOUBLE
typedef double real;
#else
typedef float real;
#endif

/* this is non-standard; but "uint" is easier to type than any of
   C99's uint8_t, uint16_t, uint32_t, and there is no need here for
   any specific size; just an "unsigned int" */
typedef unsigned int uint;

/*
** enum rndType: the possible voxel data types supported in a rndVolume
** (whereas in the previous project, the only type on which convolution might
** happen was float).  Real-world scanned volume data (from CT or MRI
** scanners) is commonly stored in 16-bit shorts (even if not all 16 bits are
** actually used), and simulations or other experiments may generate float
** values. Floating point volume data is stored as single-precision "float"
** even if for internal computation a "real" is "double" (as per above).
*/
typedef enum {
  rndTypeUnknown,            // (0) no type known
  rndTypeShort,              // (1) short
  rndTypeFloat,              // (2) float
} rndType;

/*
** another enum, but for (output, 2D) image pixel types
*/
typedef enum {
  rndITypeUnknown,            // (0) no type known
  rndITypeUChar,              // (1) uchar
  rndITypeFloat,              // (2) float
  rndITypeDouble,             // (3) double
  rndITypeReal,               // (4) whatever "real" is
} rndIType;

/*
** enum rndSpace: the spaces mapped to by (transfer) functions defined
** by control points.  Either a 3-component color space (like mprColor
** in last project), or opacity (a scalar).
*/
typedef enum {
  rndSpaceUnknown,
  rndSpaceRGB,
  rndSpaceHSV,
  rndSpaceAlpha,
} rndSpace;

/*
** enum rndProbe: the different things that can be computed by rndRaySample at
** each sample along a ray (i.e. rndProbe uses "probe" as a noun rather than a
** verb).  The rndProbe values below are in dependency order: the larger
** probes (the probes with a larger value of the rndProbe enum) may depend on
** some of the smaller probes; this may inform how to structure the code to
** compute them.
**
** In particular, every probe *lower* than rndProbeInside is just about the
** ray position, independent of the volume being rendered.
**
** rndProbeInside depends on the convolution, but it is simply the value of
** rndConvo->inside after rndConvoEval(). The value or gradient from the the
** convolution is irrelevent for rndProbeInside.
**
** For rndProbeInside and lower probes, the compositing computed by
** rndRayBlend is over the sampling results at *all* positions along the ray
** where rndRaySample was called.  This makes sense because the value of these
** probes is defined everywhere.
**
** NOTE: For all the probes *higher* than rndProbeInside, the blending is only
** over those sampling results where rndConvo->inside is true (non-zero).
** This makes sense since values, gradients, transfer function evaluations,
** etc, are only meaningful when convolution has worked.
**
** Note that rndProbeLen(probe) returns the length of the probe value (e.g. 3
** for rndProbePosWorld, 1 for rndProbeValue, 4 for rndProbeRgba).
*/
typedef enum {
  rndProbeUnknown,   /* (0) */
  rndProbePosView,   /* (1) 3-vector: view-space position of probe
                        (for debugging ray-casting geometry */
  rndProbePosWorld,  /* (2) 3-vector: world-space position of probe
                        (for debugging ray-casting geometry */
  rndProbePosIndex,  /* (3) 3-vector: index-space position of probe
                        (for debugging ray-casting geometry */
  rndProbeInside,    /* (4) scalar: 1 if kernel support entirely inside data,
                        0 if outside (just like mprCtx->outside
                        being non-zero in Project 2) */
  rndProbeValue,     /* (5) scalar: value measured by convolution */
  rndProbeGradVec,   /* (6) 3-vector: gradient (in world-space) measured by
                        convolution */
  rndProbeGradMag,   /* (7) scalar: gradient magnitude */
  rndProbeRgba,      /* (8) 4-vector: RGBA from the transfer function: the
                        univariate LUT set by rndCtxSetTxf(). Also, if
                        rndCtxSetLevoy() is called, the opacity is also
                        multiplied by the Levoy opacity functions. The RGB
                        part of this is the "material color" of FSV 5.2 */
  rndProbeRgbaLit,   /* (9) 4-vector: rgba from transfer function (above)
                        shaded according to Blinn-Phong lighting, and then
                        multiplied by the lerp'd depth cueing color */
} rndProbe;

/*
** enum rndBlend: the different ways of compositing the information learned by
** probing. Not all blendings are useful for all probes, but probing and
** blending are orthogonal enough that specifying them separately helps
** simplify implementation.
**
** Your code may have switch or if statements for branching based on all
** possible probes, and other conditionals for branching based on all possible
** blends, but these should *not* be nested. Your code should reflect the
** (near) orthogonality of blending and probing or else you will lose style
** points.
*/
typedef enum {
  rndBlendUnknown,  /* (0) blending not known */
  rndBlendMax,      /* (1) max of all values (in vectors, per-component) */
  rndBlendSum,      /* (2) sum of all values (in vectors, per-component) */
  rndBlendMean,     /* (3) mean of all values (in vectors, per-component);
                       the sum divided by the number of samples; where "the
                       number of samples" depends on the probe: for
                       rndProbeInside or less: the number of sampling
                       locations on the ray. But for probes greater than
                       rndProbeInside: it is the number of samples inside (as
                       per rndConvo->inside) the volume */
  rndBlendOver,     /* (4) "over" operator of RGB color and alpha. This
                       breaks the orthogonality of blending and probing;
                       rndBlendOver only makes sense for rndProbeRgba
                       and rndProbeRgbaPhong) */
} rndBlend;

/*
** struct rndCamera: A container for information about how a camera is viewing
** some part of world-space.  The look-from, look-at, and up vectors are in
** world-space coordinates (using an orthonormal basis).
**
** All the fields are named here the same as in FSV Section 5.3, except
** for "WtoV" instead of "M_view".
*/
typedef struct {
  // "input" camera parameters, set either directly or via rndCameraSet
  real fr[3],      // look-from point
    at[3],         // look-at point
    up[3],         // up
    nc, fc,        /* near,far clip plane distances,
                      relative to the look-at point */
    FOV;           // vertical field-of-view, in degrees
  uint size[2];    /* # horz,vert samples of image plane; this
                      determines aspect ratio ("ar") */
  int ortho;       /* if non-zero, use orthographic instead of
                      perspective projection. This does not affect how
                      any variables below are calculated, but the
                      rndCamera is best place to store this info. */

  // "output" parameters set by rndCameraUpdate(), based on the above
  real ar,         // ((real)size[0])/size[1]
    d,             // distance between fr and at
    u[3],          // right-ward basis vector of view space
    v[3],          // up-ward basis vector of view space
    n[3],          // back-ward (into eye) basis vector of view space
    WtoV[16],      // world-to-view homogeneous coordinate transform
    VtoW[16],      // view-to-world transorm = inverse(WtoV)
    ncv, fcv,      /* (positive) distances, from eye to near and far
                      clipping planes */
    hght, wdth;    // height and width of visible image plane
} rndCamera;

/*
** struct rndVolume: A container for the oriented volume data that is
** processed by rendr.  Orientation is defined in terms of a particular
** orthonormal world-space basis: {left, posterior, superior}, which has
** significance for biomedical scans of human (or human-ish) anatomy, and can
** be understood as just some right-handed coordinate frame otherwise.
**
** Data is always scalar (we're not handling multiple values per voxel)
*/
typedef struct {
  uint size[3];             /* # samples along fastest [0] to slowest [2]
                               axes */
  real ItoW[16];            /* homogeneous coordinate mapping from
                               index-space (fast-to-slow ordering) to the
                               "left-posterior-superior" world-space */
  rndType dtype;            /* type of the data; determines which of the
                               union members below to use */
  union {                   /* union for the pointer to the image data;
                               the pointer values are all the same; this
                               is just to avoid casting.  The right union
                               member to use (data.ss vs data.fl) is known
                               at run-time from value of dtype */
    void *vd;
    signed short *ss;
    float *fl;              // really float, not "real"
  } data;
  Nrrd *nrrd;               /* nrrd struct around the same data;
                               not relevant for your work */
} rndVolume;

/*
** rndImage is a container for output images. No orientation information
** is saved. Implemented in image.c
*/
typedef struct {
  uint channel,        /* how many values are at each pixel; this always
                          fastest axis */
    size[2];           /* # of samples along faster (size[0]) and
                          slower (size[1]) spatial axes */
  rndIType dtype;      /* type of the data; determines which of the union
                          members below to use */
  union {              // union for the pointer to the image data
    void *vd;
    unsigned char *uc;
    float *fl;
    double *db;
    real *re;
  } data;
} rndImage;

/*
** The rndKernel stores everything about a reconstruction kernel. The kernel
** is non-zero only within [-support/2,support/2], for integer "support"
** which may be odd or even (but always positive). The kernels are set up at
** compile-time in such a way that each kernel knows its own derivative
** "deriv", so the derivative of kernel k is k->deriv. Implemented in
** kernel.c and kparse.c.
*/
typedef struct rndKernel_t {
  const char *name;                // short identifying string
  const char *desc;                // short descriptive string
  uint support;                    // # samples needed for convolution
  real (*eval)(real xx);           // how to evaluate the kernel
  const struct rndKernel_t *deriv; /* derivative of this kernel; will point
                                      back to itself when kernel is zero */
} rndKernel;

/*
** struct rndConvo: the buffers and state that are used to do convolution and
** to store its results. This is not in the rndCtx so that rndConvoEval can be
** thread-safe.  See go.c and ray.c for the context of how the rndConvo is
** passed to functions that may need to do convolution, like rndRaySample()
*/
typedef struct {
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real * kvalues1;
  real * kvalues2;
  real * kvalues3;
  real * the_kvalues1;
  real * the_kvalues2;
  real * the_kvalues3;
  int grad_need;
  // ( 4 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  real value,    // value reconstructed by convolution
    gradient[3]; // world-space gradient, by convolution
  int inside;    /* 1 if the last call to rndConvoEval was at a position where
                    all the required data samples (in the support of the
                    kernel) were available, --OR--, 0 if this is not the case.
                    Compared to previous projects, inside is the new
                    !outside */
} rndConvo;

struct rndCtx_t; /* forward declaration of struct, so that (pointers to)
                    functions within the struct can take a struct pointer
                    as an argument */
/*
** struct rndCtx: all the state (the "context") associated with doing volume
** rendering. See rendr_go.c to see how the context is created and set up to
** perform rendering.
**
** Note: rndCtxNeedGradient(const rndCtx *ctx) is a useful little
** function that looks at ctx->probe and ctx->blend (after they are set
** by rndCtxSetParms()) to see if the gradient (or gradient magnitude)
** needs to be computed by rndConvoEval().
*/
typedef struct rndCtx_t {
  const rndVolume *vol;   // volume to render
  const rndKernel *kern;  // kernel to reconstruct with

  int debuggingPixel;     // whether we're debugging a single pixel

  // basic rendering parameters given by rndCtxSetParms()
  real planeSep;          /* separation between sampling planes through the
                             view space volume; planeSep==0 says: "just take
                             one slice at the near clipping plane and be
                             done", OR, planeStep>0 says "take steps along
                             the ray so that all rays' N-th sample lie within
                             a plane N*planeSep behind the near clipping
                             plane". Thus, with perspective projection, the
                             step size on any single ray will be slightly
                             larger than this except at the very center of
                             the image plane */
  rndProbe probe;         /* what should we be probing or computing at each
                             sample along each ray */
  rndBlend blend;         /* how the probes along the ray wil be reduced down
                             to a single per-ray value or vector */
  real outsideValue;      /* the value that rndRayFinish should store in
                             ray->result, if there weren't any samples on the
                             ray that could be blended (e.g. when probing
                             rndProbeValue but the ray never went "inside"
                             the volume). If rndProbeLength(ctx->probe) > 1,
                             this outsideValue will be copied into each
                             channel. This can be NaN. */
  uint threadNum;         /* how many threads to render with; 0 for
                             non-threaded execution */
  int timing;             /* if non-zero: record per-ray timing (in
                             milliseconds) in last channel of output image */

  // camera and image specification, set by rndCtxSetCamera()
  rndCamera cam;

  /* For rndProbeRgba and rndProbeRgbaLit: univariate RGBA transfer function
     lookup table "rgba" and alpha threshold "alphaNear1"; set by
     rndCtxSetTxf() */
  struct {
    uint len;        // length of rgba lut
    real vmin, vmax; // min and max values represented by lut
    real *rgba;      // lookup table data, rgba on faster axis
    real unitStep;   /* the "unit" length to use when computing the opacity
                        correction as a function of ray step size; No
                        correction is needed when the ray step size is
                        exactly unitStep. By default this is 1.0 */
    real alphaNear1; /* with rndProbeRgba and rndProbeRgbaLit, finish ray if
                        its opacity exceeds this value < 1.0, or, use 1.0 to
                        say "no early ray termination" */
  } txf;

  /* Additional *optional* opacity function, which multiples the opacity
     generated by the univariate txf above: Levoy's "Isovalue contour
     surface" opacity functions (from 1988 "Display of Surfaces from Volume
     Data" paper; the triangular shapes in value/gradmag space). Parameters
     are set by rndCtxSetLevoy(). The evaluation of the opacity functions
     is done by rndTxfLevoy(). Whether to use this transfer function is
     testable by levoy.num being non-zero */
  struct {
    uint num;         /* number of surfaces; logical length of "vra" array.
                         0 means "no Levoy opacity functions in use" */
    real *vra;        /* 2-D array of parameters; v,r,a on faster axis;
                         logically a 1-D array of 3-vectors, where
                         each element of the 3-vector is:
                         v: isovalue, ("f_v" in Levoy paper),
                         r: fuzzy isocontour thickness,
                         a: max opacity of isocontour ("a_v" in paper) */
  } levoy;

  /* Specification of multiple directional lights. num, rgb, dir, and vsp are
     all set by rndCtxSetLight(); xyz is allocated by rndCtxSetLight() but the
     final values are set by rndCtxUpdateLight() (which is called in
     rndRender()).  Thus, by the time you have to compute lighting in ray.c,
     only the "num", "rgb", and "xyz" fields matter. */
  struct {
    uint num;         /* number of lights, logical length of each of the
                         arrays below */
    real *rgb;        /* light ii has color (rgb + 3*ii)[0,1,2]. Colors
                         are only in RGB space, not HSV */
    real *dir;        /* light direction (i.e. the direction *towards* the
                         light); this is not necessarily a normalized vector,
                         and it can be either in view-space or world-space */
    int *vsp;         /* light ii is in view-space if vsp[ii] if is non-zero,
                         else the light is in world-space */
    real *xyz;        /* normalized light direction in world space;
                         set by rndCtxUpdateLight() */
  } light;

  /* Blinn-Phong and depth-cueing lighting parameters that are
     specific to rndProbeRgbaLit; set by rndCtxSetLparm(). If light.num
     (above) is zero, the only lighting is ambient (controlled by ka) and
     depth cueing.  */
  struct {
    real ka, kd, ks, p;   /* Blinn-Phong parameters */
    real dcn[3], dcf[3];  /* depth-cueing: as the last step of computing
                             rgbaLit probe values, the RGB color values are
                             multiplied (per-channel) by a lerp (maybe with
                             V3_LERP) between dcn and dcf as the ray sampling
                             position varies between the near and far
                             clipping planes; Setting both dcn and dcf to
                             (1,1,1) effectively turns off depth-cueing */
  } lparm;

  /* The following section is for information determined solely by the volume
     and the kernel, independent of where convolution is being requested.
     Some of the state that you put in mprCtx in the last project should
     now be in rndConvo (above). Variables specific to multi-threaded
     operation should go here. */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real* WtoI;
  int upper;
  int lower;
  int odd_support;
  real* deriv_ItoW;
  const rndKernel * the_kern;

  // ( 9 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
} rndCtx;

/*
** struct rndTxf: stores a univariate transfer function ("txf"), consisting of
** "num" control points, all specified in output space "range", which may a
** 3-component color (rndSpaceRGB or rndSpaceHSV) or a scalar opacity
** (rndSpaceAlpha).  The control points are in the "data" array, which is
** logically:
**
** * if rndSpaceAlpha != range: a 4-by-num (faster-by-slower) array of reals,
**   or a length-"num" array of (X,A,B,C) 4-vectors, where X is location of
**   control point, and (A,B,C) are color coordinates.
**
** * if rndSpaceAlpha == range: a 2-by-num array of reals, or a length-"num"
**   array of (X,A) 2-vectors, where X is location of control point, and A is
**   opacity
**
** In either case, the X values have to strictly increase with control point
** index. In this way, the rndTxf is a lot like a mprCmap from Project 2.
** However, the rndTxf has no notion of a special color for being "outside"
** the volume, nor can a rndTxf contain a LUT (that's in rndCtx->txf.rgba)
*/
typedef struct {
  rndSpace space;      // what quantity is stored at each control point
  uint num;            // number of control points
  real *data;          // control point data
} rndTxf;

/*
** struct rndRay: all the state for a single ray being cast, including the
** results of sampling and intermediate blending while the ray is being cast,
** and the results of final blending when the ray is done. This also includes
** a record of whether *any* samples were "inside" the volume (as per
** rndConvo->inside)
**
** Note that there is no rndRayNew() and rndRayNix(). Users should be able to
** declare "rndRay ray" and then call rndRayStart(&ray, ...).  Thus,
** rndRayStart must completely initialize/reset the rndRay struct.  If
** rndRayStart does any dynamic allocation (though the reference
** implementation does none), then rndRayFinish must clean it up.
**
** Note that the largest possible value of rndProbeLen(ctx->probe) is 4.
** Knowing that, storage for probes (and their blending) need not be
** dynamically allocated according to the probe learned at run-time.
*/
typedef struct {
  uint hi, vi;     /* which pixel are we rendering: fast (horizontal)
                      and slow (vertical) index */
  /* suggestions for state to store here: what is the current position of the
     ray, what is the step size along the ray (for opacity correction), how
     many samples have been taken, and how many of those were inside the
     volume, what is the value and gradient that were learned by convolution,
     and what are the RGBA values from the transfer function, and from
     lighting. In other words: storage for every possible probe, as well as
     the extra information you'll need to correctly do all blendings. */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
    real startx,starty,startz;// start world position of the ray
  
  real uvn[3];//current view position of the ray
  int num_sample;
  int num_inside;
  int inside;

  real stepu,stepv,stepn;
  real step_size;

  real xyz_w[3];//current world position of the ray
  real xyz_i[3];//current index position of the ray
  real convo_val;
  real convo_grad[3];
  real convo_gradm;
  real txf_rgba[4];
  real lit_rgba[4];
  
  real current[4];
  real accum[4];
 
  // indicate whether the higher probes are initialized.
  int initialized;
  // ( 18 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  int finished;    // non-zero if rndRayFinished has been called
  real result[4];  /* final result of blending all probes along the ray, to
                      save into output image.  Number of values set here
                      should be rndProbeLen(ctx->probe), and the values
                      should be set by rndRayFinish */
  double time0,    /* when did this ray start (as double, not real, for
                      maximal precision */
    time1,         /* when did this ray finish */
    time;          /* time from ray start to ray finish (accounting if
                      necessary for multiple rays computed simultaneously by
                      one thread) */
} rndRay;

/* misc.c: miscellaneous things. Nothing for you to implement, but you
   can also add here functions that you want to use in other files.
   Such functions should be named starting with "rnd". */
#ifndef NDEBUG
extern int rndVerbose;
#else
#define rndVerbose 0
#endif
extern const int rndRealIsDouble;
extern const char *rndBiffKey;
#define RND rndBiffKey  // identifies library in error messages
extern const char *rndKeySpace;
extern int rndRendrMgo; // we're here via "rendr mgo"
extern void rndHSVfromRGB(real *hsv, const real *rgb);
extern void rndRGBfromHSV(real *rgb, const real *hsv);
extern real rndNan(unsigned short payload);
extern unsigned short rndNanPayload(real nval);
extern const airEnum *const rndType_ae;
extern const airEnum *const rndIType_ae;
extern const airEnum *const rndSpace_ae;
extern const airEnum *const rndBlend_ae;
extern const airEnum *const rndProbe_ae;
extern uint rndProbeLen(rndProbe probe);
extern int rndArray2dCheck(const Nrrd *nin, uint size0, int cent1);
extern int rndArrayToReal(Nrrd *nout, const Nrrd *nin);
extern uint rndSprintLenMax(uint oldmax, const real *vv, uint vvNum);
// new "public" functions here
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
// ( 0 lines in reference implementation )
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

/* kernel.c: gives compile-time definition of some reconstruction
   kernels and their derivatives (nothing for you to implement) */
extern const rndKernel *const rndKernelZero;
extern const rndKernel *const rndKernelBox;
extern const rndKernel *const rndKernelTent;
extern const rndKernel *const rndKernelBspln2;
extern const rndKernel *const rndKernelBspln3;
extern const rndKernel *const rndKernelCtmr;
extern const rndKernel *const rndKernelSpark;
extern const rndKernel *const rndKernelLuna;
extern const rndKernel *const rndKernelCelie;
extern const rndKernel *const rndKernelAll[];

/* kparse.c: for identifying kernels by strings (nothing for you to
   implement) */
extern const rndKernel *rndKernelParse(const char *kstr);
extern hestCB *rndKernelHest;

/* math.c: misc math things. NOTE: You finish rnd_m4_affine_inv() */
extern void rnd_m4_affine_inv(real inv[16], const real mat[16]);
extern void rnd_m4_print(FILE *file, const real mm[16]);

/* volume.c: for representing input volume data (nothing for you to do) */
extern int rndDataTypeRtoN(rndType dtype);
extern rndType rndDataTypeNtoR(int ntype);
extern rndVolume *rndVolumeNew(void);
extern rndVolume *rndVolumeNix(rndVolume *vol);
extern int rndVolumeLoad(rndVolume *vol, const char *fname);

/* image.c: for representing output image data (nothing for you to do) */
extern rndImage *rndImageNew();
extern rndImage *rndImageNix(rndImage *img);
extern int rndImageAlloc(rndImage *img, uint channel,
                         uint size0, uint size1,
                         rndIType dtype);
extern int rndImageNrrdWrap(Nrrd *nout, const rndImage *img);
extern int rndImageSave(const char *fname, const rndImage *img, real duration);
extern int rndImageLoad(rndImage *img, const char *fname);
extern int rndImageMinMax(real minmax[2], const rndImage *img);

/* camera.c: for working with the rndCamera to control the view of a volume,
   and the number of pixels in the output image.
   NOTE: You finish rndCameraUpdate() */
extern void rndCameraInit(rndCamera *cam);
extern rndCamera *rndCameraNew(void);
extern rndCamera *rndCameraNix(rndCamera *cam);
extern void rndCameraHestOptAdd(hestOpt **Ahopt, rndCamera *cam);
extern int rndCameraUpdate(rndCamera *cam);
extern int rndCameraSet(rndCamera *cam,
                        const real fr[3],
                        const real at[3],
                        const real up[3],
                        real nc, real fc,
                        real FOV,
                        uint size0, uint size1,
                        int ortho);
extern int rndCameraCheck(const rndCamera *cam);

/* txf.c: for working with transfer functions;
   NOTE: You finish rndTxfEval() and rndTxfLevoy() */
extern rndTxf *rndTxfNew(void);
extern rndTxf *rndTxfNix(rndTxf *txf);
extern int rndTxfAlloc(rndTxf *txf, rndSpace space, uint pointNum);
extern int rndTxfSave(const char *fname, const rndTxf *txf);
extern int rndTxfLoad(rndTxf *txf, const char *fname);
extern void rndTxfEval(real *out, const rndTxf *txf, real vv);
extern int rndTxfLutGenerate(Nrrd *nlut,
                             const rndTxf *ctxf, int rescaleC,
                             const rndTxf *atxf, int rescaleA,
                             uint num, real vmin, real vmax);
extern real rndTxfLevoy(const real *vra, uint num,
                        real val, real gradmag);

/* ctx.c: for setting up and using the rndCtx. You can augment rndCtxNew()
   and rndCtxNix(), as well as the rndCtx struct itself (above) */
extern rndCtx *rndCtxNew(const rndVolume *vol, const rndKernel *kern);
extern rndCtx *rndCtxNix(rndCtx *ctx);
extern int rndCtxCheckParms(real planeSep,
                            rndProbe probe, rndBlend blend,
                            uint tnum);
extern int rndCtxSetParms(rndCtx *ctx, real planeSep,
                          rndProbe probe, rndBlend blend,
                          real outsideVal,
                          uint tnum, int timing);
extern int rndCtxSetCamera(rndCtx *ctx, const rndCamera *cam);
extern int rndCtxSetTxf(rndCtx *ctx, const Nrrd *nlut,
                        real unitStep, real alphaNear1);
extern int rndCtxSetLevoy(rndCtx *ctx, const Nrrd *nlev);
extern int rndCtxSetLight(rndCtx *ctx, const Nrrd *nlight);
extern int rndCtxSetLparm(rndCtx *ctx,
                          real ka, real kd, real ks, real p,
                          const real dcn[3], const real dcf[3]);
extern int rndCtxUpdateLight(rndCtx *ctx);
extern int rndCtxNeedGradient(const rndCtx *ctx);

/* convo.c: convolution! NOTE: you implement rndConvoEval(), and augment
   rndConvoNew() and rndConvoNix(), as well as the definition of the rndConvo
   struct (above) */
extern rndConvo *rndConvoNew(const rndCtx *ctx);
extern rndConvo *rndConvoNix(rndConvo *cnv);
extern void rndConvoEval(real xw, real yw, real zw,
                         rndConvo *cnv, const rndCtx *ctx);

/* ray.c: breaks down the low-level work of volume rendering into various
   operations on the rndRay and the things learned by sampling along a ray.
   NOTE: you finish all functions here. */
extern void rndBlinnPhong(real rgbaOut[4], const real rgbaIn[4],
                          const real grad[3], const real Vdir[3],
                          const rndCtx *ctx);
extern void rndRayStart(rndRay *ray, uint hi, uint vi,
                        rndConvo *cnv, const rndCtx *ctx);
extern void rndRayBlendInit(rndRay *ray, const rndConvo *cnv, const rndCtx *ctx);
extern void rndRaySample(rndRay *ray, rndConvo *cnv, const rndCtx *ctx);
extern void rndRayBlend(rndRay *ray, const rndConvo *cnv, const rndCtx *ctx);
extern void rndRayStep(rndRay *ray, rndConvo *cnv, const rndCtx *ctx);
extern void rndRayFinish(rndRay *ray, const rndCtx *ctx);

/* go.c: the outer loop of the volume renderer; calls the rndRay functions
   above to cast one ray per pixel of output image.  */
extern int rndRender(rndImage *iout, rndCtx *ctx,
                     int debugPixH, int debugPixV);

/*
** The "rendr" executable is compiled from rendr.c, which relies on
** the various rnd_XCmd objects (each defined in rendr_X.c) for the
** rendr commands X={about, klist, ...}  None of those objects need
** declaring here, because rendr is the only thing using them, and
** rendr.c generates (via C pre-processor tricks) its own extern
** declarations of them.
**
** You should run "rrendr" to see all the commands, then run "rrendr
** about" to learn more about what work is required to finish this
** project. For rendr command X, you can look at rendr_X.c to see how
** it works, and to see the context in which the rnd library
** functions above are called.
*/

#ifdef __cplusplus
}
#endif
#endif /* RND_HAS_BEEN_INCLUDED */
