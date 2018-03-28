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

#include "mpr.h"
#include "mprPrivate.h"

/* You may (if you want) define new functions or #define macros here;
   but our intent is only that they would be used in this file.  Thus,
   prefix any function declarations here with "static" */
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
// ( 0 lines in reference implementation )
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

/*
** mprCtxNew creates the context to contain all resources and state
** associated with computing on a given image "img", reconstruction kernel
** "kern", to be used for the purpose described by "mode".  Depending on
** "mode", this may also need a colormap "cmap", which works over values in
** the range [imm[0],imm[1]].  None of these things change over the lifetime
** of this mprCtx. Note from the error handling that the image must be a
** scalar (1-channel) image of floats, so that is the only case to handle
** inside your mprConvoEval code.
*/
mprCtx *
mprCtxNew(const mprImage *img, const mprKernel *kern,
          mprMode mode, const mprCmap *cmap, const real *imm) {
  if (!(img && kern)) { // cmap, imm can be NULL
    biffAddf(MPR, "%s: got NULL pointer (%p,%p)", __func__,
             (void*)img, (void*)kern);
    return NULL;
  }
  if (1 != img->channel) {
    biffAddf(MPR, "%s: only works on scalar (not %u-channel) images",
             __func__, img->channel);
    return NULL;
  }
  if (mprTypeFloat != img->dtype) {
    biffAddf(MPR, "%s: only works on %s (not %s) pixels", __func__,
             airEnumStr(mprType_ae, mprTypeFloat),
             airEnumStr(mprType_ae, img->dtype));
    return NULL;
  }
  if (airEnumValCheck(mprMode_ae, mode)) {
    biffAddf(MPR, "%s: %s %d not valid", __func__, mprMode_ae->name, mode);
    return NULL;
  }
  if (mprModeValueCmap == mode
      || mprModeValueCmapShaded == mode) {
    if (!(cmap && imm)) {
      biffAddf(MPR, "%s: need non-NULL cmap and imm (not %p,%p) "
               "for mode %s\n", __func__,
               (void*)cmap, (void*)imm, airEnumStr(mprMode_ae, mode));
      return NULL;
    }
    if (!( isfinite(imm[0]) && isfinite(imm[1]) )) {
      biffAddf(MPR, "%s: got non-finite image range [%g,%g]",
               __func__, imm[0], imm[1]);
      return NULL;
    }
    if (imm[0] == imm[1]) {
      biffAddf(MPR, "%s: need non-zero image range (not [" RCS "," RCS "])",
               __func__, imm[0], imm[1]);
      return NULL;
    }
  } else {
    if (cmap || imm) {
      biffAddf(MPR, "%s: got non-NULL cmap or imm (%p,%p) but "
               "not needed for mode %s\n", __func__,
               (void*)cmap, (void*)imm,
               airEnumStr(mprMode_ae, mode));
      return NULL;
    }
  }
  mprCtx *ctx = MALLOC(1, mprCtx);
  assert(ctx);
  ctx->verbose = 0;
  ctx->image = img;
  ctx->kern = kern;
  ctx->mode = mode;
  ctx->cmap = cmap; // could be NULL
  ctx->imgMinMax[0] = imm ? imm[0] : mprNan(0);
  ctx->imgMinMax[1] = imm ? imm[1] : mprNan(0);
  V3_SET(ctx->ldir, mprNan(0), mprNan(0), mprNan(0));
  ctx->shading = mprNan(0);
  ctx->zscl = 1;
  ctx->fuzzyIsoVal = mprNan(0);
  ctx->fuzzyIsoThick = mprNan(0);

  /* Here, you can set up any additional resources (setting variables or
     allocating arrays) that will be used in mprConvoEval(), mprPictureRGB(),
     or mprIsocontour.  In mprCtxNix(), free() anything malloc()ed here */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  ctx->WtoI = MALLOC(9,real);
  real tmp;
  M3_INVERSE(ctx->WtoI, ctx->image->ItoW, tmp);
  int support = (signed)(ctx->kern->support);
  ctx->odd_support = (support%2);
  ctx->upper = ctx->odd_support?(support-1)/2:support/2;
  ctx->lower = ctx->odd_support?(1-support)/2:1-support/2;
  ctx->deriv_ItoW = MALLOC(4,real);
  M23_INVERSE_TRANSPOSE(ctx->deriv_ItoW, ctx->image->ItoW, tmp);  
  // ( 26 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

  return ctx;
}

mprCtx *
mprCtxNix(mprCtx *ctx) {

  if (ctx) {
    // free anything that was allocated (on the heap) by mprCtxNew
    // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
    free(ctx->WtoI);
    free(ctx->deriv_ItoW);
    // ( 10 lines in reference implementation )
    // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
    free(ctx);
  }
  return NULL;
}

/*
** mprCtxModeValueCmapShadedParmSet: sets the additional parameters needed
** for mode mprModeValueCmapShaded.
**
** Nothing for you to implement.
*/
int
mprCtxModeFuzzyIsoParmSet(mprCtx *ctx, real isoval, real thickness) {
  if (!ctx) {
    biffAddf(MPR, "%s: got NULL pointer", __func__);
    return 1;
  }
  if (mprModeFuzzyIso != ctx->mode) {
    biffAddf(MPR, "%s: ctx->mode already set to %s; we need %s", __func__,
             airEnumStr(mprMode_ae, ctx->mode),
             airEnumStr(mprMode_ae, mprModeFuzzyIso));
    return 1;
  }
  if (!(isfinite(isoval) && isfinite(thickness))) {
    biffAddf(MPR, "%s: got non-finite isoval %g or thickness %g", __func__,
            isoval, thickness);
    return 1;
  }
  if (!( thickness > 0 )) {
    biffAddf(MPR, "%s: need thickness > 0 (not %g)", __func__, thickness);
    return 1;
  }
  ctx->fuzzyIsoVal = isoval;
  ctx->fuzzyIsoThick = thickness;
  return 0;
}

/*
** mprCtxModeValueCmapShadedParmSet: sets the additional parameters needed
** for mode mprModeValueCmapShaded.  Normalizes the light direction vector.
**
** Nothing for you to implement.
*/
int
mprCtxModeValueCmapShadedParmSet(mprCtx *ctx, const real _ldir[3],
                                 real shading, real zscl) {
  if (!(ctx && _ldir)) {
    biffAddf(MPR, "%s: got NULL pointer (%p,%p)", __func__,
             (void*)ctx, (void*)_ldir);
    return 1;
  }
  if (mprModeValueCmapShaded != ctx->mode) {
    biffAddf(MPR, "%s: ctx->mode already set to %s; we need %s", __func__,
             airEnumStr(mprMode_ae, ctx->mode),
             airEnumStr(mprMode_ae, mprModeValueCmapShaded));
    return 1;
  }
  real ldir[3];
  V3_COPY(ldir, _ldir);
  real llen = V3_LEN(ldir);
  if (!llen) {
    biffAddf(MPR, "%s: got zero-length light direction", __func__);
    return 1;
  }
  V3_SCALE(ctx->ldir, 1/llen, ldir);
  ctx->shading = shading;
  ctx->zscl = zscl;
  return 0;
}
