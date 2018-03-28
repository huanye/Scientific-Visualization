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
static real lerp(real omin, real omax, real imin, real xx, real imax) {
  real ret=0;
  real alpha = (xx-imin)/(imax-imin);
  ret = (1-alpha)*omin+alpha*omax;
  return ret;
}

static void V3Quantize(unsigned char *rgbOut,real* min, real const* source, real* max, uint num){
  for(int i=0;i<3;i++)
    rgbOut[i] = mprQuantize(min[i],source[i],max[i],num);
}
// ( 11 lines in reference implementation )
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

/*
** These two functions convert between RGB and the HSV color single hexcone.
** http://en.wikipedia.org/wiki/HSL_and_HSV
**
** All inputs and outputs should be in the range [0,1]; no range checking or
** clamping is no done.  The Hue values from 0 to 1 are interpreted as angles
** from 0 to 360 degrees. The RGB value for both HSV=(0,1,1) and HSV=(1,1,1)
** is (1,0,0).
*/
void
mprHSVfromRGB(real *hsv, const real *rgb) {
  float _hsv[3];
  dyeRGBtoHSV(_hsv+0, _hsv+1, _hsv+2,
              rgb[0], rgb[1], rgb[2]);
  V3_COPY(hsv, _hsv);
}

void
mprRGBfromHSV(real *rgb, const real *hsv) {
  float _rgb[3];
  dyeHSVtoRGB(_rgb+0, _rgb+1, _rgb+2,
              hsv[0], hsv[1], hsv[2]);
  V3_COPY(rgb, _rgb);
}

void
mprToRGB(real *rgb, const real *col, mprColor space) {
  if (mprColorRGB == space) {
    V3_COPY(rgb, col);
  } else if (mprColorHSV == space){
    mprRGBfromHSV(rgb, col);
  } else {
    V3_SET_NAN(rgb);
  }
  return;
}
/*
** mprPictureItoW sets the ItoW orientation matrix of a cell-centered sampling
** grid, based on the following specification of a viewing rectangle:
**
** size0, size1: number samples along faster, slower axes of grid
** (centX, centY): center of grid
** fov: extent along the slower axis, assuming |scl1|==1 and shear0==0
** scl0, scl1: non-isotropic scaling (can be negative to flip
**   direction of axis traversal) on faster, slower axes
** angle: angle (in degrees) between faster axis and world X axis,
**   assuming shear0==0
** shear0, shear1: amount by which faster, slower axes should be sheared
**   along the other axis
**
** This also does some useful error checking on all these parameters.
*/
int
mprPictureItoW(real ItoW[9],
               uint size0, uint size1,
               real centX, real centY, real fov,
               real scl0, real scl1,
               real angle, real shear0, real shear1) {
  if (!ItoW) {
    biffAddf(MPR, "%s: got NULL pointer", __func__);
    return 1;
  }
  if (!( size0 > 0 && size1 > 0 )) {
    biffAddf(MPR, "%s: need two positive sizes (not %u,%u)",
             __func__, size0, size1);
    return 1;
  }
  if (!( isfinite(centX) && isfinite(centY) && isfinite(fov) )) {
    biffAddf(MPR, "%s: center X,Y (%g,%g) or FOV (%g) bad", __func__,
             centX, centY, fov);
    return 1;
  }
  if (!( isfinite(scl0) && scl0 != 0 &&
         isfinite(scl1) && scl1 != 0 )) {
    biffAddf(MPR, "%s: sampling scalings (%g,%g) bad", __func__, scl0, scl1);
    return 1;
  }
  if (!( isfinite(angle) && isfinite(shear0) && isfinite(shear1) )) {
    biffAddf(MPR, "%s: angle (%g) or shears (%g,%g) bad",
             __func__, angle, shear0, shear1);
    return 1;
  }

  real cc = cos(M_PI*angle/180);
  real ss = sin(M_PI*angle/180);
  real spc = fov/size1;    // nominal spacing between samples
  real _edge0[2] = {spc*cc, spc*ss};
  real _edge1[2] = {spc*ss, spc*(-cc)};
  // edge0 = _edge0 + shear0*_edge1
  // edge1 = _edge1 + shear1*_edge0
  real edge0[2] = {_edge0[0] + shear0*_edge1[0], _edge0[1] + shear0*_edge1[1]};
  real edge1[2] = {_edge1[0] + shear1*_edge0[0], _edge1[1] + shear1*_edge0[1]};
  // at first, don't worry about center; put first sample at origin
  M3_SET(ItoW,
         scl0*edge0[0], scl1*edge1[0], 0,
         scl0*edge0[1], scl1*edge1[1], 0,
         0,             0,             1);
  // half-way through index space
  real ipos[3] = {(size0-1)/2.0, (size1-1)/2.0, 1};
  real wpos[3];
  MV3_MUL(wpos, ItoW, ipos);
  // then fix origin so half-way through index-space lands at center
  ItoW[2] = centX-wpos[0];
  ItoW[5] = centY-wpos[1];
  return 0;
}

/*
******** mprPictureRGB: Color convolution results when ctx->mode is either
** mprModeValueCmap or mprModeValueCmapShaded).
**
******** INPUT: mprCtx *ctx, just after a call to mprConvoEval
**
******** OUTPUT: unsigned char *rgbOut, the final 8-bit pixel color assigned
** to the convolution result, determined by colormapping, possibly shading
** (with mprModeValueCmapShaded), and quantizing to 8 bits.
**
** Let V0 = ctx->value.  If V0 is not finite (because mprConvoEval could not
** compute the convolution sum), rgbOut is ctx->cmap->outside, converted to
** RGB, and quantized to 8-bits, and we're done (regardles of whether
** ctx->mode is mprModeValueCmap or mprModeValueCmapShaded).
**
** Otherwise, we lerp V0 from [ctx->imgMinMax[0], ctx->imgMinMax[1]]
** (range of values to colormap) to the domain of cmap, to get a new
** scalar V1, and then we get initial color rgb1 by colormapping V1.
** -- If cmap->lut is NULL: the colormap domain is between the first and
**    last control point locations in cmap, and rgb1 is set by
**    mprCmapEval(rgb1, ctx->cmap, V).
** -- If cmap->lut is not NULL (someone already called mprCmapLutGen), the
**    colormap domain is [cmap->lutMin,cmap->lutMax] (which should actually
**    match the domain from the control points), and rgb1 is set by
**    quantizing V1 to the colormap domain divided into cmap->lutLen
**    segments, and then setting rgb1 copying the right RGB color from
**    cmap->lut.
** If ctx->mode is mprModeValueCmap, we save into rgbOut the result of
** quantizing rgb1 to 8-bits, and we're done.
**
** Otherwise (with mprModeValueCmapShaded), we compute from the 2-D
** world-space image gradient (call it G), a world-space vector N that is
** always perpendicular to (i.e. normal to) a graph of the convolution
** results.  Specifically, define a 3D surface S as the plot of function
** Z(X,Y)=(ctx->zscl)*f(X,Y) where f(X,Y) is ctx->value after
** mprConvoEval(ctx, X, Y). N is a unit-length vector logically located at
** (X,Y,Z) (i.e. a point on S), which is perpendicular to S.  If G=(0,0) then
** N=(0,0,1). Note that the scaling of S along Z due to ctx->zscl affects the
** shading, but not the colormapping by ctx->cmap to set rgb1. You have to
** determine a formula for N as a function of G and ctx->zscl; use the
** Pythagorean theorem.
**
** Then we compute diffuse shading D = max(0,dot(N, L)), where
** L=ctx->ldir. With full shading (1==ctx->shading), D scales the
** floating-point color rgb1 to make a new color rgb2 = D*rgb1.  With no
** shading (0==ctx->shading), the rgb1 computed for mprModeValueCmap is
** unchanged, rgb2 = rgb1.  Thus ctx->shading controls lerping between the
** not-shaded and the fully-shaded result. Save into rgbOut the result of
** quantizing rgb2 to 8-bits.
*/
void
mprPictureRGB(unsigned char *rgbOut, const mprCtx *ctx) {
  const mprCmap *cmap = ctx->cmap;
  V3_SET(rgbOut, 100, 80, 0); // initialize
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real v = ctx->value;
  real lmin = cmap->data[0];
  real lmax = cmap->data[0+4*(cmap->num-1)];

  real min[3]; 
  V3_SET(min,0,0,0);
  real max[3]; 
  V3_SET(max,1,1,1);
  real rgb1[3];
  if(isnan(v)){
    // quantizing cmap->outside to 8-bit pixel as rgbOut
    mprToRGB(rgb1, cmap->outside, cmap->space);
    V3Quantize(rgbOut,min, rgb1, max, 256);
  }
  
  else{
     real v1 = lerp(lmin, lmax, ctx->imgMinMax[0], v, ctx->imgMinMax[1]);
     if(cmap->lut){
       int index = mprQuantize(cmap->lutMin,v1,cmap->lutMax,cmap->lutLen);
       for(int i=0;i<3;i++)
        rgb1[i] = (cmap->lut+3*index)[i];
     }
     else{
       mprCmapEval(rgb1, cmap, v1);
     }

     if(ctx->mode==mprModeValueCmap){  
     // quantizing rgb1 to 8-bit pixel as rgbOut
      V3Quantize(rgbOut,min, rgb1, max, 256);
      }
     else if(mprModeValueCmapShaded == ctx->mode){
       real N[3];
       real tmp;
       // normal vector N before normalizing
       N[0] = -ctx->zscl*ctx->gradient[0];
       N[1] = -ctx->zscl*ctx->gradient[1];
       N[2] = 1;
       //normalize N vector N after normalizing
       V3_NORM(N, N, tmp);
       real dot_tmp = V3_DOT(N, ctx->ldir);
       real D=((dot_tmp>=0)?dot_tmp:0);
       //lerp from 1 to D as the coefficient of rgb1
       real coeff = lerp(1, D, 0, ctx->shading, 1);
       real rgb2[3];

       V3_SCALE(rgb2, coeff, rgb1);

       // quantizing rgb1 to 8-bit pixel as rgbOut
       V3Quantize(rgbOut,min, rgb2, max, 256);
       
       }
     }
  // ( 36 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return;
}

/*
** mprPictureEval: sets the pixidx-th pixel of output image omg, according to
** convolution results and colormapping/shading parameters in ctx.
*/
void
mprPictureEval(mprImage *omg, uint pixidx, const mprCtx *ctx) {
  unsigned short outside = (unsigned short)(ctx->outside);
  switch (ctx->mode) {
  case mprModeWorldPos:
    V2_COPY(omg->data.fl + 2*pixidx, ctx->wpos);
    break;
  case mprModeIndexPos:
    V2_COPY(omg->data.fl + 2*pixidx, ctx->ipos);
    break;
  case mprModeOutside:
    omg->data.fl[pixidx] = outside;
    break;
  case mprModeValue:
    if (outside) {
      omg->data.fl[pixidx] = mprNan(outside);
    } else {
      omg->data.fl[pixidx] = ctx->value;
    }
    break;
  case mprModeGradient:
    if (outside) {
      V2_SET(omg->data.fl + 2*pixidx, mprNan(outside), mprNan(outside));
    } else {
      V2_COPY(omg->data.fl + 2*pixidx, ctx->gradient);
    }
    break;
  case mprModeFuzzyIso:
    if (outside) {
      omg->data.fl[pixidx] = 0;
    } else {
      real eps = (fabs(ctx->value - ctx->fuzzyIsoVal)
                  / (V2_LEN(ctx->gradient) * ctx->fuzzyIsoThick));
      omg->data.fl[pixidx] = MAX(0, 1 - eps);
    }
    break;
  case mprModeValueCmap:
  case mprModeValueCmapShaded: {
    unsigned char rgb[3];
    mprPictureRGB(rgb, ctx);
    V3_COPY(omg->data.uc + 3*pixidx, rgb);
    break;
  }
  case mprModeUnknown:
  default:
    // error checking for this done by mprPictureSample
    break;
  }
  return;
}

/*
** mprPictureSample: Computes a final "picture": an image of convolution
** results, possibly colormapped and shaded (according to ctx->mode). The
** centX through angle parameters determine the sampling grid of the output
** image, as determined by mprPictureItoW
**
** This calls, in the inner loop, first mprConvoEval to do the convolution,
** and then mprPictureEval to figure out the final output pixel value.
*/
int
mprPictureSample(mprImage *omg, mprCtx *ctx,
                 uint size0, uint size1,
                 real centX, real centY, real fov,
                 real angle) {
  if (!(omg && ctx)) {
    biffAddf(MPR, "%s: got NULL pointer (%p,%p)", __func__,
             (void*)omg, (void*)ctx);
    return 1;
  }
  if (mprPictureItoW(omg->ItoW,
                     size0, size1,
                     centX, centY, fov,
                     1 /* scl0 */, 1 /* scl1 */,
                     angle, 0 /* shear0 */, 0 /* shear1 */)) {
    biffAddf(MPR, "%s: problems with viewing rectangle parms", __func__);
    return 1;
  }

  /* depending on ctx->mode, allocate output image (and do error checking on
   having a colormap, if its needed) */
  uint channel;
  mprType mtype;
  switch (ctx->mode) {
  case mprModeWorldPos:
  case mprModeIndexPos:
  case mprModeGradient:
    channel = 2;
    mtype = mprTypeFloat;
    break;
  case mprModeValue:
  case mprModeOutside:
  case mprModeFuzzyIso:
    channel = 1;
    mtype = mprTypeFloat;
    break;
  case mprModeValueCmap:
  case mprModeValueCmapShaded:
    if (!ctx->cmap) {
      biffAddf(MPR, "%s: can't generate %s-mode picture without a colormap",
               __func__, airEnumStr(mprMode_ae, ctx->mode));
      return 1;
    }
    if (mprModeValueCmapShaded == ctx->mode) {
      if (!( V3_ISFINITE(ctx->ldir) && isfinite(ctx->shading) )) {
        biffAddf(MPR, "%s: light direction and shading not set", __func__);
        return 1;
      }
    }
    channel = 3;
    mtype = mprTypeUChar;
    break;
  case mprModeUnknown:
  default:
    biffAddf(MPR, "%s: mode %d unknown!", __func__, ctx->mode);
    return 1;
  }

  // allocate output image, and set orientation
  if (mprImageAlloc(omg, channel, size0, size1, mtype)) {
    biffAddf(MPR, "%s: trouble allocating output image", __func__);
    return 1;
  }
  // do sampling of final output raster picture
  double time0 = airTime();
  for (uint i1=0; i1<size1; i1++) {
    for (uint i0=0; i0<size0; i0++) {
      real ipos[3], wpos[3];
      V3_SET(ipos, i0, i1, 1);
      MV3_MUL(wpos, omg->ItoW, ipos);
      mprConvoEval(ctx, wpos[0], wpos[1]);
      mprPictureEval(omg, i0 + size0*i1, ctx);
    }
  }
  double time1 = airTime();
  // save sampling rate in kHz
  ctx->srate = size0*size1/(1000*(time1 - time0));
  return 0;
}
