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
// ( 1 lines in reference implementation )
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

/*
******** mprCmapEval: use the control points in a colormap to
** color a single scalar value
**
******** INPUT:
** cmap: the colormap.  Even if cmap has a lut set inside
** (via mprCmapLutGen), the colormapping is done via sequence of
** control points that defines the true colormap.
**
** val: the value to colormap. May or may not be finite, but you
** only handle the finite case.
**
******** OUTPUT:
** rgb: an RGB color determined by mapping val through the colormap, as
** follows. First, clamp V to the interval determined by the first and last
** control points.  Then, the pair of colormap control points bracketing
** (possibly clamped) V must be found, and then lerp() should be used to
** interpolate colors *component-wise* between the colors at the two relevant
** control points. If mprColorHSV == ctx->cmap->space, then interpolated
** colors must be converted (with mprToRGB) to RGB color space. So the order
** of operations is: interpolate the colors (per component), and then convert
** to RGB.
**
** <For 33710 students> interpolating HSV colors: Care must be taken
** when interpolating hue angles (interpreting hue 0 to 1 as angle 0 to 2*pi
** radians within a circle), so that between any two hues, hue0 and hue1, the
** interpolation happens along the *shorter* of the two circular segments
** connecting hue0 and hue1 (i.e take the short way around the circle). For
** example, interpolating hues 0.1 and 0.9 should pass through 0.0==1.0,
** rather than through 0.5. </For 33710 students>
*/
void
mprCmapEval(real *rgb, const mprCmap *cmap, real val) {
  if (!rgb || !cmap) {
    return;
  }
  if (!isfinite(val)) {
    mprToRGB(rgb, cmap->outside, cmap->space);
    return;
  }
  // else have a finite value to colormap
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  uint index = 0;
  real alpha = 0.0;

  real lmin = cmap->data[0];
  real lmax = cmap->data[0 + 4*(cmap->num-1)];
  
  if(val>=lmax){
    index = cmap->num-1;
    alpha = 0.0;
  }
  else if(val<=lmin){
     index = 0;
     alpha = 0.0;
  }

  else{
   for(int i=0;i<(signed)cmap->num-1;i++)
    if(cmap->data[4*(i+1)]>val){
          index = i; 
          alpha = (val-cmap->data[4*i])/(cmap->data[4*(i+1)]-cmap->data[4*i]);
          break;
   }
 }
  real color[3];
  if(cmap->space == mprColorHSV){
    real hue0 = cmap->data[4*index+1];
    real hue1 = cmap->data[1+4*(index+1)];

    if(hue0>=hue1 && hue0 - hue1>0.5)
         hue1+=1.0;
    else if(hue0<hue1 && hue1 - hue0>0.5)
         hue0+=1.0;
    color[0] = (1-alpha)*hue0+alpha*hue1;
    if(color[0]>=1.0) 
      color[0]-=1.0;
  }
  else if(cmap->space == mprColorRGB){
    color[0] = (1-alpha)*cmap->data[4*index+1]+alpha*cmap->data[1+4*(index+1)];
  }
  color[1] = (1-alpha)*cmap->data[4*index+2]+alpha*cmap->data[2+4*(index+1)];
  color[2] = (1-alpha)*cmap->data[4*index+3]+alpha*cmap->data[3+4*(index+1)];

  mprToRGB(rgb, color, cmap->space);
  
  // ( 42 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return;
}

/*
** Equip a cmap with a look-up table.  To be consistent with how we name
** things, this also could have been called mprCmapSample, given how it
** repeatedly called mprCmapEval.
*/
int
mprCmapLutGen(mprCmap *cmap, uint len) {
  if (!cmap) {
    biffAddf(MPR, "%s: got NULL pointer", __func__);
    return 1;
  }
  if (cmap->lut) { // free any existing LUT
    free(cmap->lut);
    cmap->lut = NULL;
    cmap->lutMin = cmap->lutMax = mprNan(0);
    cmap->lutLen = 0;
  }
  if (len) {
    cmap->lut = MALLOC(3*len, real);
    if (!cmap->lut) {
      biffAddf(MPR, "%s: couldn't allocate length-%u lut", __func__, len);
      return 1;
    }
    cmap->lutLen = len;
    real lmin = cmap->lutMin = cmap->data[0];
    real lmax = cmap->lutMax = cmap->data[0 + 4*(cmap->num-1)];
    for (uint li=0; li<len; li++) {
      real aa = ((real)li + 0.5)/((real)len);
      mprCmapEval(cmap->lut + 3*li, cmap, (1-aa)*lmin + aa*lmax);
    }
  }
  return 0;
}

/*
** render a cmap->lut as an 8-bit RGB image
*/
int
mprCmapLutDraw(mprImage *img, const mprCmap *cmap, uint height) {
  if (!(img && cmap && height)) {
    biffAddf(MPR, "%s: got NULL pointers (%p,%p) or zero height (%u)",
             __func__, (void*)img, (void*)cmap, height);
    return 1;
  }
  if (!cmap->lut) {
    biffAddf(MPR, "%s: cmap doesn't have lut (call mprCmapLutGen)", __func__);
    return 1;
  }
  uint width = cmap->lutLen;
  if (mprImageAlloc(img, 3, width, height, mprTypeUChar)) {
    biffAddf(MPR, "%s: couldn't allocate output image", __func__);
    return 1;
  }
  for (uint hi=0; hi<width; hi++) {
    unsigned char rgbUC[3];
    const real *rgb = cmap->lut + 3*hi;
    V3_SET(rgbUC,
           (unsigned char)mprQuantize(0, rgb[0], 1, 256),
           (unsigned char)mprQuantize(0, rgb[1], 1, 256),
           (unsigned char)mprQuantize(0, rgb[2], 1, 256));
    for (uint vi=0; vi<height; vi++) {
      V3_COPY(img->data.uc + 3*(hi + width*vi), rgbUC);
    }
  }
  return 0;
}


// for storing colormap extras in NRRD key/value pairs
#define KEY_SPACE "space"
#define KEY_OUTSIDE "outside"

mprCmap *
mprCmapNew(void) {
  mprCmap *cmap = MALLOC(1, mprCmap);
  assert(cmap);
  cmap->space = mprColorUnknown;
  V3_SET_NAN(cmap->outside);
  cmap->num = 0;
  cmap->data = NULL;
  cmap->lut = NULL;
  cmap->lutMin = cmap->lutMax = mprNan(0);
  cmap->lutLen = 0;
  return cmap;
}

mprCmap *
mprCmapNix(mprCmap *cmap) {
  if (cmap) {
    if (cmap->data) free(cmap->data);
    if (cmap->lut) free(cmap->lut);
    free(cmap);
  }
  return NULL;
}

/* inside cmap, allocate pointNum control points in color-space
   space and outside color ocol */
int
mprCmapAlloc(mprCmap *cmap, mprColor space,
             const real ocol[3],
             uint pointNum) {
  if (!( cmap && ocol )) {
    biffAddf(MPR, "%s: got NULL pointer (%p,%p)", __func__,
             (void*)cmap, (void*)ocol);
    return 1;
  }
  if (!pointNum) {
    biffAddf(MPR, "%s: need non-zero number of points", __func__);
    return 1;
  }
  if (!( mprColorUnknown < space && space < mprColorLast )) {
    biffAddf(MPR, "%s: space %d not valid", __func__, space);
    return 1;
  }
  if (!V3_ISFINITE(ocol)) {
    biffAddf(MPR, "%s: outside color %g,%g,%g not set", __func__,
             ocol[0], ocol[1], ocol[2]);
    return 1;
  }
  if (cmap->data) {
    free(cmap->data);
  }
  cmap->data = MALLOC(4*pointNum, real);
  if (!(cmap->data)) {
    biffAddf(MPR, "%s: couldn't alloc data for %u real 4-vecs",
             __func__, pointNum);
    return 1;
  }
  cmap->space = space;
  V3_COPY(cmap->outside, ocol);
  cmap->num = pointNum;
  return 0;
}

/* static because even if this is generally useful;
   it is only needed inside this file */
static int
mprNrrdCmapCheck(const Nrrd *nin) {
  if (nrrdCheck(nin)
      || nrrd1DIrregMapCheck(nin)) {
    biffMovef(MPR, NRRD, "%s: control point array not ok", __func__);
    return 1;
  }
  if (!(nrrdTypeFloat == nin->type
        || nrrdTypeDouble == nin->type)) {
    biffAddf(MPR, "%s: need %s or %s type data, not %s", __func__,
             airEnumStr(nrrdType, nrrdTypeFloat),
             airEnumStr(nrrdType, nrrdTypeDouble),
             airEnumStr(nrrdType, nin->type));
    return 1;
  }
  if (4 != nin->axis[0].size) {
    biffAddf(MPR, "%s: axis[0].size should be 4, not %u", __func__,
             (uint)nin->axis[0].size);
    return 1;
  }
  if (!isfinite(nrrdFLookup[nin->type](nin->data, 0))) {
    biffAddf(MPR, "%s: first control point location not set", __func__);
    return 1;
  }
  /* checking the space:= and outside:= keys can be
     done by whomever needs the info */
  return 0;
}

// load a colormap from a file
int
mprCmapLoad(mprCmap *cmap, const char *fname) {
  if (!(cmap && fname)) {
    biffAddf(MPR, "%s: got NULL pointer (%p,%p)", __func__,
             (void*)cmap, (void*)fname);
    return 1;
  }
  airArray *mop = airMopNew();
  Nrrd *nin = nrrdNew();
  airMopAdd(mop, nin, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdLoad(nin, fname, NULL)) {
    biffMovef(MPR, NRRD, "%s: trouble reading file", __func__);
    airMopError(mop);
    return 1;
  }
  if (mprNrrdCmapCheck(nin)) {
    biffAddf(MPR, "%s: given array doesn't conform to a mapr colormap", __func__);
    airMopError(mop);
    return 1;
  }
  char *spcS = airToUpper(nrrdKeyValueGet(nin, KEY_SPACE));
  if (!spcS) {
    biffAddf(MPR, "%s: color space not set via \"%s:=<space>\" "
             "key/value pair", __func__, KEY_SPACE);
    airMopError(mop);
    return 1;
  }
  airMopAdd(mop, spcS, airFree, airMopAlways);
  int cspace;
  if (!strcmp("RGB", spcS)) {
    cspace = mprColorRGB;
  } else if (!strcmp("HSV", spcS)) {
    cspace = mprColorHSV;
  } else {
    biffAddf(MPR, "%s: color space %s not recognized", __func__, spcS);
    airMopError(mop);
    return 1;
  }
  char *outS = nrrdKeyValueGet(nin, KEY_OUTSIDE);
  airMopAdd(mop, outS, airFree, airMopAlways);
  real outside[3];
  double outsideD[3];
  if (!( outS && 3 == sscanf(outS, "%lf %lf %lf",
                             outsideD+0, outsideD+1, outsideD+2) )) {
    biffAddf(MPR, "%s: outside color not set via \"%s:=%s\" key/value pair",
             __func__, KEY_OUTSIDE, outS ? outS : "(null)");
    airMopError(mop);
    return 1;
  }
  outside[0] = (real)outsideD[0];
  outside[1] = (real)outsideD[1];
  outside[2] = (real)outsideD[2];
  if (mprCmapAlloc(cmap, cspace, outside,
                   (uint)nin->axis[1].size)) {
    biffAddf(MPR, "%s: trouble setting up cmap", __func__);
    airMopError(mop);
    return 1;
  }
  float  *idataF = (nrrdTypeFloat  == nin->type ? (float*)nin->data  : NULL);
  double *idataD = (nrrdTypeDouble == nin->type ? (double*)nin->data : NULL);
  for (uint pi=0; pi<cmap->num; pi++) {
    for (uint ci=0; ci<4; ci++) {
      cmap->data[ci + 4*pi] = (idataF
                               ? idataF[ci + 4*pi]
                               : idataD[ci + 4*pi]);
    }
  }
  airMopOkay(mop);
  return 0;
}

