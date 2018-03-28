/*
** rnd and rendr: SciVis-2018 Project 3
** Copyright (C)  2018 University of Chicago. All rights reserved.
**
** This file is distributed to students in the 2018 CMSC 23710/33710
** ("SciVis") class, and no one else, for use in the SciVis class. This
** file is not open-source licensed, nor licensed for any other distribution.
** You should not allow this file to be copied or downloaded by someone who
** is not a 2018 SciVis student, as that would violate the copyright
** held by University of Chicago.
*/

#include "rnd.h"
#include "rndPrivate.h"

rndTxf *
rndTxfNew(void) {
  rndTxf *txf = MALLOC(1, rndTxf); assert(txf);
  txf->space = rndSpaceUnknown;
  txf->num = 0;
  txf->data = NULL;
  return txf;
}

rndTxf *
rndTxfNix(rndTxf *txf) {
  if (txf) {
    free(txf->data);
    free(txf);
  }
  return NULL;
}

static uint
colorSize0(rndSpace space) {
  uint ret=0;
  if (rndSpaceRGB == space || rndSpaceHSV == space) {
    ret = 4;
  } else if (rndSpaceAlpha == space) {
    ret = 2;
  }
  return ret;
}

int
rndTxfAlloc(rndTxf *txf, rndSpace space, uint pointNum) {
  if (!txf) {
    biffAddf(RND, "%s: got NULL pointer", __func__);
    return 1;
  }
  uint size0 = colorSize0(space);
  if (!size0) {
    biffAddf(RND, "%s: space %d not valid", __func__, space);
    return 1;
  }
  if (!pointNum) {
    biffAddf(RND, "%s: need non-zero number of points", __func__);
    return 1;
  }
  free(txf->data);
  txf->data = MALLOC(size0*pointNum, real);
  if (!(txf->data)) {
    biffAddf(RND, "%s: couldn't alloc data for %u real %u-vecs",
             __func__, pointNum, size0);
    return 1;
  }
  txf->space = space;
  txf->num = pointNum;
  return 0;
}

int
rndTxfSave(const char *fname, const rndTxf *txf) {
  if (!(fname && txf && txf->data)) {
    biffAddf(RND, "%s: got NULL pointer (%p,%p,%p)", __func__,
             (void*)fname, (void*)txf, (void*)txf->data);
    return 1;
  }
  uint size0 = colorSize0(txf->space);
  if (!size0) {
    biffAddf(RND, "%s: space %d not valid", __func__, txf->space);
    return 1;
  }
  airArray *mop = airMopNew();
  // (can't nrrdWrap because its void *data isn't const)
  Nrrd *nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdAlloc_va(nout, nrrdTypeReal, 2,
                   AIR_CAST(size_t, size0),
                   AIR_CAST(size_t, txf->num))) {
    biffMovef(RND, NRRD, "%s: trouble allocating output nrrd", __func__);
    airMopError(mop);
    return 1;
  }
  memcpy(nout->data, txf->data, sizeof(real)*size0*txf->num);
  if (nrrd1DIrregMapCheck(nout)) {
    biffMovef(RND, NRRD, "%s: control point data not ok", __func__);
    airMopError(mop);
    return 1;
  }
  if (nrrdKeyValueAdd(nout, rndKeySpace,
                      airEnumStr(rndSpace_ae, txf->space))) {
    biffAddf(RND, "%s: trouble setting key/value pair", __func__);
    airMopError(mop);
    return 1;
  }
  if (nrrdSave(fname, nout, NULL)) {
    biffMovef(RND, NRRD, "%s: trouble saving", __func__);
    airMopError(mop);
    return 1;
  }
  airMopOkay(mop);
  return 0;
}

/* static because even if this is generally useful;
   it is only needed inside this file */
static int
rndNrrdTxfCheck(const Nrrd *nin) {
  if (nrrdCheck(nin)
      || nrrd1DIrregMapCheck(nin)) {
    biffMovef(RND, NRRD, "%s: control point array not ok", __func__);
    return 1;
  }
  if (nrrdTypeReal != nin->type) {
    biffAddf(RND, "%s: need %s type data, not %s", __func__,
             airEnumStr(nrrdType, nrrdTypeReal),
             airEnumStr(nrrdType, nin->type));
    return 1;
  }
  if (!isfinite(((real*)(nin->data))[0])) {
    biffAddf(RND, "%s: first control point location not set", __func__);
    return 1;
  }
  /* can be done by caller: checking the "space:=" key, and making
     sure that nin->axis[0].size agrees with that */
  return 0;
}

int
rndTxfLoad(rndTxf *txf, const char *fname) {
  if (!(txf && fname)) {
    biffAddf(RND, "%s: got NULL pointer (%p,%p)", __func__,
             (void*)txf, (void*)fname);
    return 1;
  }
  airArray *mop = airMopNew();
  Nrrd *_nin = nrrdNew();
  airMopAdd(mop, _nin, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdLoad(_nin, fname, NULL)) {
    biffMovef(RND, NRRD, "%s: trouble reading file", __func__);
    airMopError(mop);
    return 1;
  }
  Nrrd *nin = nrrdNew();
  airMopAdd(mop, nin, (airMopper)nrrdNuke, airMopAlways);
  if (rndArrayToReal(nin, _nin)
      || rndNrrdTxfCheck(nin)) {
    biffAddf(RND, "%s: given array doesn't look like a transfer function",
             __func__);
    airMopError(mop);
    return 1;
  }
  char *spcS = nrrdKeyValueGet(nin, rndKeySpace);
  airMopAdd(mop, spcS, airFree, airMopAlways);
  if (!spcS) {
    biffAddf(RND, "%s: color/alpha space not set in \"%s:=\" key/value pair",
             __func__, rndKeySpace);
    airMopError(mop);
    return 1;
  }
  rndSpace space = airEnumVal(rndSpace_ae, spcS);
  if (rndSpaceUnknown == space) {
    biffAddf(RND, "%s: \"%s:=%s\" key/value pair not recognized",
             __func__, rndKeySpace, spcS);
    airMopError(mop);
    return 1;
  }
  uint size0 = colorSize0(space);
  if (size0 != nin->axis[0].size) {
    biffAddf(RND, "%s: axis[0].size should be %u, not %u", __func__,
             size0, (uint)nin->axis[0].size);
    return 1;
  }
  if (rndTxfAlloc(txf, space, (uint)nin->axis[1].size)) {
    biffAddf(RND, "%s: trouble allocating txf", __func__);
    airMopError(mop);
    return 1;
  }
  memcpy(txf->data, nin->data, sizeof(real)*size0*txf->num);
  airMopOkay(mop);
  return 0;
}

/* You may define here new functions or #define to be used in this
   file.  Function declarations here must begin with "static" */
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
// ( 8 lines in reference implementation )
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

/*
******** rndTxfEval: evaluate given transfer function txf (maps to either
** colors or opacities) at given value vv, and save result to out.
**
** INPUT: transfer function txf maps to colors if txf->space is rndSpaceRGB
** or rndSpaceHSV; it maps to opacities if txf->space is rndSpaceAlpha (note
** how colorSize0 is used below; this determines how many values are on the
** faster axis of the 2D array inside the transfer function).
**
** OUTPUT: In the case of a color transfer function, your code should save
** R,G,B values (converting from HSV if needed) to out[0], out[1],
** out[2]. With an opacity transfer function, save opacity to *out.  Like
** mprCmapEval in project 2, rndTxfEval should clamp the given value to the
** transfer function domain.
**
** For this project, all students (not just 33710 students) have to correctly
** handle hue interpolation: hue H in [0,1] is interpreted as angle 2*pi*H,
** on the circle, and the interpolation of H0 and H1 should traverse the
** shorter of the two possible paths on the circle between angles H0 and H1.
** You can assume that all hues are in the range [0,1].
**
** No additional error checking is needed; this will only be called in settings
** that have already done the required error checking.
*/
void
rndTxfEval(real *out, const rndTxf *txf, real vv) {
  uint size0 = colorSize0(txf->space);
  // can assume that 2==size0 (opacity txf) or 4==size0 (color txf)
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  uint index = 0;
  real alpha = 0.0;
  real lmin = txf->data[0];
  real lmax = txf->data[0 + size0*(txf->num-1)]; 
  if(vv>=lmax){
    index = txf->num-1;
    alpha = 0.0;
  }
  else if(vv<=lmin){
     index = 0;
     alpha = 0.0;
  }
  else{
   for(int i=0;i<(signed)txf->num-1;i++)
    if(txf->data[size0*(i+1)]>vv){
          index = i; 
          alpha = (vv-txf->data[size0*i])/(txf->data[size0*(i+1)]-txf->data[size0*i]);
          break;
   }
  }
  if(txf->space == rndSpaceAlpha){
    out[0] = (1-alpha)*txf->data[2*index+1]+alpha*txf->data[1+2*(index+1)];
    return; 
  }
  real color[3];
  color[1] = (1-alpha)*txf->data[4*index+2]+alpha*txf->data[2+4*(index+1)];
  color[2] = (1-alpha)*txf->data[4*index+3]+alpha*txf->data[3+4*(index+1)];
  if(txf->space == rndSpaceHSV){
    real hue0 = txf->data[4*index+1];
    real hue1 = txf->data[1+4*(index+1)];
    if(hue0>=hue1 && hue0 - hue1>0.5)
         hue1+=1.0;
    else if(hue0<hue1 && hue1 - hue0>0.5)
         hue0+=1.0;
    color[0] = (1-alpha)*hue0+alpha*hue1;
    if(color[0]>=1.0) 
       color[0]-=1.0;
    rndRGBfromHSV(out,color);
  }
  else if(txf->space == rndSpaceRGB){
  color[0] = (1-alpha)*txf->data[4*index+1]+alpha*txf->data[1+4*(index+1)];
  V3_COPY(out,color);
  }
  // ( 40 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return;
}

/*
** rndTxfLutGenerate: generate a lookup table (or "LUT") to that samples the
** transfer functions for color ("ctxf") and opacity ("atxf") over the value
** space ["vmin", "vmax"], with "len" samples. NOTE that either ctxf, or
** atxf, or both!, can be NULL, which is the way of saying "default": if ctxf
** is NULL, the rgb color is always (1,1,1); if atxf is NULL, the opacity is
** always 1.
*/
int
rndTxfLutGenerate(Nrrd *nlut,
                  const rndTxf *ctxf, int rescaleC,
                  const rndTxf *atxf, int rescaleA,
                  uint len, real vmin, real vmax) {
  if (!nlut) {
    biffAddf(RND, "%s: got NULL pointer", __func__);
    return 1;
  }
  if (!( ctxf || atxf )) {
    fprintf(stderr, "%s: Warning: got neither color or opacity txf\n", __func__);
    fprintf(stderr, "%s: LUT will have constant r,g,b,a=1,1,1,1\n", __func__);
  }
  if (ctxf) {
    if (!( rndSpaceRGB == ctxf->space || rndSpaceHSV == ctxf->space )) {
      biffAddf(RND, "%s: expecting ctxf with space %s or %s (not %s)",
               __func__, airEnumStr(rndSpace_ae, rndSpaceRGB),
               airEnumStr(rndSpace_ae, rndSpaceHSV),
               airEnumStr(rndSpace_ae, ctxf->space));
      return 1;
    }
  }
  if (atxf) {
    if (rndSpaceAlpha != atxf->space) {
      biffAddf(RND, "%s: expecting atxf with space %s (not %s)", __func__,
               airEnumStr(rndSpace_ae, rndSpaceAlpha),
               airEnumStr(rndSpace_ae, atxf->space));
      return 1;
    }
  }
  if (!len) {
    biffAddf(RND, "%s: need a postive lut length", __func__);
    return 1;
  }
  if (!( isfinite(vmin) && isfinite(vmax) && vmin < vmax )) {
    biffAddf(RND, "%s: need finite vmin < vmax (not " RCS " " RCS ")",
             __func__, vmin, vmax);
    return 1;
  }
  if (nrrdMaybeAlloc_va(nlut, nrrdTypeReal,
                        (size_t)2, (size_t)4, (size_t)len)) {
    biffMovef(RND, NRRD, "%s: output allocation problem", __func__);
    return 1;
  }
  nlut->axis[0].kind = nrrdKindRGBAColor;
  nlut->axis[1].min = vmin;
  nlut->axis[1].max = vmax;
  nlut->axis[1].center = nrrdCenterCell;
  real cvmin, cvmax;
  if (ctxf) {
    if (rescaleC) {
      cvmin = ctxf->data[0];
      cvmax = ctxf->data[0 + 4*(ctxf->num-1)];
    } else {
      cvmin = vmin;
      cvmax = vmax;
    }
  } else {
    if (rescaleC) {
      biffAddf(RND, "%s: can't rescale to domain of non-existent color txf", __func__);
      return 1;
    }
    cvmin = cvmax = rndNan(0);
  }
  real avmin, avmax;
  if (atxf) {
    if (rescaleA) {
      avmin = atxf->data[0];
      avmax = atxf->data[0 + 2*(atxf->num-1)];
    } else {
      avmin = vmin;
      avmax = vmax;
    }
  } else {
    if (rescaleA) {
      biffAddf(RND, "%s: can't rescale to domain of non-existent alpha txf", __func__);
      return 1;
    }
    avmin = avmax = rndNan(0);
  }
  real *rgba = (real*)(nlut->data);
  for (uint ii=0; ii<len; ii++) {
    real aa = ((real)ii + 0.5)/((real)len);
    real rgb[3];
    if (ctxf) {
      rndTxfEval(rgb, ctxf, (1-aa)*cvmin + aa*cvmax);
    } else {
      V3_SET(rgb, 1, 1, 1);
    }
    real opac;
    if (atxf) {
      rndTxfEval(&opac, atxf, (1-aa)*avmin + aa*avmax);
    } else {
      opac = 1;
    }
    V4_SET(rgba + 4*ii, rgb[0], rgb[1], rgb[2], opac);
  }
  return 0;
}

/*
** rndTxfLevoy: evaluate a list, length "num", of Levoy fuzzy isosurface
** functions; each one defined by three parameters on the faster axis of the
** given "vra" array. The three parameters are described in rnd.h, in the
** section about the "levoy" struct inside the rndCtx. The given (data value,
** gradient magnitude) coordinates ("val","gradmag") may overlap the support
** of multiple fuzzy isosurfaces defined this way, so this function should
** accumulate the opacities assigned to ("val","gradmag") (by each of the
** fuzzy isosurfaces specified) by *multiplying* *their* *transparencies*.
** If gradmag == 0, then the opacity for the isosurface defined by (V,R,A)
** should be A if val==V, else 0.
*/
real
rndTxfLevoy(const real *vra, uint num,
            real val, real gradmag) {
  real alpha=0;
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real trans_accum = 1;
  for(int i=0;i<(signed)num;i++){
      real a_v = vra[i*3+2];
      real r = vra[i*3+1];
      real f_v = vra[i*3];
      real diminish;
      if(gradmag == 0){
         if(val == f_v)
          diminish = 1;
        else
          diminish = 0;
      }
      else if(fabsf(f_v - val)<r*gradmag)
        diminish = 1-fabsf(f_v - val)/(r*gradmag);
      else
        diminish = 0;
      trans_accum*=(1-a_v*diminish);
  }
  alpha = 1-trans_accum;
  // ( 16 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return alpha;
}
