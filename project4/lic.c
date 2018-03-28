#include "vcr.h"
#include "vcrPrivate.h"

/* You may define here new functions or #define to be used in this
   file.  Function declarations here must begin with "static". */
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
static real lerp(real omin, real omax, real imin, real xx, real imax) {
  real alpha = (xx-imin)/(imax-imin);
  return (1-alpha)*omin+alpha*omax;
}
static int min(int a, int b){return a<b?a:b;}
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (5L in ref)

/*
** vcrPropCalc computes from the given vector vec and the
** Jacobian jac the property requested.
*/
real
vcrPropCalc(int prop, const real vec[2], const real jac[4]) {
  USED(jac);
  real ret=vcrNan(0);
  if (!airEnumValCheck(vcrProp_ae, prop)) {
    switch (prop) {
    case vcrPropOne:
      ret = 1;
      break;
    case vcrPropLength:
      ret = V2_LEN(vec);
      break;
    case vcrPropVort: // vorticity
      // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
      ret = jac[2]-jac[1];
      // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (1L in ref)
      break;
    case vcrPropDiverg: // divergence; trace of Jacobian
      // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
      ret = jac[0] + jac[3]; 
      // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (1L in ref)
      break;
    case vcrPropJacFrob: // Frobenius norm of Jacobian
      // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
      ret = sqrt(jac[0]*jac[0]+jac[1]*jac[1]+jac[2]*jac[2]+jac[3]*jac[3]);
      // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (2L in ref)
      break;
    case vcrPropJacDet: // determinant of Jacobian
      // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
       ret = jac[0]*jac[3] - jac[1]*jac[2];
      // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (1L in ref)
      break;
    }
  }
  return ret;
}

/*
******** vcrLICEval: evaluate LIC at one seed point
**
** INPUT:
** - vcrSline *sln: to contain the streamline vertices
** - wx,wy: seedpoint for streamline integration
** - wkern: kernel for weighting noise samples along streamline
** - rnd: the noise image
** - rndLinterp: if non-zero and if you're in 33710: do linear interpolation
**   into the noise texture. Else, use nearest neighbor interpolation
** - hh: integration step size
** - normalize: integrate not through the vector field V=ctx->vfl, but
**   through the normalized V/|V|
** - intg: how to integrate the streamline
** - ctx: the computational context, used for convolution
**
** OUTPUT:
** - the LIC result is saved in *result.
** - The computed streamline is saved in sln as part of computing LIC,
**   but if the streamline wasn't seeded inside the vector field then
**   sln->seedInside will be 0.  The streamline geometry in sln->pos
**   is used in this function, but will likely be over-written with
**   later calls to this same function.
**
** Returns 1 in case of error (described w/ biff), else 0 if all is well.
*/
int
vcrLICEval(real *result, vcrSline *sln, real wx, real wy,
           const vcrKernel *wkern, const vcrImage *rmg, int rmgLinterp,
           uint halfLen, real hh, int normalize, int intg,
           vcrCtx *ctx) {
  USED(wkern); USED(rmg); USED(rmgLinterp);
  /* error checking we assume has already been done (e.g. by vcrLIC):
     -  result, sln, wkern, rmg, ctx non-NULL
     -  wkern != vcrKernelZero
     -  rmg image is scalar (1==rmg->channel) float (vcrTypeFloat==rmg->dtype)
     -  hh > 0
     -  intg is a valid value in the vcrIntg enum */
  if (vcrSlineTrace(sln, wx, wy, halfLen, hh, normalize, intg, ctx)) {
    biffAddf(VCR, "%s: trouble computing streamline", __func__);
    return 1;
  }
  *result = vcrNan(0);
  if (!sln->seedInside) {
    return 0;
  }
  /* Your code below should conform to the instructions in vcrLIC on
     how to evaluate LIC. There is no need for further error checking
     (i.e. the only reason vcrLICEval returns non-zero is if
     vcrSlineTrace returns non-zero) */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real WtoI_noise[9]; 
  real tmp;
  M3_INVERSE(WtoI_noise, rmg->ItoW, tmp);
  real tsum =0.0;
  real wsum = 0.0;
  int incount = 0;
  uint posIdx;
  for (posIdx = halfLen - sln->backNum;
        posIdx <= halfLen + sln->forwNum;
           posIdx++){
  real wpos[3];
  real ipos[3];
  V3_SET(wpos,sln->pos[2*posIdx],sln->pos[2*posIdx+1],1);
  MV3_MUL(ipos, WtoI_noise, wpos);
  real ix = ipos[0];
  real iy = ipos[1];
  real rr;
  if(ix>=0 && ix<=rmg->size[0]-1 && iy>=0 && iy<=rmg->size[1]-1){
    if(ix==floor(ix)&&iy==floor(iy)){
      rr = rmg->data.fl[(int)iy*(rmg->size[0])+(int)ix];
    }
   else if(rmgLinterp){
    // bilinear lerp
    real alphax = ix-floor(ix);
    real alphay = iy-floor(iy);
    //first lerp on along x axis two times
    real lerpx1 = lerp(rmg->data.fl[(int)floor(iy)*(rmg->size[0])+(int)floor(ix)],
         rmg->data.fl[(int)floor(iy)*(rmg->size[0])+min((int)(floor(ix)+1),rmg->size[0])],
         0,alphax,1);
    real lerpx2 = lerp(rmg->data.fl[(min((int)(floor(iy)+1),rmg->size[1]))*(rmg->size[0])+(int)floor(ix)],
         rmg->data.fl[min((int)(floor(iy)+1),rmg->size[1])*(rmg->size[0])+min((int)(floor(ix)+1),rmg->size[0])],
         0,alphax,1);
    // then lerp along y axis
    rr = lerp(lerpx1,lerpx2,0,alphay,1);
    }
   else{
    if(ix-floor(ix)<0.5)
      ix = floor(ix);
    else
      ix = floor(ix+1);

    if(iy-floor(iy)<0.5)
      iy = floor(iy);
    else
      iy = floor(iy+1);
    rr = rmg->data.fl[(int)iy*(rmg->size[0])+(int)ix];
   }
  // lerp for kp
  real kp = lerp(-0.4999*wkern->support,0.4999*wkern->support,
                 0,posIdx,2*halfLen);
  real ww = wkern->eval(kp);
  wsum += ww*ww;
  tsum += ww*rr;
  incount += 1;
  }
}
if (wsum) 
     tsum /= sqrt(wsum);

tsum /= sqrt(incount);
result[0] = tsum;
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (88L in ref)
  return 0;
}

/*
** vcrLIC: compute Line Integral Convolution (LIC) of vector field in given
** vcrCtx "ctx", by starting one streamline for each pixel in given noise
** texture image "rmg", saving the LIC result in "lmg" (Lic iMaGe).
** Depending on the value "prop", different things about the field or its
** Jacobian are saved into the corresponding pixel of "pmg" (Property
** iMaGe). The samples of the noise texture along the streamline are found by
** nearest neighbor interpolation (except if you're a 33710 student and
** rmgLinterp is non-zero), and then weighted by kernel "wkern" in the
** convolution sum (along the streamline).  Parameters "halfLen", "hh",
** "normalize", and "intg" passed directly to vcrSlineTrace.
**
** Note that the only kind of noise texture image you need need to handle is a
** scalar (not a multi-component) image of floats, and the output images lmg
** and pmg are also created as scalar images of floats.  For all of these,
** you get at the pixel values via vcrImage->data.fl
**
** The caller should have created "rmg" so that its world-space domain fits
** entirely within the domain of the vector field, so that even at the corners
** of "rmg", streamlines have somewhere to go.
*/
int
vcrLIC(vcrImage *lmg, vcrImage *pmg, int prop,
       const vcrKernel *wkern, const vcrImage *rmg, int rmgLinterp,
       uint halfLen, real hh, int normalize, int intg,
       vcrCtx *ctx) {
  USED(rmgLinterp);
  if (!(lmg && pmg && wkern && rmg && ctx)) {
    biffAddf(VCR, "%s: got NULL pointer (%p,%p,%p,%p,%p)", __func__,
             (void*)lmg, (void*)pmg, (void*)wkern, (void*)rmg, (void*)ctx);
    return 1;
  }
  if (vcrKernelZero == wkern) {
    biffAddf(VCR, "%s: can't use %s kernel for LIC weighting",
             __func__, vcrKernelZero->name);
    return 1;
  }
  if (!( 1 == rmg->channel && vcrTypeFloat == rmg->dtype )) {
    biffAddf(VCR, "%s: need 1-channel (not %u) image of %s (not %s) for noise",
             __func__, rmg->channel, airEnumStr(vcrType_ae, vcrTypeFloat),
             airEnumStr(vcrType_ae, rmg->dtype));
    return 1;
  }
  if (!(halfLen < 1024)) {
    biffAddf(VCR, "%s: halfLen %u seems unreasonable", __func__, halfLen);
    return 1;
  }
  if (!(hh > 0)) {
    biffAddf(VCR, "%s: stepsize (hh) %g not > 0", __func__, hh);
    return 1;
  }
  if (airEnumValCheck(vcrIntg_ae, intg)) {
    biffAddf(VCR, "%s: %d not a valid %s", __func__, intg,
             vcrIntg_ae->name);
    return 1;
  }
  if (airEnumValCheck(vcrProp_ae, prop)) {
    biffAddf(VCR, "%s: %d not a valid %s", __func__, prop,
             vcrProp_ae->name);
    return 1;
  }

  airArray *mop = airMopNew();
  vcrSline *sln = vcrSlineNew();
  airMopAdd(mop, sln, (airMopper)vcrSlineNix, airMopAlways);
  if (vcrSlineAlloc(sln, halfLen)
      || vcrImageAlloc(lmg, 1, rmg->size[0], rmg->size[1], vcrTypeFloat)
      || vcrImageAlloc(pmg, 1, rmg->size[0], rmg->size[1], vcrTypeFloat)) {
    biffAddf(VCR, "%s: trouble allocating output", __func__);
    airMopError(mop);
    return 1;
  }
  M3_COPY(lmg->ItoW, rmg->ItoW);
  M3_COPY(pmg->ItoW, rmg->ItoW);
  {
    /* do trial streamline, to minimize risk from student code not
       doing its own error checking. Using braces to make ipos and
       wpos local to this block. */
    real ipos[3], wpos[3];
    V3_SET(ipos, (rmg->size[0]-1)/2, (rmg->size[1]-1)/2, 1);
    MV3_MUL(wpos, rmg->ItoW, ipos);
    if (vcrSlineTrace(sln, wpos[0], wpos[1], halfLen, hh,
                      normalize, intg, ctx)) {
      biffAddf(VCR, "%s: trouble doing test trace", __func__);
      airMopError(mop);
      return 1;
    }
  }

  /* Todo:
  **
  ** Learn inverse of rmg->ItoW, so that you can query the noise image at each
  ** vertex along each streamline.
  **
  ** For each pixel of noise image "rmg":
  ** -- Find world-space location of pixel, call it wpos
  ** -- Compute streamline starting there by calling
  **    vcrSlineTrace(sln, wpos[0], wpos[1], halfLen, hh, normalize, intg, ctx);
  ** -- If (!sln->seedInside), set the output pixel in both lmg and pmg to 0
  ** -- Else, save vcrPropCalc(prop, sln->vecSeed, sln->jacSeed) to pmg, and
  **    do LIC for this pixel, with the following.  Initialize:
  **       real tsum=0, wsum=0;
  **       int incount=0;
  **    and then for each streamline vertex, indexed by posIdx:
  **    (see vectr_sline.c for an example of such a for-loop; note that
  **     posIdx varies all the way from 0 to 2*halfLen if neither end of
  **     the streamline was stopped by the vector field boundary):
  **    -- convert the world-space streamline vertex position to "rmg"'s
  **       (continuous) index-space (i,j)
  **    -- If (i,j) is inside "rmg" (that is, you you can do the required
  **       noise texture interpolation at (i,j)), then
  **          rr = rmg at (i,j), via nearest-neighbor or linear interpolation,
  **               (mimicking the behavior of separable convolution with
  **               the box or tent kernels, respectively)
  **          ww = wkern->eval(kp)
  **          wsum += ww*ww
  **          tsum += ww*rr
  **          incount += 1
  **       where kp varies from -0.4999*wkern->support to
  **       0.4999*wkern->support as posIdx varies from 0 to 2*halfLen
  **       (using 0.4999 instead of 0.5 to ensure that with wkern=box,
  **        we really get the first and last ww to be 1.0)
  **    -- Else (i,j NOT inside rmg): wsum, tsum, incount are unchanged
  ** -- Normalize the texture convolution sum tsum according to the
  **    "wkern" kernel weights:
  **       if (wsum) tsum /= sqrt(wsum);
  **    This is taking the RMS (root mean square) of the kernel weights,
  **    assuming that they were non-zero.
  ** -- Then normalize the texture convolution sum tsum according to the
  **    number of samples: the standard deviation of the sum of N samples
  **    of noise (drawn from a normal distribution with stdv==1) is
  **    sqrt(N) (https://en.wikipedia.org/wiki/Signal_averaging),
  **    so we can put convolution output into a normalized range by:
  **       tsum /= sqrt(incount);
  ** -- Save final tsum to lmg output image.
  **
  ** This algorithm is slightly different (mainly in how results are
  ** normalized) than what was presented in class.
  **
  ** No further error checking (biff usage) is needed.
  */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  

  for(uint i2=0;i2<rmg->size[1];i2++)
    for(uint i1=0;i1<rmg->size[0];i1++){
        real ipos[3], wpos[3];
        V3_SET(ipos, i1,i2, 1);
        MV3_MUL(wpos, rmg->ItoW, ipos);
        real result[1];
        vcrLICEval(result, sln, wpos[0], wpos[1],wkern,rmg,rmgLinterp,halfLen,
          hh,normalize,intg,ctx);
        if (!sln->seedInside){
          pmg->data.fl[i2*(rmg->size[0])+i1] = 0;
          lmg->data.fl[i2*(rmg->size[0])+i1] = 0;
        }
        else{
          pmg->data.fl[i2*(rmg->size[0])+i1] = 
            vcrPropCalc(prop, sln->vecSeed, sln->jacSeed);
          lmg->data.fl[i2*(rmg->size[0])+i1] = result[0];
          
        }

  }

  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (22L in ref)

  airMopOkay(mop);
  return 0;
}
