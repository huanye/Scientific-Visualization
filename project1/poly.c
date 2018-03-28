#include "plt.h"

/*
******** pltPolyEval: evaluates at xx a polynomial with coefficients pc
**
** INPUT:
** xx: position at which to evaluate polynomial
**
** pc, pclen: pclen is the length of the array pc; valid values are
** pc[0], pc[1], ... pc[pclen-1].  The polynomial to evaluate is:
** The pclen coefficients of pc are stored in increasing degree; the
** polynomial to evaluate is
** p(x) = pc[0] + x*pc[1] + (x^2)*pc[2] ... + (x^(pcn-1))*pc[pcn-1]
**
** OUTPUT:
** Returns p(x). There is no error checking. Assume pc is non-NULL.
**
** NOTE: Use Horner's method to evaluate the polynomial; you will
** lose points if you don't. Certainly, do not use pow() for
** exponentiation; it is not needed.
*/
float
pltPolyEval(const float *pc, uint pclen, float xx) {
  float ret=0;
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  for(int i=(signed)(pclen-1);i>=0;i--){
     ret = pc[i]+ret*xx;
  }
  // ( 3 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return ret;
}

/*
******** pltPolySample: evaluate polynomial p (given by pc, pclen) with len
** samples in the interval [imin, imax], either cell-centered or
** node-centered depending on center.
**
** INPUT:
** len: the number of samples at which to evaluate the polynomial
**
** imin, imax: len samples are placed inside interval [imin, imax]
** (with either cell- or node-centered sampling), and the polynomial
** is evaluated at these samples.
**
** center: either pltCenterCell or pltCenterNode, to request either cell-
** or node-centered sampling.
**
** pc, pclen: defines a polynomial p(x), same as with pltPolyEval
**
** OUTPUT:
** data: the given code below will initialize this container to store len
** samples.  The results of sampling p(x) (via repeated calls to
** pltPolyEval) should be stored in data->vv[]. data->vv[0] should hold
** the results of evaluating p(x) for x nearest imin, regardless of
** whether imin is less than or greater than imax.
**
** RETURNS:
** 0 if all is well, 1 if there is an error, which is detected with the
** given code.  No additional error handling should be in the student code.
*/
int
pltPolySample(pltData *data,
              uint len, float imin, float imax, pltCenter center,
              const float *pc, uint pclen) {
  if (!(data && pc)) {
    biffAddf(PLT, "%s: got NULL pointer (%p,%p)", __func__,
             (void*)data, (void*)pc);
    return 1;
  }
  if (!pclen) {
    biffAddf(PLT, "%s: got zero polynomial coeffs?", __func__);
    return 1;
  }
  if (pltDataInit(data, len, imin, imax, center)) {
    biffAddf(PLT, "%s: problem initializing data", __func__);
    return 1;
  }

  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  int i;
  for(i=0;i<(signed)len;i++){
      data->vv[i] = pltPolyEval(pc, pclen, pltItoW(data, i));
  }
  
  // ( 4 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

  return 0;
}
