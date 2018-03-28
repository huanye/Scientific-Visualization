#include "plt.h"

// sanity check on order of derivative
#define DMAX 10

/*
******** pltConvoEval: compute (or try to compute) convolution-based
** reconstruction of data with kern, or some derivative of the
** reconstruction, at world-space location xw.  This is the function that
** implements the abstraction of having a continuous function defined over
** some segment of world-space. If convolution can't be performed because
** some data values are missing, then record how many were missing.
**
** INPUT:
** xw: the world-space position at which to evaluate the convolution
**
** deriv: evaluate convolution with the deriv-th derivative of the given
** kernel kern.  For no differentiation, deriv=0. To measure a first
** derivative, deriv=1; second derivative, deriv=2, etc. The differentiation
** is in *world* space! It is your responsibility to account for difference
** between differentiation in index space (convolving with a derivative
** kernel) vs in world space (taking into account space between samples).
**
** data: the discretely sampled data to convolve with
**
** kern: the kernel to convolve with
**
** OUTPUT:
** result: The convolution result, when it can be computed, is stored in
** *result. Otherwise (when convolution can't be computed), *result should
** left unchanged.
**
** outside: If the convolution sum could not be computed (that is, because
** some terms (associated with the support of the kernel) required data
** values data->vv[i] for which i was outside the valid range of data indices
** [0,data->len-1]), then *result is not set, and *outside records the number
** missing data values (the number of required data indices that were outside
** [0,data->len-1]). Missing some of the required data values is NOT an
** actual error, in the sense of a non-zero return value and the subsequent
** error handling.  When the convolution result can be computed (all the data
** sample indices within the kernel support were valid data indices), then
** set *outside to 0.
**
** didx: **if non-NULL**, store here the kern->support data indices (ordered from
** low to high) required for the convolution sum, some of which may in fact
** be outside the valid index bounds of the data. Your convolution implementation
** should look at no more than these data values. Assume that a non-NULL didx
** has been allocated for at least kern->support values. **If didx is NULL**,
** that is the caller's way of saying "I don't care to learn this information".
**
** kpos: **if non-NULL**, store here the kern->support kernel evaluation
** locations that were used in the convolution sum, such that each term of
** the convolution sum is data->vv[didx[i]]*kern->eval(kpos[i]).  Assume that
** a non-NULL didx has been allocated for at least kern->support values.
** **If kpos is NULL**, the caller is saying "I don't care to learn this
** information".
**
** The point of didx and kpos is to provide a way for the graders to learn
** about your understanding of how to implement the convolution formula. When
** either is non-NULL, you must set all kern->support elements (with indices
** 0 through kern->support-1) of them. Your code MUST NOT assume, however,
** that didx and kpos are non-NULL. Regardless of getting NULL didx or kpos,
** you will still need to know the data indices and kernel positions in order
** to correctly compute compute the convolution, and the way that *result and
** *outside are set is unrelated to whether didx and/or kpos is NULL.
**
** Note that in practical use, this function might may be called multiple
** times, first with deriv=0, then deriv=1, ..., with the same xw, in which
** case there is certainly redundant computation being done. For the sake of
** this project, it is worth taking that efficiency hit in order to simplify
** the implementation.
**
** RETURNS:
** 0 if all is well, else
** 1 in case of error (using biff), because of missing pointers or
** nonfinite xw.  No additional error checking is required in your code.
*/
int
pltConvoEval(float *result, unsigned short *outside,
             int *didx, float *kpos, float xw, uint deriv,
             const pltData *data, const pltKernel *kern) {
  if (!(result && outside && data && kern)) {
    biffAddf(PLT, "%s: got NULL pointer (%p,%p,%p,%p)", __func__,
             (void*)result, (void*)outside, (void*)data, (void*)kern);
    return 1;
  }
  // it is ok for didx and kpos to be NULL
  if (!isfinite(xw)) {
    biffAddf(PLT, "%s: given world-space position %g not finite", __func__, xw);
    return 1;
  }
  if (deriv > DMAX) {
    biffAddf(PLT, "%s: deriv %u > reasonable limit %u", __func__, deriv, DMAX);
    return 1;
  }

  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  // L is the world space interval distance between two sample points
  int denom = data->center==pltCenterNode?(signed)(data->len)-1:(signed)data->len;
  float L = (data->max-data->min)/denom;
  float xi = pltWtoI(data, xw);
  // dealing with high order derivative
  const pltKernel * the_kern  = kern;
  for(int i=0;i<(signed)deriv;i++)
    the_kern = the_kern->deriv;
  int s = the_kern->support;
  // declare 3 varaibles taking different values 
  // which depends on if the support is even or odd.
  int n, upper, lower;
  if(s%2==0){
      n = (int)floor(xi);
      upper = s/2;
      lower = 1-s/2;
  }
  else{
      n = (int)floor(xi+0.5);
      upper = (s-1)/2;
      lower = (1-s)/2;
  }
  float alpha = xi-n;
  //initialize some variables for required results 
  *outside = 0;
  // the accummulated sum of the convolution,
  // which would be assigned to the *result if
  // no computed index is outside the data range. 
  float acc_sum = 0;
  // boolean flag recording if outside index exist.
  int no_outside = 1;
  for(int i=lower,j=0;i<=upper;i++,j++){
      if(didx)
          didx[j] = n+i;
      if(kpos)
          kpos[j] = alpha-i;  
      if(n+i<0 || n+i>(signed)(data->len)-1){
          *outside+=1;
          no_outside = 0;
          continue;
      }
      else
          acc_sum+=data->vv[n+i]*the_kern->eval(alpha-i);
  }
  // if all computed indices are in the data range,
  // compute and convert the result back to the world space
  if(no_outside){
      *result = acc_sum;    
      for(int i=0;i<(signed)deriv;i++)
          *result = *result/L;
   }
  // ( 40 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return 0;
}

int
pltConvoSample(pltData *odata,
               uint olen, float omin, float omax, pltCenter center,
               uint deriv,
               const pltData *idata, const pltKernel *kern) {
  if (!(odata && idata && kern)) {
    biffAddf(PLT, "%s: got NULL pointer (%p,%p,%p)", __func__,
             (void*)odata, (void*)idata, (void*)kern);
    return 1;
  }
  if (pltDataInit(odata, olen, omin, omax, center)) {
    biffAddf(PLT, "%s: problem initializing output data", __func__);
    return 1;
  }
  float result;
  unsigned short outside;
  if (pltConvoEval(&result, &outside, NULL, NULL,
                   (idata->max + idata->min)/2, deriv,
                   idata, kern)) {
    biffAddf(PLT, "%s: trial convolution failed", __func__);
    return 1;
  }

  /* here is where pltConvoEval is computed over all the sample positions,
     and where inability to do convolution (because of non-zero outside) is
     recorded as a NaN value, with outside as payload. */
  for (uint si=0; si<olen; si++) {
    float sw = pltItoW(odata, si);
    pltConvoEval(&result, &outside, NULL, NULL, sw, deriv, idata, kern);
    if (!outside) {
      odata->vv[si] = result;
    } else {
      odata->vv[si] = pltNan(outside);
    }
  }

  return 0;
}
