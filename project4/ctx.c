#include "vcr.h"
#include "vcrPrivate.h"

/* You may (if you want) define new functions or #define macros here;
   but our intent is only that they would be used in this file.  Thus,
   prefix any function declarations here with "static" */
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (0L in ref)

/*
** vcrCtxNew creates a context for doing convolution (and streamlines, and
** LIC) on given vector field "vfl" with kernel "kern". See comments within.
*/
vcrCtx *
vcrCtxNew(const vcrImage *vfl, const vcrKernel *kern) {
  if (!(vfl && kern)) {
    biffAddf(VCR, "%s: got NULL pointer (%p,%p)", __func__,
	     (void*)vfl, (void*)kern);
    return NULL;
  }
  if (2 != vfl->channel) {
    biffAddf(VCR, "%s: only works on vector (not %u-channel) images",
             __func__, vfl->channel);
    return NULL;
  }
  if (vcrTypeFloat != vfl->dtype) {
    biffAddf(VCR, "%s: only works on %s (not %s) images", __func__,
             airEnumStr(vcrType_ae, vcrTypeFloat),
             airEnumStr(vcrType_ae, vfl->dtype));
    return NULL;
  }

  vcrCtx *ctx = MALLOC(1, vcrCtx);
  if (!ctx) {
    biffAddf(VCR, "%s: malloc(vcrCtx) failure", __func__);
    return NULL;
  }
  ctx->vfl = vfl;
  ctx->kern = kern;

  /* Here, set up any additional resources (setting variables or allocating
     arrays) that will be used in vcrConvoEval() WHETHER OR NOT the Jacobian
     (the first derivative) is required; that will be known only by the
     "getjac" argument to vcrConvoEval. Be sure to free any
     dynamically-allocated memory in vcrCtxNix().  Also, note from the code
     above that you only have to worry about doing convolution in 2-channel
     images of float; that is the only kind of data accepted into a vcrCtx. */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  ctx->WtoI = MALLOC(9,real);
  real tmp;
  M3_INVERSE(ctx->WtoI, ctx->vfl->ItoW, tmp);
  int support = (signed)(ctx->kern->support);
  ctx->odd_support = (support%2);
  ctx->upper = ctx->odd_support?(support-1)/2:support/2;
  ctx->lower = ctx->odd_support?(1-support)/2:1-support/2;
  ctx->deriv_ItoW = MALLOC(4,real);
  M23_INVERSE_TRANSPOSE(ctx->deriv_ItoW, ctx->vfl->ItoW, tmp);
  ctx->kvalues1 = MALLOC(ctx->kern->support+1,real);
  ctx->kvalues2 = MALLOC(ctx->kern->support+1,real);
  ctx->the_kern = ctx->kern->deriv;
  ctx->the_kvalues1 = MALLOC(ctx->the_kern->support+1,real);
  ctx->the_kvalues2 = MALLOC(ctx->the_kern->support+1,real);
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (26L in ref)
  return ctx;
}

vcrCtx *
vcrCtxNix(vcrCtx *ctx) {

  if (ctx) {
    // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
    free(ctx->WtoI);
    free(ctx->deriv_ItoW);
    free(ctx->kvalues1);
    free(ctx->kvalues2);
    free(ctx->the_kvalues1);
    free(ctx->the_kvalues2);
    // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (6L in ref)
    free(ctx);
  }
  return NULL;
}
