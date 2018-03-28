#include "vcr.h"
#include "vcrPrivate.h"

/*
** vcrConvoEval: perform (or try to) convolution-based reconstruction
** of input vector field ctx->vfl with ctx->kern at world-space location (xw, yw).
** If "getjac" is non-zero, then the first derivative of the field (the Jacobian)
** should be measured as well.  Both the reconstructed vector and Jacobian
** are saved within the context, in ctx->vec and ctx->jac, respectively.
**
** Converting (xw,yw) to index space (xi,yi) may produce a position
** for which the separable kernel support, centered at (xi,yi), is
** non-zero for data data indices that are out of bounds, in which
** case the convolution sum cannot be computed.  Note that the given
** code initializes the output to represent this case; when you can
** do the convolution, set ctx->inside to 1, and set ctx->vec (and
** possibly ctx->jac).
**
** This function does not really need error checking; since we're
** called in contexts that should do sufficient error checking for us
*/
void
vcrConvoEval(vcrCtx *ctx, real xw, real yw, int getjac) {
  USED(xw); USED(yw);
  // initialize output to represent failure to reconstruct
  ctx->inside = 0;
  real nn = vcrNan(0);
  V2_SET(ctx->vec, nn, nn);
  if (getjac) {
    M2_SET(ctx->jac, nn, nn, nn, nn);
  }
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real WPOS[3];
  real IPOS[3];
  WPOS[0] = xw;
  WPOS[1] = yw;
  WPOS[2] = 1;
  
  MV3_MUL(IPOS, ctx->WtoI, WPOS);
  real ix = IPOS[0];
  real iy = IPOS[1];
  int n1 = ctx->odd_support?floor(ix+0.5):floor(ix);
  int n2 = ctx->odd_support?floor(iy+0.5):floor(iy);
  real alpha1 = ix - n1;
  real alpha2 = iy - n2;

  for(int i = ctx->lower,j=0;i<=ctx->upper;i++,j++){
      ctx->kvalues1[j] = ctx->kern->eval(alpha1-i);
      ctx->kvalues2[j] = ctx->kern->eval(alpha2-i);
  }
  ctx->inside = 1;
  for(int i = ctx->lower;i<=ctx->upper;i++){
    if(i+n2<0||i+n2>(signed)(ctx->vfl->size[1])-1){
            ctx->inside=0;
      }
    if(i+n1<0||i+n1>(signed)(ctx->vfl->size[0])-1){
            ctx->inside=0;
      }
  }

  if(ctx->inside){
   real cov_sum0 = 0.0;
   real cov_sum1 = 0.0;
   for(int i2 = ctx->lower,j2=0;i2<=ctx->upper;i2++,j2++)
      for(int i1 = ctx->lower,j1=0;i1<=ctx->upper;i1++,j1++){
       
        int index_1D_0 = (n2+i2)*(ctx->vfl->size[0]*ctx->vfl->channel)+
                       (n1+i1)*ctx->vfl->channel+0;
        int index_1D_1 = index_1D_0+1;
        cov_sum0 += ctx->vfl->data.fl[index_1D_0]*ctx->kvalues1[j1]*ctx->kvalues2[j2];
        cov_sum1 += ctx->vfl->data.fl[index_1D_1]*ctx->kvalues1[j1]*ctx->kvalues2[j2];     
   }
   ctx->vec[0] = cov_sum0;
   ctx->vec[1] = cov_sum1;
  }
  if(getjac&&ctx->inside){
     
    for(int i = ctx->lower,j=0;i<=ctx->upper;i++,j++){
      ctx->the_kvalues1[j] = ctx->the_kern->eval(alpha1-i);
      ctx->the_kvalues2[j] = ctx->the_kern->eval(alpha2-i);
    }
     real cov_sum_deriv0 = 0;
     real cov_sum_deriv1 = 0;
     real cov_sum_deriv2 = 0;
     real cov_sum_deriv3 = 0;
     for(int i2 = ctx->lower,j2=0;i2<=ctx->upper;i2++,j2++){
       for(int i1 = ctx->lower,j1=0;i1<=ctx->upper;i1++,j1++){
        int the_index_1D_0 = (n2+i2)*(ctx->vfl->size[0]*ctx->vfl->channel)+
                            (n1+i1)*ctx->vfl->channel+0;
        int the_index_1D_1 = the_index_1D_0+1;
        cov_sum_deriv0 += ctx->vfl->data.fl[the_index_1D_0]*ctx->the_kvalues1[j1]*ctx->kvalues2[j2];
        cov_sum_deriv1 += ctx->vfl->data.fl[the_index_1D_0]*ctx->kvalues1[j1]*ctx->the_kvalues2[j2];
        cov_sum_deriv2 += ctx->vfl->data.fl[the_index_1D_1]*ctx->the_kvalues1[j1]*ctx->kvalues2[j2];
        cov_sum_deriv3 += ctx->vfl->data.fl[the_index_1D_1]*ctx->kvalues1[j1]*ctx->the_kvalues2[j2];

       }
      }

    // convert jacobian from index space to world space
    real g[2];
    g[0] = cov_sum_deriv0;
    g[1] = cov_sum_deriv1;
    real gw[2];
    MV2_MUL(gw, ctx->deriv_ItoW, g);
    ctx->jac[0] = gw[0];
    ctx->jac[1] = gw[1];
    g[0] = cov_sum_deriv2;
    g[1] = cov_sum_deriv3;
    MV2_MUL(gw, ctx->deriv_ItoW, g);
    ctx->jac[2] = gw[0];
    ctx->jac[3] = gw[1];
  }

  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (83L in ref)

  return;
}
