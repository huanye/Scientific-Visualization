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
******** mprConvoEval: perform (or try to perform) convolution-based
** reconstruction of input image ctx->image with kernel ctx->kern (and maybe
** also with its first derivative) at world-space location (xw, yw).  The
** convolution is separable: k(x,y) = k(x)*k(y). If convolution can't be
** performed because some data values are missing, then indicate how many were
** missing.
**
******** INPUT:
**
** mprContext *ctx: contains all the state required for doing the
** convolution, in particular the image ctx->image, the kernel ctx->kern,
** ctx->mode which indicates what kind of computation will be required, and
** anything you set or allocated in mprCtxNew as a function of image, kern,
** and mode. The ctx->cmap is not relevant for this function.
**
** THE ONLY KIND OF mprImage YOU SHOULD HANDLE HERE IS SCALAR FLOATS.
** That is, you should assume
** 1 == ctx->image->channel, and
** mprTypeFloat == ctx->image->dtype, and so the image data values are in
** ctx->image->data.fl
** Also, even if we are compiled with make CFLAGS="-DMPR_REAL_IS_DOUBLE=1"
** (so that mprRealIsDouble is non-zero), the values in ctx->image->data.fl[]
** are still of type (single-precision) float.
**
** real xw, yw: the *world-space* position at which to do convolution
**
******** OUTPUT:
**
** There is no function return as such (the return is void); the result of
** calling this is to set other fields in the mprContext *ctx:
**
** The given code below saves in ctx->wpos[] the given (wx,wy).  You set
** ctx->ipos[] to the result of converting (wx,wy) to index-space. this is
** something you'll need to compute in any case to do convolution; saving the
** index-space position in ctx->ipos[] just provides a way to externally
** check that this first step is being done correctly (mapr wtoi should rely
** on the same coordinate transform ability).
**
** If, centered around index space position, the kernel support overlaps with
** invalid image data indices, then the convolution can't be computed.  You
** set ctx->outside to the *sum* of the number of out-of-bounds indices along
** the faster axis, plus those along the slower axis.  Otherwise, when
** convolution can be performed, set ctx->outside to 0, and set ctx->value to
** the convolution result (regardless of ctx->mode). Also, if
** mprModeNeedsGradient(ctx->mode), then save the *world-space* gradient (a
** 2-vector) of the convolution result in ctx->gradient.
**
** You should not call malloc() here in mprConvoEval.
**
** This function does not need its own error checking since it is called from
** functions that already do sufficient error checking.
*/
void
mprConvoEval(mprCtx *ctx, real xw, real yw) {
  // figure out if you need to also measure the 1st derivative
  int needgrad = mprModeNeedsGradient(ctx->mode);
  // initialize output
  ctx->wpos[0] = xw;
  ctx->wpos[1] = yw;
  /* NOTE: you need to (in code below) set ctx->ipos to
     inverse of ctx->image->ItoW, multiplied by ctx->wpos */
  ctx->outside = 0;
  ctx->value = mprNan(0);
  if (needgrad) {
    ctx->gradient[0] = ctx->gradient[1] = mprNan(0);
  }

  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real WPOS[3];
  real IPOS[3];
  WPOS[0] = xw;
  WPOS[1] = yw;
  WPOS[2] = 1;
  
  MV3_MUL(IPOS, ctx->WtoI, WPOS);
  real ix = ctx->ipos[0] = IPOS[0];
  real iy = ctx->ipos[1] = IPOS[1];
  int n1 = ctx->odd_support?floor(ix+0.5):floor(ix);
  int n2 = ctx->odd_support?floor(iy+0.5):floor(iy);
  real alpha1 = ix - n1;
  real alpha2 = iy - n2;
  // store the kernal value in two array2 before convolution
  real kvalues1[ctx->kern->support+1];
  real kvalues2[ctx->kern->support+1];
  for(int i = ctx->lower,j=0;i<=ctx->upper;i++,j++){
      kvalues1[j] = ctx->kern->eval(alpha1-i);
      kvalues2[j] = ctx->kern->eval(alpha2-i);
  }
  

  for(int i2 = ctx->lower;i2<=ctx->upper;i2++){
    if(i2+n2<0||i2+n2>(signed)(ctx->image->size[1])-1){
            ctx->outside+=1;
      }
  }
  for(int i1 = ctx->lower;i1<=ctx->upper;i1++){
    if(i1+n1<0||i1+n1>(signed)(ctx->image->size[0])-1){
            ctx->outside+=1;
      }
  }
  if(!ctx->outside){
   real cov_sum = 0;
   for(int i2 = ctx->lower,j2=0;i2<=ctx->upper;i2++,j2++)
      for(int i1 = ctx->lower,j1=0;i1<=ctx->upper;i1++,j1++){
       
        int index_1D = (n2+i2)*(ctx->image->size[0])+n1+i1;
        cov_sum += ctx->image->data.fl[index_1D]*kvalues1[j1]*kvalues2[j2];     
   }
   ctx->value = cov_sum;
  }

  if(needgrad&&!ctx->outside){
     const mprKernel* the_kern = ctx->kern->deriv;
     // store the kernal value in two array2 before convolution
    real the_kvalues1[the_kern->support+1];
    real the_kvalues2[the_kern->support+1];
    for(int i = ctx->lower,j=0;i<=ctx->upper;i++,j++){
      the_kvalues1[j] = the_kern->eval(alpha1-i);
      the_kvalues2[j] = the_kern->eval(alpha2-i);
    }

     real cov_sum_deriv1 = 0;
     real cov_sum_deriv2 = 0;
     for(int i2 = ctx->lower,j2=0;i2<=ctx->upper;i2++,j2++){
       for(int i1 = ctx->lower,j1=0;i1<=ctx->upper;i1++,j1++){
        int the_index_1D = (n2+i2)*(ctx->image->size[0])+n1+i1;
        cov_sum_deriv1 += ctx->image->data.fl[the_index_1D]*the_kvalues1[j1]*kvalues2[j2];
       }
      }

      for(int i2 = ctx->lower,j2=0;i2<=ctx->upper;i2++,j2++){
       for(int i1 = ctx->lower,j1=0;i1<=ctx->upper;i1++,j1++){
        int the_index_1D = (n2+i2)*(ctx->image->size[0])+n1+i1;
        cov_sum_deriv2 += ctx->image->data.fl[the_index_1D]*kvalues1[j1]*the_kvalues2[j2];
       }
      }
    // convert gradients from index space to world space
    real g[2];
    g[0] = cov_sum_deriv1;
    g[1] = cov_sum_deriv2;
    real gw[2];
    MV2_MUL(gw, ctx->deriv_ItoW, g);
    ctx->gradient[0] = gw[0];
    ctx->gradient[1] = gw[1];
  }
  

  // ( 83 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

  return;
}
