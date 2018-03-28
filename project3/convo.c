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

/*
** rndConvoNew: create an rndConvo based on state set up in
** rndCtx; this will be called after the rndCtx has been fully
** set up, and just prior to rndRayStart.  This is the place to
** allocate buffers (e.g. arrays as long as the kernel support)
** that will be needed during convolution, but need to be allocated
** per thread; the rndCtx is shared between threads.
**
** Anything allocated here should be freed by rndConvoNix
**
** No biff-based error reporting is necessary: you should use assert()
** to test the success of these allocations
*/
rndConvo *
rndConvoNew(const rndCtx *ctx) {
  rndConvo *cnv = MALLOC(1, rndConvo); assert(cnv);
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  cnv->kvalues1 = MALLOC(ctx->upper-ctx->lower+2,real);assert(cnv->kvalues1);
  cnv->kvalues2 = MALLOC(ctx->upper-ctx->lower+2,real);assert(cnv->kvalues2);
  cnv->kvalues3 = MALLOC(ctx->upper-ctx->lower+2,real);assert(cnv->kvalues3);
  cnv->the_kvalues1 = MALLOC(ctx->upper-ctx->lower+2,real);assert(cnv->the_kvalues1);
  cnv->the_kvalues2 = MALLOC(ctx->upper-ctx->lower+2,real);assert(cnv->the_kvalues2);
  cnv->the_kvalues3 = MALLOC(ctx->upper-ctx->lower+2,real);assert(cnv->the_kvalues3);
  cnv->grad_need = rndCtxNeedGradient(ctx);

  // ( 4 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return cnv;
}

rndConvo *
rndConvoNix(rndConvo *cnv) {
  if (cnv) {
    // Free anything dynamically allocated in this rndConvo
    // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
    free(cnv->kvalues1);
    free(cnv->kvalues2);
    free(cnv->kvalues3);
    free(cnv->the_kvalues1);
    free(cnv->the_kvalues2);
    free(cnv->the_kvalues3);
    // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
    free(cnv);
  }
  return NULL;
}

/*
******** rndConvoEval: 3D convolution of float or short data
**
** INPUT:
** world-space positions (xw,yw,zw) at which to evaluate convolution
** rndConvo *cnv and rndCtx *ctx structs
**
** OUTPUT:
** Set cnv->inside to 1 if convolution could be computed (no need
** for data values outside valid data indices).
** Set cnv->value to value reconstructed by convolution
** If rndCtxNeedGradient(ctx) (though avoid calling that here),
** set cnv->gradient to world-space gradient of value
*/
void
rndConvoEval(real xw, real yw, real zw,
             rndConvo *cnv, const rndCtx *ctx) {
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real WPOS[4];
  real IPOS[4];
  WPOS[0] = xw;
  WPOS[1] = yw;
  WPOS[2] = zw;
  WPOS[3] = 1;
  
  MV4_MUL(IPOS, ctx->WtoI, WPOS);

  


  int n1 = ctx->odd_support?floor(IPOS[0]+0.5):floor(IPOS[0]);
  int n2 = ctx->odd_support?floor(IPOS[1]+0.5):floor(IPOS[1]);
  int n3 = ctx->odd_support?floor(IPOS[2]+0.5):floor(IPOS[2]);

  real alpha1 = IPOS[0] - n1;
  real alpha2 = IPOS[1] - n2;
  real alpha3 = IPOS[2] - n3;

  

  cnv->inside = 1;

  for(int i = ctx->lower;i<=ctx->upper;i++){
    if(i+n3<0||i+n3>(signed)(ctx->vol->size[2])-1
      ||i+n2<0||i+n2>(signed)(ctx->vol->size[1])-1
      ||i+n1<0||i+n1>(signed)(ctx->vol->size[0])-1){
            cnv->inside = 0;
      }
  }
 
  real cov_sum = 0;  
  real cov_sum_deriv1 = 0;
  real cov_sum_deriv3 = 0;
  real cov_sum_deriv2 = 0;

  if(cnv->inside){
   
   for(int i = ctx->lower,j=0;i<=ctx->upper;i++,j++){
      cnv->kvalues1[j] = ctx->kern->eval(alpha1-i);
      cnv->kvalues2[j] = ctx->kern->eval(alpha2-i);
      cnv->kvalues3[j] = ctx->kern->eval(alpha3-i);
  }
   if (cnv->grad_need){
   
   for(int i = ctx->lower,j=0;i<=ctx->upper;i++,j++){
      cnv->the_kvalues1[j] = ctx->the_kern->eval(alpha1-i);
      cnv->the_kvalues2[j] = ctx->the_kern->eval(alpha2-i);
      cnv->the_kvalues3[j] = ctx->the_kern->eval(alpha3-i);
     }
  }

   for(int i3 = ctx->lower,j3=0;i3<=ctx->upper;i3++,j3++)
    for(int i2 = ctx->lower,j2=0;i2<=ctx->upper;i2++,j2++)
      for(int i1 = ctx->lower,j1=0;i1<=ctx->upper;i1++,j1++){
       
        int index_1D = (n3+i3)*(ctx->vol->size[1]*ctx->vol->size[0])+
                        (n2+i2)*(ctx->vol->size[0])+n1+i1;
        if(ctx->vol->dtype == rndTypeFloat)
          cov_sum += ctx->vol->data.fl[index_1D]*cnv->kvalues1[j1]*cnv->kvalues2[j2]*cnv->kvalues3[j3]; 
        else if(ctx->vol->dtype == rndTypeShort)
          cov_sum += ctx->vol->data.ss[index_1D]*cnv->kvalues1[j1]*cnv->kvalues2[j2]*cnv->kvalues3[j3]; 
         if (cnv->grad_need){
          if(ctx->vol->dtype == rndTypeFloat){
          cov_sum_deriv1 += ctx->vol->data.fl[index_1D]*cnv->the_kvalues1[j1]*cnv->kvalues2[j2]*cnv->kvalues3[j3];
          cov_sum_deriv2 += ctx->vol->data.fl[index_1D]*cnv->kvalues1[j1]*cnv->the_kvalues2[j2]*cnv->kvalues3[j3];
          cov_sum_deriv3 += ctx->vol->data.fl[index_1D]*cnv->kvalues1[j1]*cnv->kvalues2[j2]*cnv->the_kvalues3[j3];  
          } else if(ctx->vol->dtype == rndTypeShort){
          cov_sum_deriv1 += ctx->vol->data.ss[index_1D]*cnv->the_kvalues1[j1]*cnv->kvalues2[j2]*cnv->kvalues3[j3];
          cov_sum_deriv2 += ctx->vol->data.ss[index_1D]*cnv->kvalues1[j1]*cnv->the_kvalues2[j2]*cnv->kvalues3[j3];
          cov_sum_deriv3 += ctx->vol->data.ss[index_1D]*cnv->kvalues1[j1]*cnv->kvalues2[j2]*cnv->the_kvalues3[j3]; 
        }
    
    }

  }
  if (cnv->grad_need){
    real g[3];
    g[0] = cov_sum_deriv1;
    g[1] = cov_sum_deriv2;
    g[2] = cov_sum_deriv3;

    

    real gw[3];
    
    MV3_MUL(gw, ctx->deriv_ItoW, g);
    //V3_NORM(gw, gw, tmp);
    cnv->gradient[0] = gw[0];
    cnv->gradient[1] = gw[1];
    cnv->gradient[2] = gw[2];
  }
    cnv->value = cov_sum;
}
  

  
  
  // ( 106 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return;
}
