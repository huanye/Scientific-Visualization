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

/* You may (if you want) define new #define macros and functions here; but
   they should only be used in this file. Thus, prefix any function
   declarations here with "static" */
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
static void addSegment(mprArray *segArr,real x1,real y1,real x2,real y2){
   uint segIdx = mprArrayLenIncr(segArr, 1);
   real *ss = (real*)(segArr->data) + 4*segIdx;
   ss[0] = x1;
   ss[1] = y1;
   ss[2] = x2;
   ss[3] = y2;
}

static real lerp(real omin, real omax, real imin, real xx, real imax) {
  real ret=0;
  real alpha = (xx-imin)/(imax-imin);
  ret = (1-alpha)*omin+alpha*omax;
  return ret;
}
// ( 12 lines in reference implementation )
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

/*
******** mprIsocontour: run Marching Squares on given mprImage to compute
** isocontour at isovalue "isoval", by setting segArr to contain the line
** segments composing the isocontour.
**
******** INPUT:
**
** mprCtx *ctx: context containing the image (ctx->image) to isocontour
**
** isoval: the isovalue at which to compute the isocontour
**
** skipOutOfView: 23710 students can ignore and the remaining arguments.
** 33710 students: if this is non-zero, then the size0, size1, centX, centY,
** fov, and angle arguments define a rectangular viewing window, which mean
** the same as for mprPictureSample.  size0 and size1 don't deterine a raster
** sampling (as for mprPictureSample), only the aspect ratio of the window.
** When traversing cells for Marching Squares, skip over the cells that fall
** outside the viewing window, as determined by whether the *center* of the
** cell (the center of the parallelogon bounded by four vertices) is inside
** the view rectangle.  Otherwise, if skipOutOfView is zero, then all cells
** may contribute to the isocontour.  HINT: use mprPictureItoW: it takes the
** same arguments to create an ItoW matrix, which you can invert. However
** you solve this, add comments (for style points) to document how it works.
**
******** OUTPUT:
**
** At completion, segArr->data should (cast as a real*) point to what is
** effectively a 4 (fast) by N (slow) array of world-space vertex coordinates,
** giving the end-points of the N line segments making up the isocontour.  The
** 4 values along the faster axis are x0 y0 x1 y1, to represent the segment
** from (x0,y0) to (x1,y1). For 33710 students: if skipOutOfView is non-zero,
** then this list should not include isocontour segments falling outside the
** viewing rectangle.
**
** The topology of the isocontour should match that produced by true bilinear
** interpolation, though naturally the isocontour elements are straight line
** segments rather than hyperbolas.
**
** If there are no isocontour segments at isoval, then leave segArr unchanged.
*/
int
mprIsocontour(mprArray *segArr, const mprImage *img,
              real isoval, int skipOutOfView,
              uint size0, uint size1,
              real centX, real centY, real fov,
              real angle) {
  if (!(segArr && img)) {
    biffAddf(MPR, "%s: got NULL pointer (%p,%p)", __func__,
             (void*)segArr, (void*)img);
    return 1;
  }
  if (1 != img->channel) {
    biffAddf(MPR, "%s: only works on scalar (not %u-channel) images",
             __func__, img->channel);
    return 1;
  }
  if (mprTypeFloat != img->dtype) {
    biffAddf(MPR, "%s: only works on %s (not %s) pixels", __func__,
             airEnumStr(mprType_ae, mprTypeFloat),
             airEnumStr(mprType_ae, img->dtype));
    return 1;
  }
  if (4*sizeof(real) != segArr->unit) {
    biffAddf(MPR, "%s: need mprArray->unit == 4*sizeof(real) == %u, "
             "not %u", __func__, (uint)(4*sizeof(real)),
             (uint)(segArr->unit));
    return 1;
  }
  if (!isfinite(isoval)) {
    biffAddf(MPR, "%s: got non-finite isoval %g", __func__, isoval);
    return 1;
  }
  /* whether or not any previous segments had been stored in segArr->data, now
     we're starting a new isocontour, so we reset the nominal array length.
     We don't support accumulating segments from multiple calls (e.g. from
     multiple distinct isocontours) into one mprArray */
  segArr->len = 0;
  /* Now, to add a new segment to the isocontour, you:
       uint segIdx = mprArrayLenIncr(segArr, 1);
       real *ss = (real*)(segArr->data) + 4*segIdx;
       ss[0] = ... 1st endpoint X coord ...;
       ss[1] = ... 1st endpoint Y coord ...;
       ss[2] = ... 2nd endpoint X coord ...;
       ss[3] = ... 2nd endpoint Y coord ...;
     where all endpoint coordinates are in *world-space*.  Because of
     re-allocation, the (pointer) value of segArr->data may change after
     mprArrayLenIncr(). It may be convenient to write a new (static)
     function above for adding one segment. */

  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  uint fast = img->size[0];
  uint slow = img->size[1];
  for(int i2=0;i2<(signed)slow-1;i2++)
    for(int i1=0;i1<(signed)fast-1;i1++){

      if(skipOutOfView){
         // create a ItoW matrix using mprPictureItoW
        real localItoW[9];
        mprPictureItoW(localItoW,size0,size1,centX,centY,fov,1,1,angle,0,0);
        // inverse the ItoW matrix to get the WtoI matrix;
        real localWtoI[9];
        real tmp;
        M3_INVERSE(localWtoI, localItoW, tmp);
        //convert the current point indices to the location in the world space
        real ipos[3],wpos[3];
        V3_SET(ipos,i1,i2,1);
        MV3_MUL(wpos, img->ItoW, ipos);
        //convert the current location in the world space to the local indices
        real localipos[3];
        MV3_MUL(localipos, localWtoI, wpos);
        // check if the localipos is within the range defined by size0 and size1
        if(localipos[0]<0||localipos[0]>size0-1
           ||localipos[1]<0||localipos[1]>size1-1)
          continue;
      }

          int index1 = i2*fast+i1;
          int index2 = i2*fast+i1+1;
          int index3 = (i2+1)*fast+i1;
          int index4 = (i2+1)*fast+i1+1;

          real a = (real)(img->data.fl[index1]);
          real b = (real)(img->data.fl[index2]);
          real c = (real)(img->data.fl[index3]);
          real d = (real)(img->data.fl[index4]);

          // if one of a, b, c and d exactly equals to isoval, 
          // subtract that corner value by epsilon = 0.00001
          real epsilon = 0.00001;
          if(a==isoval)
            a-=epsilon;
          else if(b==isoval)
            b-=epsilon;
          else if(c==isoval)
            c-=epsilon;
          else if(d==isoval)
            d-=epsilon;

          int idx = 8*(a>isoval)+4*(b>isoval)+2*(c>isoval)+(d>isoval);
          if(idx==0||idx==15)
            continue;
          // precalculate all "zeros" no matter which case
          real x_ab = lerp(i1, i1+1, a, isoval, b);
          real y_ab = lerp(i2, i2, a, isoval, b);
          real i_ab[3];
          V3_SET(i_ab, x_ab, y_ab, 1);
          real x_bd = lerp(i1+1, i1+1, b, isoval, d);
          real y_bd = lerp(i2, i2+1, b, isoval, d);
          real i_bd[3];
          V3_SET(i_bd, x_bd, y_bd, 1);
          real x_dc = lerp(i1+1, i1, d, isoval, c);
          real y_dc = lerp(i2+1, i2+1, d, isoval, c);
          real i_dc[3];
          V3_SET(i_dc, x_dc, y_dc, 1);
          real x_ca = lerp(i1, i1, c, isoval, a);
          real y_ca = lerp(i2+1, i2, c, isoval, a);
          real i_ca[3];
          V3_SET(i_ca, x_ca, y_ca, 1);
          // the world space coordinates of all 4 zero crossings
          real w_ab[3];
          real w_bd[3];
          real w_dc[3];
          real w_ca[3];
          // convert from index space to world space
          MV3_MUL(w_ab, img->ItoW, i_ab);
          MV3_MUL(w_bd, img->ItoW, i_bd);
          MV3_MUL(w_dc, img->ItoW, i_dc);
          MV3_MUL(w_ca, img->ItoW, i_ca);
   
          if(idx==9||idx==6){// four ZCs
               // estimate the central data value by average the four corner values
               real central = (a+b+c+d)/4.0;
               if(central>isoval){
                   if(idx==9){
                       addSegment(segArr,w_ab[0],w_ab[1],w_bd[0],w_bd[1]);
                       addSegment(segArr,w_dc[0],w_dc[1],w_ca[0],w_ca[1]);
                   }
                   else{
                        addSegment(segArr,w_ab[0],w_ab[1],w_ca[0],w_ca[1]);
                        addSegment(segArr,w_bd[0],w_bd[1],w_dc[0],w_dc[1]);
                   }
               }
          }
          // the rest cases only have 2 zero crossings
          else if(idx==8||idx==7){ 
               addSegment(segArr,w_ab[0],w_ab[1],w_ca[0],w_ca[1]);
          }
          else if(idx==4||idx==11){
              addSegment(segArr,w_ab[0],w_ab[1],w_bd[0],w_bd[1]);
          }
          else if(idx==1||idx==14){
              addSegment(segArr,w_bd[0],w_bd[1],w_dc[0],w_dc[1]);
          }
          else if(idx==2||idx==13){
              addSegment(segArr,w_ca[0],w_ca[1],w_dc[0],w_dc[1]);
          }
          else if(idx==12||idx==3){
              addSegment(segArr,w_ca[0],w_ca[1],w_bd[0],w_bd[1]);
          }
          else if(idx==5||idx==10){
              addSegment(segArr,w_ab[0],w_ab[1],w_dc[0],w_dc[1]);
          }
    }
  // ( 129 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

  return 0;
}
