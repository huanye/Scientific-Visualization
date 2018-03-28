/*
** plt and plotr: SciVis-2018 Project 1
** Copyright (C)  2018 University of Chicago. All rights reserved.
**
** This file is distributed to students in the 2018 CMSC 23710/33710
** ("SciVis") class, and no one else, for use in the SciVis class. This
** file is not open-source licensed, nor licensed for any other distribution.
** You should not allow this file to be copied or downloaded by someone who
** is not a 2018 SciVis student, as that would violate the copyright
** held by University of Chicago.
*/

#include "plt.h"

/*
******** pltLerp(omin,omax,imin,xx,imax): return value yy that ranges
** from omin to omax as xx ranges from imin to imax. Or:
**
**     yy - omin        xx - imin
**    -----------   =  -----------
**    omax - omin      imax - imin
**
** It is the caller's responsibility to make sure imax != imin.
** No error checking (for imax == imax, or any non-finite values)
** is required here.
*/
float
pltLerp(float omin, float omax, float imin, float xx, float imax) {
  float ret=0;
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  float alpha = (xx-imin)/(imax-imin);
  ret = (1-alpha)*omin+alpha*omax;
  // ( 2 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return ret;
}

/*
******** pltItoW: convert from index to world space, for given data
******** pltWtoI: convert from world to index space, for given data
**
** No error is checking needed for either pltItoW or pltWtoI.
** You can assume that:
** - data is non-NULL,
** - data->center is either pltCenterCell or pltCenterNode,
** - data->len is large enough for given centering, and
** - data->min != data->max
**
** Note that you will not actually be looking at data values in data->vv.
** Your code will use data->center, data->len, data->min, and data->max.
**
** You will lose style points if your pltItoW and pltWtoI code does not
** rely on pltLerp. You really just need two ways of calling pltLerp,
** depending on the value of data->center.
*/
float
pltItoW(const pltData *data, float indexPos) {
  float ret=0;
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  if(data->center == pltCenterNode)
    ret = pltLerp(data->min,data->max,0,indexPos,(data->len)-1);
  else
    ret = pltLerp(data->min,data->max,-0.5,indexPos,(data->len)-0.5);
  // ( 5 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return ret;
}

float
pltWtoI(const pltData *data, float worldPos) {
  float ret=0;
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  if(data->center == pltCenterNode)
    ret = pltLerp(0,(data->len)-1,data->min,worldPos,data->max);
  else
    ret = pltLerp(-0.5,(data->len)-0.5,data->min,worldPos,data->max);
  // ( 5 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return ret;
}
