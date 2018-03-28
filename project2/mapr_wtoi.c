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

#define INFO "Convert from world-space to index-space"
static char *infoLong =
  (INFO " coordinates. This relies on whatever you do with the "
   "mprImage->ItoW matrix in order to create an mprCtx in mprCtxNew; "
   "you'll have to be doing something with mprImage->ItoW so that "
   "mprConvoEval can take word-space positions and compute a "
   "convolution in index-space.");

static int
wtoiMain(int argc, const char **argv, const char *me,
         hestParm *hparm) {
  airArray *mop = airMopNew();
  hestOpt *hopt = NULL;
  char *fname;
  hestOptAdd(&hopt, NULL, "data", airTypeString, 1, 1, &fname, NULL,
             "filename of image data");
  real wpos2[2];
  hestOptAdd(&hopt, NULL, "w0 w1", airTypeReal, 2, 2, wpos2, NULL,
             "2-vector position in world-space, "
             "to convert to index-space of data");
  hestParseOrDie(hopt, argc, argv, hparm,
                 me, infoLong, AIR_TRUE, AIR_TRUE, AIR_TRUE);
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  // . . . . . . . . . . . . . . . . . useful code here
  mprImage *img = mprImageNew();
  airMopAdd(mop, img, (airMopper)mprImageNix, airMopAlways);
  if (mprImageLoad(img, fname)) {
    ERROR(MPR, "trouble loading image data");
  }
  real wpos[3];
  wpos[0] = wpos2[0];
  wpos[1] = wpos2[1];
  wpos[2] = 1;
  mprCtx *ctx = mprCtxNew(img, mprKernelTent,
                          mprModeValue, NULL, NULL);
  if (!ctx) {
    ERROR(MPR, "trouble setting up context");
  }

  /* As part of creating the new mprCtx *ctx, you should have done
     something with img->ItoW so that later calls to mprConvoEval can
     take a *world*-space position and compute the convolution in
     *index*-space.  Either by copying the information you computed
     and stored in the mprCtx, or by computing something here,
     store in the WtoI variable here the matrix that inverts img->ItoW.
     Written in the convention layout, the elements would be:

          WtoI[0]   WtoI[1]   WtoI[2]
          WtoI[3]   WtoI[4]   WtoI[5]
          WtoI[6]   WtoI[7]   WtoI[8]

     You do need to set information in WtoI, because it is used
     below for printing out the matrix-vector multiplication of
     ipos = WtoI*wpos (and also, because WtoI is initialized
     here to a bad matrix). */
  real WtoI[9];
  M3_SET_NAN(WtoI);
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real tmp;
  M3_INVERSE(WtoI, ctx->image->ItoW, tmp);
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

  real ipos[3];
  MV3_MUL(ipos, WtoI, wpos);
  // have computed ipos, now print it
  uint nlen = mprSprintLenMax(0, WtoI, 9);
  nlen = mprSprintLenMax(nlen, ipos, 3);
  nlen = mprSprintLenMax(nlen, wpos, 3);
  char buff[512];
  printf("ipos = WtoI * wpos:\n");
  sprintf(buff, "%% %u." RCD " g    %%c    %% %u." RCD " g  %% %u." RCD
          "g  %% %u." RCD "g    %%c    %% %u." RCD " g  \n",
          nlen, nlen, nlen, nlen, nlen);
  printf(buff, ipos[0], ' ', WtoI[0], WtoI[1], WtoI[2], ' ', wpos[0]);
  printf(buff, ipos[1], '=', WtoI[3], WtoI[4], WtoI[5], '*', wpos[1]);
  printf(buff, ipos[2], ' ', WtoI[6], WtoI[7], WtoI[8], ' ', wpos[2]);
  // ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' useful code here

  airMopOkay(mop);
  return 0;
}

unrrduCmd mpr_wtoiCmd = { "wtoi", INFO, wtoiMain, AIR_FALSE };
