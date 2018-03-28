#include "rnd.h"
#include "rndPrivate.h"

// any structs or functions used to manage multi-threaded operation
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
static void procRayRow(real* odata, uint probLen, uint len, int idRow,rndCtx *ctx,rndConvo *cnv){
  rndRay ray;
  for (uint ii=0; ii<len; ii++) {
          cnv = rndConvoNew(ctx);
          rndRayStart(&ray, ii, idRow, cnv, ctx);
          while (!ray.finished) {
            rndRayStep(&ray, cnv, ctx);
          }
          for (uint pi=0; pi<probLen; pi++) {
            odata[pi] = ray.result[pi];
          }
          if (ctx->timing) {
            odata[probLen] = ray.time;
          }
          odata += probLen+ctx->timing;
  }
  
}
// the following code fragments are modified based on the hw3 pthread reference solution.
typedef struct{
  pthread_t thread;
  pthread_mutex_t *pmut;
  uint *pnext;
  uint size[2];
  real *dout;
  int muterror;
  rndCtx *ctx;
  rndConvo *cnv;
  uint probLen;
  uint tidx;
} threadArg;

static void * threadBody(void * _arg){
  threadArg *arg = (threadArg* ) _arg;

  int pthe;
  while(1){
    pthe = pthread_mutex_lock(arg->pmut);
    if(pthe){
      fprintf(stderr,"pthread_mutex_lock error: %s(%d)\n",strerror(pthe),pthe);
      arg->muterror = pthe;
      return arg;
    }
    uint jj = *(arg->pnext);
    
    *(arg->pnext) = *(arg->pnext) + 1;
    pthe = pthread_mutex_unlock(arg->pmut);
    if(pthe){
       fprintf(stderr,"pthread_mutex_unlock error:%s(%d)\n",strerror(pthe),pthe);
      arg->muterror = pthe;
      return arg;
    }
    if(jj>=arg->size[1]){
      break;
    }
    uint sz0 = arg->size[0];
  
    procRayRow(arg->dout+jj*(arg->probLen+arg->ctx->timing)*sz0, arg->probLen, sz0, jj, arg->ctx, arg->cnv);
  }
  arg->muterror = 0;
  return arg;
}

// ( 66 lines in reference implementation )
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

static int
checkTxf(const rndCtx *ctx) {
  if (rndProbeRgba == ctx->probe
      || rndProbeRgbaLit == ctx->probe) {
    if (!(ctx->txf.len
          && ctx->txf.vmin < ctx->txf.vmax
          && ctx->txf.rgba
          && 0 < ctx->txf.alphaNear1
          && ctx->txf.alphaNear1 <= 1)) {
      biffAddf(RND, "%s: don't seem to have txf set (was rndCtxSetTxf "
               "called? \"rendr go\" may need \"-lut\" option)", __func__);
      return 1;
    }
  }
  return 0;
}

static int
checkLparm(const rndCtx *ctx) {
  if (rndProbeRgbaLit == ctx->probe) {
    if (!(isfinite(ctx->lparm.ka)
          && isfinite(ctx->lparm.kd)
          && isfinite(ctx->lparm.ks)
          && isfinite(ctx->lparm.p)
          && V3_ISFINITE(ctx->lparm.dcn)
          && V3_ISFINITE(ctx->lparm.dcf))) {
      biffAddf(RND, "%s: lighting parms don't seem to be set "
             "(rndCtxSetLparm called?)", __func__);
      return 1;
    }
  }
  return 0;
}

/*
** rndRender: the main function that computes a volume rendering.
**
** If _dhi >= 0 and _dvi >= 0: instead of doing the usual image
** rendering, just render the ray through pixel (dhi,dvi),
** which may be helpful for debugging. In this case, the passed
** img should be NULL.
*/
int
rndRender(rndImage *img, rndCtx *ctx, int _dhi, int _dvi) {
  if (!ctx) {
    biffAddf(RND, "%s: got NULL rndCtx", __func__);
    return 1;
  }
  uint dbp[2]={UINT_MAX,UINT_MAX};
  if (_dhi >= 0 && _dvi >= 0) {
    ctx->debuggingPixel = 1;
    dbp[0] = (uint)_dhi;
    dbp[1] = (uint)_dvi;
    if (img) {
      biffAddf(RND, "%s: debugging ray [%u,%u] but got non-NULL rndImage",
               __func__, dbp[0], dbp[1]);
      return 1;
    }
    if (ctx->threadNum) {
      biffAddf(RND, "%s: debugging ray [%u,%u] but asked for %u threads "
               "(can only do a single ray single-threaded)",
               __func__, dbp[0], dbp[1], ctx->threadNum);
      return 1;
    }
    printf("%s: will render single ray through pixel [%u,%u]\n", __func__,
           dbp[0], dbp[1]);
  } else {
    ctx->debuggingPixel = 0;
    if (!img) {
      biffAddf(RND, "%s: got NULL output rndImage", __func__);
      return 1;
    }
  }
  if (rndCtxCheckParms(ctx->planeSep, ctx->probe,
                       ctx->blend, ctx->threadNum)
      || checkTxf(ctx)
      // the Levoy txfs are completely optional
      || checkLparm(ctx)
      || rndCameraUpdate(&(ctx->cam))
      // find world-space coords for view-space lights
      || rndCtxUpdateLight(ctx)) {
    biffAddf(RND, "%s: problems", __func__);
    return 1;
  }

  uint plen = rndProbeLen(ctx->probe);
  uint timing = !!(ctx->timing);
  if (!ctx->debuggingPixel) {
    /* allocate output image */
    if (rndImageAlloc(img, plen + timing,
                      ctx->cam.size[0],
                      ctx->cam.size[1], rndITypeReal)) {
      biffAddf(RND, "%s: trouble allocating output", __func__);
      return 1;
    }
  }
  airArray *mop = airMopNew();
  rndConvo *cnv = rndConvoNew(ctx); // for single-threaded
  airMopAdd(mop, cnv, (airMopper)rndConvoNix, airMopAlways);
  real *odata = img->data.re;
  if (!ctx->threadNum) {
    rndRay ray; /* local ray */
    if (ctx->debuggingPixel) {
      rndRayStart(&ray, dbp[0], dbp[1], cnv, ctx);
      uint stepi = 0;
      while (!ray.finished) {
        rndRayStep(&ray, cnv, ctx);
        stepi++;
        if (!ray.finished) {
          if (rndVerbose) {
            printf("%s: still going after %u step%s\n",
                   __func__, stepi, stepi > 1 ? "s" : "");
          }
        } else {
          printf("%s: FINISHED! after %u step%s\n",
                 __func__, stepi, stepi > 1 ? "s" : "");
        }
      }
      for (uint pi=0; pi<plen; pi++) {
        printf("%s: " RCS " = ray result[%u]\n", __func__,
               ray.result[pi], pi);
      }
    } else {
      printf("%s: starting without threads\n", __func__);
      char doneStr[13];
      uint npix = ctx->cam.size[0]*ctx->cam.size[1];
      // controls if and when the progress indication is printed
      uint tick = npix/200;
      if (tick) {
        printf("%s: rendering ...       ", __func__); fflush(stdout);
      }
      uint pixi=0;
      for (uint jj=0; jj<ctx->cam.size[1]; jj++) {
        for (uint ii=0; ii<ctx->cam.size[0]; ii++) {
          if (tick && !(pixi % tick)) {
            printf("%s", airDoneStr(0, pixi, npix, doneStr)); fflush(stdout);
          }
          /* the inner loop for volume rendering: start the ray, and then keep
             taking steps along the ray until the ray is finished */
          rndRayStart(&ray, ii, jj, cnv, ctx);
          while (!ray.finished) {
            rndRayStep(&ray, cnv, ctx);
          }
          // copy ray result into output image
          for (uint pi=0; pi<plen; pi++) {
            odata[pi] = ray.result[pi];
          }
          if (timing) {
            odata[plen] = ray.time;
          }
          odata += plen + timing;
          pixi++;
        }
      }
      if (tick) {
        printf("%s\n", airDoneStr(0, pixi, npix, doneStr)); fflush(stdout);
      }
    }
  } else {
    // multi-threaded operation
    // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
    // the following code fragments are modified based on the hw3 pthread reference solution.
    uint nextRow = 0;
    pthread_mutex_t mutex;
    threadArg *tharg = (threadArg *)calloc(ctx->threadNum,sizeof(threadArg));
    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr,PTHREAD_MUTEX_ERRORCHECK);
    int pthe =  pthread_mutex_init(&mutex,&attr);
    if (pthe){
      fprintf(stderr,"%s: failed to create mutex: %s(%d)\n", __func__, strerror(pthe),pthe);
      return 1;
    }
    for(uint threadIdx=0;threadIdx<ctx->threadNum;threadIdx++){
      threadArg *tha = tharg+threadIdx;
      tha->pmut = &mutex;
      tha->pnext = &nextRow;
      tha->dout = odata;
      tha->size[0] = ctx->cam.size[0];
      tha->size[1] = ctx->cam.size[1];
      tha->probLen = plen;
      tha->ctx = ctx;
      tha->cnv = cnv;
      tha->tidx = threadIdx;
      pthe = pthread_create(&(tha->thread),NULL,threadBody,(void *)tha);
      if(pthe){
        fprintf(stderr,"%s: failed to create thread: %s(%d)\n", __func__, strerror(pthe),pthe);
        return 1;
      }
    }
    for(uint threadIdx = 0; threadIdx<ctx->threadNum;threadIdx++){
      threadArg *tha = tharg + threadIdx;
      pthe = pthread_join(tha->thread,NULL);
      if(pthe){
        fprintf(stderr,"%s: failed to join thread: %s(%d)\n", __func__, strerror(pthe),pthe);
        return 1;
      }
      if(tha->muterror){
        fprintf(stderr,"%s: thread had mutex problems\n", __func__);
        return 1;
      }
    }
    free(tharg);
    pthread_mutex_destroy(&mutex);
    // ( 36 lines in reference implementation )
    // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  }

  airMopOkay(mop);
  return 0;
}
