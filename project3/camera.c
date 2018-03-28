#include "rnd.h"
#include "rndPrivate.h"

void
rndCameraInit(rndCamera *cam) {
  assert(cam);
  V3_SET_NAN(cam->fr);
  V3_SET_NAN(cam->at);
  V3_SET_NAN(cam->up);
  cam->nc = rndNan(0);
  cam->fc = rndNan(0);
  cam->FOV = rndNan(0);
  cam->size[0] = cam->size[1] = 0;
  cam->ortho = 0;

  cam->d = rndNan(0);
  cam->ar = rndNan(0);
  M4_SET_NAN(cam->WtoV);
  M4_SET_NAN(cam->VtoW);
  cam->ncv = rndNan(0);
  cam->fcv = rndNan(0);
  cam->hght = rndNan(0);
  cam->wdth = rndNan(0);
  return;
}

rndCamera *
rndCameraNew() {
  rndCamera *cam = MALLOC(1, rndCamera); assert(cam);
  rndCameraInit(cam);
  return cam;
}

rndCamera *
rndCameraNix(rndCamera *cam) {
  if (cam) {
    free(cam);
  }
  return NULL;
}

/*
** Adds hest command-line options for all the camera parameters
** Does no error checking.
*/
void
rndCameraHestOptAdd(hestOpt **Ahopt, rndCamera *cam) {
  hestOptAdd(Ahopt, "fr", "x y z", airTypeReal, 3, 3, cam->fr, NULL,
             "look-from point");
  hestOptAdd(Ahopt, "at", "x y z", airTypeReal, 3, 3, cam->at, NULL,
             "look-at point");
  hestOptAdd(Ahopt, "up", "x y z", airTypeReal, 3, 3, cam->up, NULL,
             "up direction");
  hestOptAdd(Ahopt, "nc", "dist", airTypeReal, 1, 1, &(cam->nc), NULL,
             "at-relative near clipping distance");
  hestOptAdd(Ahopt, "fc", "dist", airTypeReal, 1, 1, &(cam->fc), NULL,
             "at-relative far clipping distance");
  hestOptAdd(Ahopt, "fov", "angle", airTypeReal, 1, 1, &(cam->FOV), NULL,
             "vertical field-of-view, in degrees. Full vertical "
             "extent of image plane subtends this angle.");
  hestOptAdd(Ahopt, "sz", "s0 s1", airTypeUInt, 2, 2, &(cam->size), NULL,
             "# samples (horz vert) of image plane. ");
  hestOptAdd(Ahopt, "ortho", NULL, airTypeInt, 0, 0, &(cam->ortho), NULL,
             "use orthographic instead of (the default) "
             "perspective projection ");
  return;
}

/* some really ugly macros for repetitive checking of whether
   fields of rndCamera have been set */
#define CHECK(F)                                        \
  if (!isfinite(cam->F)) {                              \
    biffAddf(RND, "%s: cam->" #F " not set", __func__); \
    return 1;                                           \
  }
#define CHECKVEC(V)                                     \
  if (!V3_ISFINITE(cam->V)) {                           \
    biffAddf(RND, "%s: cam->" #V " not set", __func__); \
    return 1;                                           \
  }
#define CHECKMAT(M)                                     \
  if (!M3_ISFINITE(cam->M)) {                           \
    biffAddf(RND, "%s: cam->" #M " not set", __func__); \
    return 1;                                           \
  }

/*
** rndCameraUpdate: implement the computations of FSV Section 5.3 to set the
** "output" camera parameters based on the "input" parameters.  The input
** parameters may have been set via rndCameraSet(), or by directly setting the
** struct members (as with rndCameraHestOptAdd() and then hestParseOrDie())
**
** Repeatedly calling this should have no adverse effect (to facilitate its
** use for cautious error checking)
*/
int
rndCameraUpdate(rndCamera *cam) {
  if (!cam) {
    biffAddf(RND, "%s: got NULL pointer", __func__);
    return 1;
  }
  CHECKVEC(fr);
  CHECKVEC(at);
  CHECKVEC(up);
  CHECK(nc);
  CHECK(fc);
  CHECK(FOV);
  if (!( 0 < cam->FOV && cam->FOV < 160 )) {
    biffAddf(RND, "%s: need 0 < FOV > 160 degrees (not " RCS ")",
             __func__, cam->FOV);
    return 1;
  }
  if (!( cam->size[0] && cam->size[1] )) {
    biffAddf(RND, "%s: need two non-zero cam->size[] (not %u %u)",
             __func__, cam->size[0], cam->size[1]);
    return 1;
  }
  cam->ar = ((real)cam->size[0])/cam->size[1];
  /* To-do here: set cam fields d, u, v, n, WtoV, VtoW, ncv, fcv, hght, wdth
     (yes, all of these!) based on previous fields (see rnd.h). No further
     error handling is required; the user may get NaNs or garbage in return
     for bad camera parameters (e.g. fr == at, or up aligned with
     fr-at). Reference code does do error checking for these conditions */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  real atof[3];
  V3_SUB(atof, cam->fr, cam->at);
  cam->d = V3_LEN(atof);
  real tmpp;
  V3_NORM(cam->n, atof, tmpp);
  real u[3];
  V3_CROSS(u, cam->up, cam->n);
  V3_NORM(cam->u, u, tmpp);
  V3_CROSS(cam->v, cam->n, cam->u);
  M4_SET(cam->WtoV,cam->u[0],cam->u[1],cam->u[2],-V3_DOT(cam->fr,cam->u),
                   cam->v[0],cam->v[1],cam->v[2],-V3_DOT(cam->fr,cam->v),
                   cam->n[0],cam->n[1],cam->n[2],-V3_DOT(cam->fr,cam->n),
                   0,0,0,1);
  rnd_m4_affine_inv(cam->VtoW, cam->WtoV);
  
  cam->ncv = cam->nc+cam->d;
  cam->fcv = cam->fc+cam->d;
  cam->hght = 2*tan(cam->FOV/2*M_PI/180)*cam->d;
  cam->wdth = cam->hght*((real)cam->size[0])/cam->size[1];
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return 0;
}

/*
** rndCameraSet: sets up an rndCamera
*/
int
rndCameraSet(rndCamera *cam,
             const real fr[3],
             const real at[3],
             const real up[3],
             real nc, real fc,
             real FOV,
             uint size0, uint size1,
             int ortho) {
  if (!(cam && fr && at && up)) {
    biffAddf(RND, "%s: got NULL pointer", __func__);
    return 1;
  }
  V3_COPY(cam->fr, fr);
  V3_COPY(cam->at, at);
  V3_COPY(cam->up, up);
  cam->nc = nc;
  cam->fc = fc;
  cam->FOV = FOV;
  cam->size[0] = size0;
  cam->size[1] = size1;
  cam->ortho = !!ortho;
  if (rndCameraUpdate(cam)) {
    biffAddf(RND, "%s: trouble with camera parms", __func__);
    return 1;
  }
  return 0;
}

/*
** rndCameraCheck: does minimal check to tell whether rndCameraUpdate() seem
** to have been called
*/
int
rndCameraCheck(const rndCamera *cam) {
  if (!cam) {
    biffAddf(RND, "%s: got NULL pointer", __func__);
    return 1;
  }
  CHECK(ar);
  CHECK(d);
  CHECKVEC(u);
  CHECKVEC(v);
  CHECKVEC(n);
  CHECKMAT(WtoV);
  CHECKMAT(VtoW);
  CHECK(ncv);
  CHECK(fcv);
  CHECK(hght);
  CHECK(wdth);
  return 0;
}

#undef CHECK
#undef CHECKVEC
#undef CHECKMAT
