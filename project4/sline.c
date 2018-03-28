#include "vcr.h"
#include "vcrPrivate.h"

static void
slineInit(vcrSline *sln) {
  real nan = vcrNan(0);
  sln->halfLen = sln->forwNum = sln->backNum = 0;
  V2_SET(sln->vecSeed, nan, nan);
  M2_SET(sln->jacSeed, nan, nan, nan, nan);
  return;
}

vcrSline *
vcrSlineNew(void) {
  vcrSline *sln = MALLOC(1, vcrSline); assert(sln);
  sln->pos = NULL;
  slineInit(sln);
  return sln;
}

/*
** vcrSlineAlloc: allocates the coordinate array according
** to halfLen, and resets the rest of the state
*/
int
vcrSlineAlloc(vcrSline *sln, uint halfLen) {
  if (!sln) {
    biffAddf(VCR, "%s: got NULL pointer", __func__);
    return 1;
  }
  int doalloc;
  if (!(sln->pos)) {
    doalloc = 1;
  } else {
    // already allocated for something */
    if (sln->halfLen != halfLen) {
      // but its not the right length
      free(sln->pos);
      sln->pos = NULL;
      doalloc = 1;
    } else {
      // already have right allocation
      doalloc = 0;
    }
  }
  if (doalloc) {
    sln->pos = MALLOC(2*(1 + 2*halfLen), real);
    if (!(sln->pos)) {
      biffAddf(VCR, "%s: couldn't allocate pos for %u points "
               "(halfLen %u)", __func__, 2*(1 + 2*halfLen), halfLen);
      return 1;
    }
  }
  slineInit(sln);
  sln->halfLen = halfLen;
  return 0;
}

vcrSline *
vcrSlineNix(vcrSline *sln) {
  if (sln) {
    free(sln->pos);
    free(sln);
  }
  return NULL;
}

/* You may define here new functions or #define to be used in this file.
   Function declarations here must begin with "static".  It might be useful,
   for example, to define a function that computes a single step in the
   integration, that is, the change-in-position vector dpos, at current i-th
   streamline vertex coordinate p, so that the (i+1)-th vertex is at p +
   dpos. This will need to take p, the integration parameters (stepsize and
   integration method), the streamline direction (up or downstream), and the
   vcrCtx. Keep in mind that for RK2 and RK4, multiple calls to vcrConvoEval
   will be required to compute a single step, and on any one of them you could
   land outside the field (so ctx->inside is zero); that situation has to be
   indicated to whatever is calling this singlestep function. FYI, the
   reference implementation does this checking within a #define macro that
   includes the call to vcrConvoEval.  With that macro in in place,
   finding the Euler step is a single line of code, RK2 is 4 lines, and RK4 is
   12 lines.  Note that for computing streamlines, the final ("getjac")
   argument to vcrConvoEval should always be 0 (false).
*/
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
static void 
vcrStep(real K[2], real X[2], real hh, int norm, int sign, vcrCtx *ctx) {
    real V[2];
    vcrConvoEval(ctx, X[0], X[1], 0);
    V2_COPY(V, ctx->vec);
    if (norm) {
        real temp;
        V2_NORM(V, V, temp);}
    V2_SCALE(K, sign*hh, V);
} // stores step into K[2] given X[2]

static void
vcrCalcDpos(real *Xout, real *X, real hh, int intg, int norm, int sign, vcrCtx *ctx) {   
    if (intg == 1) {                    // do Euler
        real K1[2];
        vcrStep(K1, X, hh, norm, sign, ctx);
        V2_ADD(Xout, X, K1);
    }
    if (intg == 2) {                    // do RK2
        real K1[2], K2[2], X1[2], K12[2];   // K12 being K1/2
        vcrStep(K1, X, hh, norm, sign, ctx);
        V2_SCALE(K12, 1.0/2, K1);
        V2_ADD(X1, X, K12);
        vcrStep(K2, X1, hh, norm, sign, ctx);
        V2_ADD(Xout, X, K2);
    }
    if (intg == 3) {                    // do RK4
        real K1[2], K2[2], K3[2], K4[2];
        real K12[2], K22[2], K16[2], K23[2], K33[2], K46[2];
        real X1[2], X2[2], X3[2];
        
        vcrStep(K1, X, hh, norm, sign, ctx);
        V2_SCALE(K12, 1.0/2, K1);
        V2_ADD(X1, X, K12);
        
        vcrStep(K2, X1, hh, norm, sign, ctx);
        V2_SCALE(K22, 1.0/2, K2);
        V2_ADD(X2, X, K22);
        
        vcrStep(K3, X2, hh, norm, sign, ctx);
        V2_ADD(X3, X, K3);
        
        vcrStep(K4, X3, hh, norm, sign, ctx);
        
        V2_SCALE(K16, 1.0/6, K1);
        V2_SCALE(K23, 1.0/3, K2);
        V2_SCALE(K33, 1.0/3, K3);
        V2_SCALE(K46, 1.0/6, K4);
        
        real tmp[2];
        V2_ADD(tmp, K16, K23);
        V2_ADD(tmp, tmp, K33);
        V2_ADD(tmp, tmp, K46);
        V2_ADD(Xout, X, tmp);
    }
}
   
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (43L in ref)

/*
** vcrSlineTrace: compute one streamline starting at ("seedX", "seedY"), for
** at most "halfLen" steps forward (downstream, following vector field) and
** backward (upstream, following negative vector field), with integration
** stepsize "hh", according to other parameters:
**
** if "normalize" is non-zero: integrate not the vector V(v) reconstructed by
** convolution but rather V(x)/|V(x)|. This is useful for LIC.
**
** integrate according to integration method "intg": vcrIntgEuler,
** vcrIntgMidpoint, or vcrIntgRK4.
**
** The given vcrCtx* "ctx" should be passed to vcrConvoEval.
*/
int
vcrSlineTrace(vcrSline *sln, real seedX, real seedY,
              uint halfLen, real hh,
              int normalize, int intg, vcrCtx *ctx) {
  USED(normalize);
  if (!(sln && ctx)) {
    biffAddf(VCR, "%s: got NULL pointer (%p,%p)", __func__,
             (void *)sln, (void *)ctx);
    return 1;
  }
  real seed[2];
  V2_SET(seed, seedX, seedY);
  if (!V2_ISFINITE(seed)) {
    biffAddf(VCR, "%s: seed location (%g,%g) not finite", __func__,
             seed[0], seed[1]);
    return 1;
  }
  if (!(hh > 0)) {
    biffAddf(VCR, "%s: given step size %g not > 0", __func__, hh);
    return 1;
  }
  if (airEnumValCheck(vcrIntg_ae, intg)) {
    biffAddf(VCR, "%s: integration %d not a valid %s", __func__, intg,
             vcrIntg_ae->name);
    return 1;
  }
  /* note that this does not reallocate with every call,
     it only reallocates when needed */
  if (vcrSlineAlloc(sln, halfLen)) {
    biffAddf(VCR, "%s: couldn't allocate streamline", __func__);
    return 1;
  }

  // initialize output
  V2_COPY(sln->pos + 2*halfLen, seed);
  // try reconstructing at seed point
  vcrConvoEval(ctx, seed[0], seed[1],
               1 /* getjac, since it might be needed by the prop
                    passed to vcrLIC, which might be calling us */);
  if (!(ctx->inside)) {
    // seed point outside field, nowhere to go
    sln->seedInside = 0;
    real nan = vcrNan(0);
    V2_SET(sln->vecSeed, nan, nan);
    M2_SET(sln->jacSeed, nan, nan, nan, nan);
    return 0;
  } else {
    sln->seedInside = 1;
    V2_COPY(sln->vecSeed, ctx->vec);
    M2_COPY(sln->jacSeed, ctx->jac);
  }

  /*
  ** Todo: Integrate the streamline. That is:
  ** for both backwards/upstream and forwards/downstream directions:
  ** -- (If at any time in the following the vector field needs to be
  **    reconstructed at a location that was "outside" (ctx->inside was set
  **    by 0 by vcrConvoEval), then that terminates the integration along
  **    that direction. Furthermore, the saved streamline vertices should all
  **    be locations for which ctx->inside would be non-zero, were
  **    vcrConvoEval to be called at that location.
  ** -- Using integration method intg with stepsize hh, find the vector, call
  **    it dpos, for stepping away from seed position in current direction
  **    (either up or downstream). If normalize, integrate V(x)/|V(x)| rather
  **    than V(x)
  ** -- For at most halfLen steps:
  **    -- update current position by dpos
  **    -- save new position to the correct index in the sln->pos array; see
  **       definition of vcrSline vcr.h for the indexing details basically in
  **       the forward direction, the index into sln->pos will increase, but
  **       will decrease for the backward direction
  **    -- if not on last step: compute next step
  ** Record the number of time the inner loop above could be computed in
  ** sln->backNum and sln->forwNum
  **
  ** No further error checking (or biff usage) is needed.
  */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  sln->forwNum = 0;
  sln->backNum = 0;
  sln->pos[2 * sln->halfLen] = seedX;
  sln->pos[2 * sln->halfLen + 1] = seedY;  
  
  uint counter = 0;
  while (counter < sln->halfLen) {
    vcrCalcDpos(sln->pos + 2*(sln->halfLen+counter+1), sln->pos + 2*(sln->halfLen+counter), hh, intg, normalize, 1, ctx);   
    real newpos[2];
    V2_SET(newpos, sln->pos[2*(sln->halfLen+counter+1)], sln->pos[2*(sln->halfLen+counter+1) + 1]);
    vcrConvoEval(ctx, newpos[0], newpos[1], 0);
    if (!ctx->inside)
        break;
    counter++;
  }
  sln->forwNum = counter;

  uint counter2 = 0;
  while (counter2 < sln->halfLen) {
    vcrCalcDpos(sln->pos + 2*(sln->halfLen - counter2-1), sln->pos + 2*(sln->halfLen-counter2), hh, intg, normalize, -1, ctx);  
    real newpos2[2];
    V2_SET(newpos2, sln->pos[2*(sln->halfLen -counter2 - 1)], sln->pos[2*(sln->halfLen - counter2 - 1) + 1]);
    vcrConvoEval(ctx, newpos2[0], newpos2[1], 0);
    if (!ctx->inside)
        break;     
    counter2++;
  }
  sln->backNum = counter2;  

  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code (59L in ref)

  return 0;
}
