#include "rnd.h"
#include "rndPrivate.h"

// for any new functions you write to be used by the rndRay functions below
// v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
real lerp(real omin, real omax, real imin, real xx, real imax) {
  real alpha = (xx-imin)/(imax-imin);
  return (1-alpha)*omin+alpha*omax;
}
real max(real a, real b) {return a>b?a:b; }
// ( 5 lines in reference implementation )
// ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

/*
******** rndBlinnPhong: compute Blinn-Phong lighting
**
** INPUT:
** real rgbaIn[4]: the intrinsic "material" RGBA color, before any lighting
**                 (and before any depth-cueing)
** real grad[3]: the world-space gradient
** real Vdir[3]: the direction towards the viewer (e.g. rndCtx->cam.n)
** rndCtx *ctx: contains needed lighting parameters
** This function is not passed a rndRay pointer (so maybe it doesn't belong
** in this file), but it will mainly be called by rndRaySample.
**
** OUTPUT:
** real rgbaOut[4]: the results of computing Blinn-Phong on rgbaIn and grad.
** rgbaOut[3] is copied from rgbaIn[3]. Do *not* compute depth-cueing here.
**
** You will need the basic Blinn-Phong parameters in
** ctx->lparm.{ka,kd,ks,p}. If !(ctx->light.num) (there are no directional
** lights) or both kd and ks are zero, then only the ambient term remains:
** rgbaOut = ka*rgbaIn.  Otherwise you need the "surface normal" to compute
** the lighting. Find the surface normal as the *negative* gradient,
** *normalized*.  If the gradient length is zero, the diffuse and specular
** terms are zero.  Else (with a meaningful surface normal), loop through the
** ctx->light.num lights described by the ctx->light struct, and add the
** diffuse and specular term for each.  For light with index ii, you need:
**
**    const real *lcol = ctx->light.rgb + 3*ii;
**    const real *ldir = ctx->light.xyz + 3*ii;
**
** lcol is the light RGB color; ldir is the *world-space* unit-length
** direction *towards* the light. Besides .num, .rgb, and .xyz, no
** other members of ctx->light are used for this function.
*/
void
rndBlinnPhong(real rgbaOut[4], const real rgbaIn[4],
              const real grad[3], const real Vdir[3],
              const rndCtx *ctx) {
  V4_COPY(rgbaOut, rgbaIn);
  if (!rgbaIn[3]) {
    /* no opacity here, so no reason to compute a color that won't matter.
       This assumes blending with rndBlendOver; results may unexpected
       for other blendings, but oh well */
    return;
  }
  // else, you must to over-write rgbaOut[0,1,2] with your lighting results
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  rgbaOut[3] = rgbaIn[3];
   if(!ctx->light.num || (ctx->lparm.kd == 0 && ctx->lparm.ks == 0 )||V3_LEN(grad)==0)
         V3_SCALE(rgbaOut, ctx->lparm.ka, rgbaIn);
    else{
        real ambient[3],sum[3];
        V3_SCALE(ambient,ctx->lparm.ka,rgbaIn);
        V3_COPY(sum,ambient);
        
        real N[3];  
        V3_SCALE(N, -1, grad);
        real tmp;
        V3_NORM(N, N, tmp);
        for(int i=0;i<(signed)ctx->light.num;i++){
            const real *lcol = ctx->light.rgb + 3*i;
            const real *ldir = ctx->light.xyz + 3*i;
            real CL[3],L[3];
            V3_SET(CL,lcol[0],lcol[1],lcol[2]);
            V3_SET(L,ldir[0],ldir[1],ldir[2]);
            
            real diffuse[3];
            real dcoeff = ctx->lparm.kd*max(0,V3_DOT(L, N));
            V3_SET(diffuse,dcoeff*CL[0]*rgbaIn[0],dcoeff*CL[1]*rgbaIn[1],dcoeff*CL[2]*rgbaIn[2]);
            
            real H[3];
            V3_ADD(H,Vdir,L);
            V3_NORM(H,H,tmp);
            real specular[3];
            real scoeff = ctx->lparm.ks*pow(max(0,V3_DOT(H, N)),ctx->lparm.p);
            V3_SET(specular,scoeff*CL[0],scoeff*CL[1],scoeff*CL[2]);
            
            V3_ADD(sum,sum,diffuse);
            V3_ADD(sum,sum,specular);
        }
          rgbaOut[0] = sum[0];
          rgbaOut[1] = sum[1];
          rgbaOut[2] = sum[2];
    }
  // ( 39 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return;
}

/*
** rndRayStart: This is called to start the ray that will render pixel
** (hi,vi) (faster,slower) in the image array.  It must determine where, in
** world-space, the ray starts (specifically, where on the near clipping
** plane), and how to find subsequent samples on the ray.  Remember that the
** geometry of the sampling is as a stack of planes, separated by
** ctx->planeSep.  If ctx->cam.ortho, all ray samples are separated by
** ctx->planeSep, but the spacing between samples will in general be *larger*
** than ctx->planeSep with perspective projection.  Still, each ray will have
** a fixed step size between samples, and a fixed direction.
**
** The reference implementation also uses this function to set up a pointer
** (inside "ray") to point to probing results (also inside "ray") that will
** later be blended, so that the blending code doesn't need a switch on
** ctx->probe.
*/
void
rndRayStart(rndRay *ray, uint hi, uint vi,
            rndConvo *cnv, const rndCtx *ctx) {
  int codeMissing=0; /* this is set to 1 in the given student code so that if
                        nothing is implemented, the ray will finish quickly */
  if (ctx->timing) {
    ray->time0 = airTime();
  }
  ray->hi = hi;
  ray->vi = vi;
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  // Ray start location
  real h_view = lerp(-ctx->cam.wdth/2.0,ctx->cam.wdth/2.0,-0.5,hi,ctx->cam.size[0]-0.5);
  real v_view = lerp(ctx->cam.hght/2.0,-ctx->cam.hght/2.0,-0.5,vi,ctx->cam.size[1]-0.5);

  real ncu_view = ctx->cam.ncv/ctx->cam.d*h_view;
  real ncv_view = ctx->cam.ncv/ctx->cam.d*v_view;
  real ncn_view = -ctx->cam.ncv;
  V3_SET(ray->uvn,ncu_view,ncv_view,ncn_view);
  real nc_view[4],nc_world[4];
  V4_SET(nc_view,ncu_view,ncv_view,ncn_view,1);
  MV4_MUL(nc_world,ctx->cam.VtoW,nc_view);
  ray->startx = nc_world[0];
  ray->starty = nc_world[1];
  ray->startz = nc_world[2];
  // ray step size
  ray->stepu = h_view/ctx->cam.ncv*ctx->planeSep;
  ray->stepv = v_view/ctx->cam.ncv*ctx->planeSep;
  if(ctx->cam.ortho){
    ray->stepu = ray->stepv = 0.0;
  }
  ray->stepn = -ctx->planeSep;
 
  ray->step_size = sqrt(ray->stepu*ray->stepu+ray->stepv*ray->stepv+ray->stepn*ray->stepn);

  ray->num_sample = ray->num_inside = 0;

  ray->initialized = 0;

  V4_SET(ray->current,0,0,0,0);

  // ( 76 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  rndRaySample(ray, cnv, ctx);
  rndRayBlendInit(ray, cnv, ctx);
  // if planeSep is zero, we're already done
  if (codeMissing || !ctx->planeSep) {
    rndRayFinish(ray, ctx);
  } else {
    ray->finished = 0;
  }
  return;
}

/*
** rndRaySample: This is the function that does all the work of learning and
** computing the probe at the current ray position.  Whatever calls this
** function is going to then call rndRayBlendInit (first time) or rndRayBlend
** (after that), so this function needs to compute and store the results of
** sampling the volume at the current ray position, and then, prepare the
** results to be blended.
**
** Because this has to handle the computation of all the different probes, it
** would be natural to have "if"s or "switch"s around ctx->probe, and to
** structure the code in a way that exploits how the numerical value of the
** probes is a topological sort of their dependency relationship.  But, there
** shouldn't be code specific to ctx->blend here; that should be in rndBlend.
**
** For rndProbeRgba: The opacity from either the univariate lookup table, or
** the Levoy opacity functions, may be outside the range [0,1] (depending on
** how those opacity functions have been set up). Be sure to clamp the
** opacity to range [0,1] prior to doing opacity correction based on step
** size, and blending.
**
** For rndProbeRgbaLit: if there is no opacity, there is no need to compute
** lighting.  If there is opacity, then use ctx->light and ctx->lparm to
** compute Blinn-Phong lighting, and then depth-cueing.  For Blinn-Phong,
** remember that if g is the gradient, then the "normal" to use is -g/|g|,
** the negative normalized gradient.  If |g|=0, then set n=0 to zero out all
** directional lighting effects.
*/
void
rndRaySample(rndRay *ray, rndConvo *cnv, const rndCtx *ctx) {
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  V4_SET(ray->current,0,0,0,0);
    if(ctx->probe<=rndProbeInside)
      ray->num_sample++;
    if(ctx->probe == rndProbeUnknown)
       return;
    V4_SET(ray->current,ray->uvn[0],ray->uvn[1],ray->uvn[2],1);
    if(ctx->probe == rndProbePosView)
        return;

    real view[4],world[4];
    V4_SET(view,ray->uvn[0],ray->uvn[1],ray->uvn[2],1);
    MV4_MUL(world, ctx->cam.VtoW, view);
    V3_SET(ray->xyz_w,world[0],world[1],world[2]);
    V4_SET(ray->current,world[0],world[1],world[2],1);
    if(ctx->probe == rndProbePosWorld)
        return;
    
    real index[4];
    V4_SET(world,ray->xyz_w[0],ray->xyz_w[1],ray->xyz_w[2],1);
    MV4_MUL(index, ctx->WtoI, world);
    V3_SET(ray->xyz_i,index[0],index[1],index[2]);
    V4_SET(ray->current,index[0],index[1],index[2],1);
    if(ctx->probe == rndProbePosIndex)
         return;
    rndConvoEval(ray->xyz_w[0],ray->xyz_w[1],ray->xyz_w[2],cnv,ctx);
    ray->inside = cnv->inside;
    ray->current[0] = cnv->inside;
    if(ctx->probe == rndProbeInside)
         return;
    if(cnv->inside){
      ray->num_inside++;
    ray->convo_val = cnv->value;
    ray->current[0] = cnv->value;
    if(ctx->probe == rndProbeValue)
         return;
    V3_SET(ray->convo_grad,cnv->gradient[0],cnv->gradient[1],cnv->gradient[2]);
    V4_SET(ray->current,cnv->gradient[0],cnv->gradient[1],cnv->gradient[2],0); 
    if(ctx->probe == rndProbeGradVec)
        return; 
    
    ray->convo_gradm = V3_LEN(ray->convo_grad);
    ray->current[0] = ray->convo_gradm;
    if(ctx->probe == rndProbeGradMag)
         return;
    
    real convo_val = ray->convo_val;
    if(convo_val<ctx->txf.vmin)
       convo_val = ctx->txf.vmin;
    else if(convo_val>ctx->txf.vmax)
       convo_val = ctx->txf.vmax;
    real interval = (ctx->txf.vmax - ctx->txf.vmin)/ctx->txf.len;
    int index = floor((convo_val - ctx->txf.vmin)/interval);
    real a = ctx->txf.rgba[index*4+3];
    if(ctx->levoy.num)
      a = a*rndTxfLevoy(ctx->levoy.vra,ctx->levoy.num,ray->convo_val,ray->convo_gradm);
    if(a<0)
      a = 0;
    else if(a>1)
      a = 1;
    //opacity correction
    a = 1-pow((1-a),ray->step_size/ctx->txf.unitStep);

    V4_SET(ray->txf_rgba,ctx->txf.rgba[index*4],ctx->txf.rgba[index*4+1],ctx->txf.rgba[index*4+2],a);
    V4_SET(ray->current,ctx->txf.rgba[index*4],ctx->txf.rgba[index*4+1],ctx->txf.rgba[index*4+2],a);
    if(ctx->probe == rndProbeRgba)
      return;

    rndBlinnPhong(ray->lit_rgba,ray->txf_rgba,ray->convo_grad,ctx->cam.n,ctx);
    real coeff = (-ctx->cam.ncv-ray->uvn[2])/(-ctx->cam.nc+ctx->cam.fc);
    real lerp_coeff[3];
    V3_LERP(lerp_coeff, ctx->lparm.dcn, ctx->lparm.dcf, coeff);
    ray->lit_rgba[0]*=lerp_coeff[0];
    ray->lit_rgba[1]*=lerp_coeff[1];
    ray->lit_rgba[2]*=lerp_coeff[2];
    V4_COPY(ray->current,ray->lit_rgba);
   
    if (ctx->probe == rndProbeRgbaLit)
       return;
    }
  // ( 118 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return;
}

/*
** rndRayBlendInit: This called by rndRayStart after it calls rndRaySample for
** the first time (for that ray); so this needs to initialize whatever buffer
** is used for accumulating the blending computations.  Conversely if the
** probe couldn't be computed (e.g. probing rndProbeValue but ray is outside
** the volume), then this function needs to remember that the blending buffer
** has yet to be initialized.
*/
void
rndRayBlendInit(rndRay *ray, const rndConvo *cnv, const rndCtx *ctx) {
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  if (ctx->probe <= rndProbeInside){
    V4_COPY(ray->accum,ray->current);
    return;
  }
  if(cnv->inside){
    V4_COPY(ray->accum,ray->current);
    ray->initialized = 1;
    return;
  }
  // ( 20 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return;
}

/*
** rndRayBlend: This is called after rndRaySample, to blend (according to
** ctx->blend) the results of latest sample with whatever has already been
** blended/accumulated in this ray.
**
** You might want to use a conditional test like:
**
**    if (ctx->probe <= rndProbeInside || cnv->inside) { ... }
**
** to control whether you do any work here (the same consideration applies to
** rndRayBlendInit): if the probe is rndProbeInside, or anything smaller
** (probes about ray position), then there is always something to blend.
** Otherwise, the probe depends on the convolution result, and there is only
** something to blend if the convolution result was good, that is, if
** cnv->inside is true.
**
** Other than the suggested usage of ctx->probe above, the code here should
** be organized in terms of ctx->blend.  The fine-grained dependence on
** ctx->probe should be in rndRaySample, and this function handles the
** (mostly) orthogonal issue of blending.
**
** Calling rndRayBlendInit from here makes sense, e.g. for rndProbeValue, when
** you detect that this is the first time on the ray that cnv->inside is true.
*/
void
rndRayBlend(rndRay *ray, const rndConvo *cnv, const rndCtx *ctx) {
  /* If blending either rndProbeRgba or rndProbeRgbaLit with rndBlendOver,
     and if opacity has exceeded ctx->txf.alphaNear1, then set opacity to
     exactly 1.0 and call rndRayFinish() */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  switch(ctx->blend){
        case rndBlendUnknown:
         ;
         break;
        case rndBlendMax:
         if (ctx->probe <= rndProbeInside) {
              for(int i=0;i<(signed)rndProbeLen(ctx->probe);i++)
                ray->accum[i] = max(ray->accum[i],ray->current[i]); 
         }
         else{
               if(cnv->inside){
                  if(ray->initialized)
                    for(int i=0;i<(signed)rndProbeLen(ctx->probe);i++)
                      ray->accum[i] = max(ray->accum[i],ray->current[i]);
                  else{
                    V4_COPY(ray->accum,ray->current);
                    ray->initialized = 1;
                  }
               }
         }
         break;
         case rndBlendSum:
         case rndBlendMean:
         if (ctx->probe <= rndProbeInside) {
              for(int i=0;i<(signed)rndProbeLen(ctx->probe);i++)
                ray->accum[i] = ray->accum[i]+ray->current[i]; 
         }
         else {
               if(cnv->inside){
                  if(ray->initialized)
                    for(int i=0;i<(signed)rndProbeLen(ctx->probe);i++)
                      ray->accum[i] = ray->accum[i]+ray->current[i];
                  else{
                    V4_COPY(ray->accum,ray->current);
                    ray->initialized = 1;
                  }
               }
         }
         break;
         case rndBlendOver:
             if(ctx->probe == rndProbeRgba || ctx->probe == rndProbeRgbaLit){
                if(ray->accum[3]>ctx->txf.alphaNear1){
                  ray->accum[3]=1.0;
                  rndRayFinish(ray,ctx);
                 // return;
                }
                if(cnv->inside){
                  if(ray->initialized){
                    real alpha0 = ray->accum[3];
                    real alpha1 = ray->current[3];
                    real alpha = alpha0+(1-alpha0)*alpha1;
                    if(alpha==0.0)
                      return;
                    else{
                    ray->accum[3] = alpha;
                    for(int i=0;i<3;i++)
                      ray->accum[i] = (alpha0*ray->accum[i]+(1-alpha0)*alpha1*ray->current[i])/ray->accum[3];
                    }
                  }
                  else{
                    V4_COPY(ray->accum,ray->current);
                    ray->initialized = 1;
                  }
                }
              }
          break;

      }    
  // ( 38 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return;
}

/*
** rndRayStep: Compute the next step along the ray: increment the ray position
** (but if has just stepped past the far clipping plane, call
** rndRayFinish()), and do whatever adjustments to ray state are needed prior
** to calling rndRaySample() and rndRayBlend().
*/
void
rndRayStep(rndRay *ray, rndConvo *cnv, const rndCtx *ctx) {
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  ray->uvn[0]+=ray->stepu;
  ray->uvn[1]+=ray->stepv;
  ray->uvn[2]+=ray->stepn;
  if(ray->uvn[2] < -(ctx->cam.d+ctx->cam.fc)){
    rndRayFinish(ray,ctx);
   // return;
  }
  // ( 15 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  rndRaySample(ray, cnv, ctx);
  rndRayBlend(ray, cnv, ctx);
  return;
}

/*
** rndRayFinish: Finish the computation of a ray. This entails doing any
** final computation of blending results (e.g. computing rndBlendMean by
** dividing the sum by the number of samples, or dividing out opacity from
** colors with pre-multiplied alpha for rndBlendOver).
**
** Finish here by setting rndProbeLen(ctx->probe) values in ray->result[];
** these will be picked up by rndRender().  Use ctx->outsideValue when there
** is no valid answer.
**
** Note that this can be called from rndRayStart (when 0==ctx->planeSep),
** rndRayStep (when the ray was gone past the far clipping plane), and
** rndRayBlend (when opacity has exceeded ctx->alphaNear1, for rndBlendOver
** blending of rndProbeRgba or rndProbeRgbaLit)
*/
void
rndRayFinish(rndRay *ray, const rndCtx *ctx) {
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  if(ctx->probe > rndProbeInside && !ray->initialized)
      ray->result[0] = ctx->outsideValue;
    else if(ctx->blend==rndBlendMean){
      int denom = ctx->probe <= rndProbeInside?ray->num_sample:ray->num_inside;
      for(int i=0;i<(signed)rndProbeLen(ctx->probe);i++)
        ray->result[i] = ray->accum[i]/denom;
    } else{
      for(int i=0;i<(signed)rndProbeLen(ctx->probe);i++)
        ray->result[i] = ray->accum[i];

    }
  // ( 42 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  if (ctx->timing) {
    ray->time1 = airTime();
    ray->time = ray->time1 - ray->time0;
    // use code block below to fix ray->time for rays computed simultaneously
  }
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  // ( 0 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  ray->finished = 1;
  return;
}
