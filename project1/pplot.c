#include "plt.h"

/*
** rgb_e enum: the different colors to use in the plot image
*/
typedef enum {
  rgbBackground,
  rgbAxes,
  rgbPolynomial,
  rgbOutsideMore,
  rgbOutsideOne,
  rgbPlot,
  rgbDot,
  rgbZeroCrossing,
} rgb_e;

static unsigned char RGBLUT[][3] = {
  {255, 255, 255}, // rgbBackground
  {220, 220, 255}, // rgbAxes
  {255, 150, 110}, // rgbPolynomial
  {20,  20,  255}, // rgbOutsideMore
  {220, 100, 255}, // rgbOutsideOne
  {0, 0, 0},       // rgbPlot
  {100, 100, 100}, // rgbDot
  {100, 255, 100}, // rgbZeroCrossing
};

/*
** You should use this function to set a given pixel in the output image, e.g.
** setrgb(pix, rgbBackground) or setrgb(pix, rgbPlot).  It takes one of the
** rgb_e enum values to lookup into the RGBLUT array; changes to the RGBLUT
** contents should change the output colors your code uses.
*/
static void
setrgb(unsigned char *pix, rgb_e which) {
  const unsigned char *rgb = RGBLUT[which];
  pix[0] = rgb[0];
  pix[1] = rgb[1];
  pix[2] = rgb[2];
}

/*
******** pltPlot: plots a convolution result into an RGB image,
** according to various parameters.
**
** INPUT:
** xsize, ysize: the output image array rgb has been pre-allocated for
** 3*xsize*ysize pixels (each of type unsigned char).  The axis ordering
** of array "rgb" is (fast to slow) is RGB, X, Y.
**
** xmin, xmax: the X (horizontal) axis of the output image should have xsize
** CELL-CENTERED samples in the world-space interval [xmin,xmax], with xmin on
** the left side of the image, and xmax on the right side. Note that RGB
** raster images are laid out with the X index increasing from left to right.
**
** ymin, ymax: the Y (vertical) axis of the output image should have ysize
** CELL-CENTERED samples in the world-space interval [ymin,ymax]. ymin is the
** vertical position of the *bottom* of the image; ymax is the vertical
** position of the *top* of the image. Because RGB raster images are laid out
** with the Y index increasing from top to bottom, increasing the Y
** (world-space) coordinate will decrease the (index-space) coordinate, which
** is opposite of how the X axis is handled.
**
** deriv: the deriv-th derivative of the convolution result should be
** plotted.  With deriv==0, the convolution itself is plotted, with deriv==1,
** the first derivative is plotted, etc.  Derivatives are with respect to
** *world* space (thanks to pltConvoEval).
**
** idata, kern: exactly as in pltConvoEval.  If idata is non-NULL, the
** convolution of idata with kern is plotted. The convolution result should
** ONLY be plotted where it is defined (as per the *outside value set by
** pltConvoEval).  If idata is NULL, then there is no data convolution to
** compute or plot; perhaps only the polynomial is to be plotted.
**
** pc, pclen: a representation of a polynomial, as with pltPolyEval. If pc is
** non-NULL and pclen > 0, then the polynomial represented by pc should also
** be plotted.  This analytic polynomial is defined for all world-space
** positions, regardless of the sampling interval of the input data. If pc is
** NULL, then no polynomial should be plotted.
**
** thicknessPlot: the *vertical* height, in pixels, of the graph of the
** convolution (and if requested, the polynomial). That is, the judgement of
** whether a pixel is inside the plot curve is based on the difference between the
** vertical (index-space) pixel location and the vertical (index-space)
** position of the plot.
**
** apcoth: (if non-zero (true), then APproximate COnstant-THickness plots of
** the convolution (and if requested, the polynomial). The per-pixel inside
** test is still based on vertical offset from the plot, but by knowing the
** slope (in index-space) of the plot at that horizontal position, and using
** the Pythagorean Theorem, the amount of the vertical offset that still
** counts as inside the plot curve becomes a function of slope: the more plot
** curve is near vertical (with either positive or negative slope), the
** greater the range of vertical offsets that are inside the plot curve.
**
** thicknessAxes: if this is non-zero: draw world-space axes for y=0 and x=0,
** with this thickness (in pixels).
**
** showOutside: if this is non-zero and there is non-NULL idata: color the
** portions of the y=0 axis that are within the sampled interval (accounting
** for the cell- vs node-centered-ness of the interval), but outside where
** you can do convolution because the kernel support extends past valid data
** indices.  If the kernel was missing more than one data value, color the
** axis by rgbOutsideMore; if only one data value was missing, color the axis
** by rgbOutsideOne.  Even when so colored, the y=0 axis still has thickness
** thicknessAxes. To implement this, note how pltConvoSample() uses pltNan()
** to produce a NaN when some samples are missing from the convolution sum.
** The NaN stores an integer "payload": the number of missing samples. The
** payload can be retrieved with pltNanPayload(), which determines whether to
** use rgbOutsideOne or rgbOutsideMore.
**
** diameterDot: if this is non-zero: indicate the input discrete data values
** with circular dots, with this diameter (in pixels).  There are smart and
** dumb ways of implementing this; dumb is great for this project.  In
** particular, it is okay if the execution of pltPlot slows down with
** increasing number of data points.
**
** heightZC: if this is non-zero, and if idata is non-NULL, draw vertical
** tickmarks of this height (in pixels) and with thicknessAxes, at those
** locations along the horizontal axis where linear interpolation (lerping) of
** successive input data values produced a zero-crossing (i.e. y=0).  This is
** independent of node- vs cell-centering, and of the kernel indicated.
**
** OUTPUT:
** rgb: the output plot image is stored here, with graphical marks
** (i.e. plots, axes, dots, tick-marks) drawn according to the parameters and
** specification above. A pixel is inside a graphical mark if the *center* of
** the pixel (recall that this image space is always sampled with
** cell-centered samples) falls within the image-space extent of that
** mark. For all the graphical marks, the *location* of the mark is defined in
** terms of world-space, but the size (thickness or diameter), and hence a
** pixel's membership in the mark, is defined in terms of image index-space.
**
** For each of the marks, there is a pre-defined RGB color to use to indicate
** that element, as follows:
**
** -- rgbBackground: color of the background (no marks present at that pixel)
** -- rgbAxes: for the X and Y axes through the origin
**        (y=0 and x=0, respectively)
** -- rgbPolynomial: for the plot of the polynomial (if requested)
** -- rgbOutsideMore: show on X (y=0) axis that X location is inside data
**        sampling interval but outside convolution support by  more than
**        one sample
** -- rgbOutsideOne: show on X axis that X location is inside data sampling
**        interval, but outside convolution support by just one sample
** -- rgbPlot: plot of result of convolving sampled data with kernel
** -- rgbDot: dots to show the input data values
** -- rgbZeroCrossing: vertical tick-marks at zero-crossings of linearly-
**        interpolated data
**
** The ordering is important: the correct image is defined as the result of
** the marks being drawn in the above order: first background, then axes,
** then polynomial, etc.  In other words, the zero-crossing tick-marks appear
** on top of any other marks. The dots to show input data values appear on
** top of the convolution plot, which appears on top of the polynomial plot,
** etc.  The ordering of marks above is also very roughly a possible order of
** implementation: the grading will start with simpler plots using only the
** earlier marks, and then advance to more complex plots using all the marks
** (such as the zero-crossing indicators).
**
** RETURNS:
** 0 if all is well, else
** 1 in case of error (using biff), because of missing pointers or
** nonfinite xw.  No additional error checking is required in the
** student code.
*/
int
pltPlot(unsigned char *rgb,
        uint xsize, uint ysize,
        float xmin, float xmax, float ymin, float ymax,
        uint deriv,
        const pltData *idata, const pltKernel *kern,
        const float *pc, uint pclen,
        float thicknessPlot, int apcoth,
        float thicknessAxes, int showOutside,
        float diameterDot, float heightZC) {
  if (!(rgb)) {
    biffAddf(PLT, "%s: got NULL rgb pointer", __func__);
    return 1;
  }
  /* idata and kern are allowed to be NULL: no plot of data or convolution
     requested, perhaps we plot just the polynomial */
  if (pc && !pclen) {
    biffAddf(PLT, "%s: got zero polynomial coeffs?", __func__);
    return 1;
  }
  /* either pc is NULL (no plot of polynomial) or
     we have 1 or more coefficients in polynomial */
  if (!(isfinite(thicknessPlot) && thicknessPlot >= 0 &&
        isfinite(thicknessAxes) && thicknessAxes >= 0 &&
        isfinite(diameterDot) && diameterDot >= 0)) {
    biffAddf(PLT, "%s: non-finite or negative graphing params (%g,%g,%g)",
             __func__, thicknessPlot, thicknessAxes, diameterDot);
    return 1;
  }
  /* NOTE: the ydata->vv[] data elements are not actually needed or set; we
     just use "ydata" to store the meta-data about the Y axis, so that you can
     use pltItoW and pltWtoI to convert between Y index- and world-space.  You
     may lose style points if you don't use ydata this way. Note also the
     order in which ymin and ymax are passed to pltDataInit: this is how we
     reconcile the raster image convention of Y index increasing as we move
     downwards, with the Cartesian coordinate convention of Y position
     decreasing as we move downwards. */
  airArray *mop = airMopNew();
  pltData *ydata = pltDataNew();
  airMopAdd(mop, ydata, (airMopper)pltDataNix, airMopAlways);
  if (pltDataInit(ydata, ysize, ymax, ymin, pltCenterCell)) {
    biffAddf(PLT, "%s: trouble setting up ydata", __func__);
    airMopError(mop);
    return 1;
  }
  pltData *xdata = pltDataNew();
  airMopAdd(mop, xdata, (airMopper)pltDataNix, airMopAlways);
  pltData *ddata;
  if (idata) {
    // we have data for which we plot the convolution
    if (!kern) {
      biffAddf(PLT, "%s: got NULL kernel", __func__);
      airMopError(mop);
      return 1;
    }
    ddata = pltDataNew();
    airMopAdd(mop, ddata, (airMopper)pltDataNix, airMopAlways);
    /* precompute a vector of convolution results, one value per horizontal
       output pixel, both for what we plot (the "deriv"-th derivative ) and
       its derivative (the ("deriv"+1)-th derivative) */
    if (pltConvoSample(xdata, xsize, xmin, xmax, pltCenterCell, deriv,
                       idata, kern) ||
        pltConvoSample(ddata, xsize, xmin, xmax, pltCenterCell, deriv+1,
                       idata, kern)) {
      biffAddf(PLT, "%s: trouble convolving with input data", __func__);
      airMopError(mop);
      return 1;
    }
  } else {
    // even if not plotting data, we can usefully set fields in xdata
    if (pltDataInit(xdata, xsize, xmin, xmax, pltCenterCell)) {
      biffAddf(PLT, "%s: trouble setting up xdata", __func__);
      airMopError(mop);
      return 1;
    }
    ddata = NULL;
  }
  /* either way, now xdata is set: YOU SHOULD USE IT as part of any and all
     conversions between index and world space on the X axis. You may lose
     style points if you don't. */
  pltData *pdata, *qdata; // for polynomial (pdata) & its derivative (qdata)
  if (pc) {
    // we have a polynomial to plot
    pdata = pltDataNew();
    airMopAdd(mop, pdata, (airMopper)pltDataNix, airMopAlways);
    if (pltPolySample(pdata, xsize, xmin, xmax, pltCenterCell, pc, pclen)) {
      biffAddf(PLT, "%s: trouble sampling polynomial", __func__);
      airMopError(mop);
      return 1;
    }
    if (apcoth) {
      // need the derivative of polynomial
      uint qclen = pclen > 1 ? pclen-1 : 1;
      // set up qc as coefficients of the derivative of polynomial
      float *qc = PLT_MALLOC(qclen, float); assert(qc);
      airMopAdd(mop, qc, airFree, airMopAlways);
      qc[0] = 0; // in case pclen==1
      for (uint pi=1; pi<pclen; pi++) {
        qc[pi-1] = pi*pc[pi];
      }
      qdata = pltDataNew();
      airMopAdd(mop, qdata, (airMopper)pltDataNix, airMopAlways);
      // sampled the polynomial derivative
      if (pltPolySample(qdata, xsize, xmin, xmax, pltCenterCell, qc, qclen)) {
        biffAddf(PLT, "%s: trouble sampling polynomial derivative", __func__);
        airMopError(mop);
        return 1;
      }
    } else {
      qdata = NULL;
    }
  } else {
    // no polynomial to plot
    pdata = NULL;
    qdata = NULL;
  }

  /*
  ** Set very first pixel to background color, to disambiguate (by looking at
  ** the output PNG) where in the image the first pixel is.
  */
  setrgb(rgb, rgbBackground);
  /* the beginning of the following student code block is where you
     can do additional experiments with setrgb() to make sure you
     understand the pointer arithmetic required to set pixels' colors,
     until you are able to write the solution code */

  /*
  ** Suggested strategy: traverse the output image index space. At each
  ** pixel, initialize its color to the background, and then find its
  ** world-space location. Use a mix of geometry and arithmetic to see if
  ** that location is inside one of the things to draw: coordinate axes, the
  ** polynomial, dots at data points, convolution result, or zero-crossing
  ** tick marks (in that order), and set the pixel color accordingly. Make
  ** use of the pre-computed (above) samples of convolution (xdata, and its
  ** derivative ddata) or polynomial (pdata, and its derivative qdata): no
  ** new work for doing convolution or polynomial evaluation should be
  ** happening below. It is up to you to do the pointer arithmetic needed to
  ** get the address of a specific pixel (other than the first pixel) to pass
  ** to setrgb().
  **
  ** Implementing this function is essentially a test of one's ability to
  ** juggle the X and Y image world-spaces, and X and Y image index-spaces,
  ** and the world- and index-spaces of the input data.
  */
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  
  // the pixel index of origin in world space
  float x_axis_yi = pltWtoI(ydata,0);
  float y_axis_xi = pltWtoI(xdata,0);

  // convert the slope from world space to index space
  // for later use by finding the linear lerping coefficient 
  // in y axis for plotting polynomials and convolution.
  // let yi = a*yw + b => a = (a*1+b) - (a*0+b)
  // or a = pltWtoI(ydata,1) - pltWtoI(ydata,0)
  float a = pltWtoI(ydata,1) - pltWtoI(ydata,0);
  // similarly, for xi = c*xw + d
  // c = pltWtoI(xdata,1) = pltWtoI(xdata,0);
  float c = pltWtoI(xdata,1) - pltWtoI(xdata,0);
  // so to translate slope from world space to index space
  // the coefficient is a/c

  // find all zeros for the convoluting curves and store 
  // their x value in world-space in an array for plotting
  // zero-crossing tick marks
  float zeros[idata->len];
  // count of zeros
  int count = 0;
  for(int i=0;i<(signed)(idata->len)-1;i++){
      float y1 = idata->vv[i];
      float x1 = pltItoW(idata,i);
      float y2 = idata->vv[i+1];
      float x2 = pltItoW(idata,i+1);
      float x =  pltLerp(x2,x1,y2,0,y1);
      if(x>=x1 && x <=x2)
        zeros[count++] = x;
  }

  //traverse the image index space and check if each
  //pixel falls into different graphic marks
  for(int i=0;i<(int)xsize;i++)
    for(int j=0;j<(int)ysize;j++){
      // initialize the color to the background
      setrgb(rgb+(j*xsize+i)*3, rgbBackground); 
       // add coordinate axis
      if(thicknessAxes){
        // determine if the pixel(i,j) falls within 
        // either x or y axis mark 
        if((i>=y_axis_xi-thicknessAxes/2.0 && i<=y_axis_xi+thicknessAxes/2.0)
         ||(j>=x_axis_yi-thicknessAxes/2.0 && j<=x_axis_yi+thicknessAxes/2.0))
             setrgb(rgb+(j*xsize+i)*3,rgbAxes);   
      }
      // add polynomials
      if(pc && pclen>0){  
          float yw = pdata->vv[i]; 
          float yi = pltWtoI(ydata,yw);
          if(!apcoth){
            if((j-yi)<=thicknessPlot/2 && (j-yi)>=-thicknessPlot/2)
              setrgb(rgb+(j*xsize+i)*3,rgbPolynomial);
          }
          else{         
            // use pythagoren theorems to figure out the range
            float side1 = (thicknessPlot/2.0);
            // the ratio a/c is to convert slope from world to index space
            float side2 = (thicknessPlot/2.0)*fabsf(qdata->vv[i]*a/c);
            if((j-yi)*(j-yi)<= side1*side1 + side2*side2)
                     setrgb(rgb+(j*xsize+i)*3,rgbPolynomial);
          } 
      }
      //add convolution result
      if(idata){
           float yw = xdata->vv[i];
           // check if the convolution result is valid
           if (isnan(yw)){
                if(showOutside){
                  rgb_e color;
                  if((signed)pltNanPayload(yw)==1)
                        color = rgbOutsideOne;
                    else if((signed)pltNanPayload(yw)>1)
                        color = rgbOutsideMore;
                  // compute the original sampled data range 
                  float xi_min = pltWtoI(xdata,idata->min);
                  float xi_max = pltWtoI(xdata,idata->max);
                  // color the x axis within the sampled data range
                  // based on the number of indices outside the range
                  if(i>=xi_min && i<=xi_max &&
                     j>=x_axis_yi-thicknessAxes/2.0 && 
                     j<=x_axis_yi+thicknessAxes/2.0)
                     setrgb(rgb+(j*xsize+i)*3,color);
                }
           }
           else{
           float yi = pltWtoI(ydata,yw);
           if(!apcoth){
            if((j-yi)<=thicknessPlot/2 && (j-yi)>=-thicknessPlot/2)
              setrgb(rgb+(j*xsize+i)*3,rgbPlot);
            }
           else{
             // use pythagoren theorems to figure out the range
            float side1 = (thicknessPlot/2.0);
            // the ratio a/c is to convert slope from world to index space
            float side2 = (thicknessPlot/2.0)*fabsf(ddata->vv[i]*a/c);
             if((j-yi)*(j-yi)<= side1*side1 + side2*side2)
                setrgb(rgb+(j*xsize+i)*3,rgbPlot);
            } 
           }       
        }
        // add dots at data point
        if(diameterDot){
             // go through each sampled data point
             for(int k=0;k<(int)idata->len;k++){
                 float xw = pltItoW(idata,k);
                 float xi = pltWtoI(xdata,xw);
                 float yw = idata->vv[k];
                 float yi = pltWtoI(ydata,yw);
                 if((xi-i)*(xi-i)+(yi-j)*(yi-j)<=(diameterDot/2.0)*(diameterDot/2.0))
                    setrgb(rgb+(j*xsize+i)*3,rgbDot);
               }
        }
        // add zero-crossing ticks
        if(heightZC && idata){
             // go through each x coordinate of zeros
             for(int k=0;k<count;k++){
                 float xw = zeros[k];
                 float xi = pltWtoI(xdata,xw); 
                 float yi = pltWtoI(ydata,0);
                 if(i>=xi-thicknessAxes/2.0 && i<=xi+thicknessAxes/2.0
                    && j>=yi-heightZC/2.0 && j<=yi+heightZC/2.0)
                     setrgb(rgb+(j*xsize+i)*3,rgbZeroCrossing); 
             }
        }
    }
  // ( 73 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code

  airMopOkay(mop);
  return 0;
}
