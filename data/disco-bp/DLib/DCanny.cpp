
#include <DCanny.h>
#include <vector>
#include <DGaussianKernel.h>

using namespace std;

template<class T>
_DPlane<unsigned char> DCanny<T>::do_canny(const _DPlane<T> &img_in, bool relative, float thresh1, float thresh2, float sigma, int filter_size)
{
  _DPlane<float> img(img_in.rows(), img_in.cols());
  int pixel_count = img_in.rows() * img_in.cols();

  const T *img_in_cp = img_in[0];
  float *img_cp = img[0];
  for(int i=0; i<pixel_count; i++)
    img_cp[i] = float(img_in_cp[i]);

  _DGaussianKernel<float> _gaussian(sigma, filter_size);

  _DMatrix<float> dog_1 = _gaussian / _gaussian.max();
  _DMatrix<float> dog_2(1, filter_size);

  int center = filter_size / 2;
  for(int j=0; j<filter_size; j++)
    {
      int x = j-center, y=0;
      dog_2[0][j] = -x * exp(-(x*x+y*y) / (2 * sigma * sigma)) / (M_PI * sigma * sigma);
    }

  img = img.cross_correlate_separable(_gaussian, _gaussian.transpose(), false);
  dx = img.cross_correlate_separable(dog_2, dog_1.transpose(), false);
  dy = img.cross_correlate_separable(dog_1, dog_2.transpose(), false);

  canny_edges = _DPlane<unsigned char>(img.rows(), img.cols());
  canny_edges = 0;

  mag = _DPlane<float>(dx.rows(), dx.cols());
  float *_mag = mag[0], *_dx = dx[0], *_dy = dy[0];
  int sz = dx.rows() * dx.cols();
  for(int i=0; i<sz; i++)
    _mag[i] = _dx[i]*_dx[i]+_dy[i]*_dy[i];


  if(thresh1 <= 0.0 && thresh2 <= 0.0)
    {
      float max_val = sqrt(mag.max());

      float inc = max_val / 64.0;
      std::vector<int> hist(64, 0);

      int count = img_in.rows() * img_in.cols();      
      for(int i=0; i<count; i++)
	{
	  hist[int(ceil((float(sqrt(_mag[i])) / inc)-0.5))]++;
	}

      int cum_sum = 0;
      int _i=0;
      for(int i=0; i<64; i++)
	{
	  cum_sum += hist[i];

	  if(cum_sum >= 0.7 * count)
	    {
	      _i = i;
	      break;
	    }
	}

      thresh1 = max_val * (_i+1) / 64.0;
      thresh2 = 0.4 * thresh1;
    }


  if (0 && relative) 
    {                    /* thresholds are relative => get noise */
      float noise;			     /* and scale them by it. */
      
      noise = cannyNoiseEstimate(mag);
      assert(noise > 0.0);

      /*  Noise is actually noise squared times 4.0, but, thresholds need to be squared
       *  since they should be thresholds for the gradient magnitude squared * 4.0.
       *  Be sure to preserve the signs of the thresholds, since thresh1 <= 0
       *  means no thresholding, and thresh2 < 0.0 means only one threshold.
       */
      
      CannyThreshold1 = thresh1*noise;	     /* save values in global vars */
      CannyThreshold2 = thresh2*noise;
      
      thresh1 = thresh1 * thresh1 * noise * ((thresh1 < 0.0) ? -1.0 : 1.0);
      thresh2 = thresh2 * thresh2 * noise * ((thresh2 < 0.0) ? -1.0 : 1.0);
    }
  else 
    {		    /* make thresholds right for grad-mag squared * 4.0 */
      CannyThreshold1 = thresh1;		     /* save values in global vars */
      CannyThreshold2 = thresh2;

      thresh1 = thresh1*thresh1;
      thresh2 = thresh2*thresh2;
      
      //      thresh1 = thresh1 * thresh1  * ((thresh1 < 0.0) ? -1.0 : 1.0); // * 4.0;
      //      thresh2 = thresh2 * thresh2  * ((thresh2 < 0.0) ? -1.0 : 1.0); // * 4.0;
    }



  if (thresh2 <= 0.0) {			     /* only one threshold specified */
    canny_edge(img_in, dx, dy, mag, thresh1, canny_edges);
  }
  else 
    {				     /* two thresholds were specified */
      if (thresh1 > thresh2) 
	{		     /* swap 'em if necessary */
	  float tmp = thresh1;
	  thresh1 = thresh2;
	  thresh2 = tmp;
	}

      canny_edge2(img_in, dx, dy, mag, thresh1, thresh2, canny_edges);
    }

  canny_edges = canny_edges.binary_thin();

  return canny_edges;  
}


/* This is based on Chapter 2 of Harry Voorhees' SM Thesis
 * "Finding Testure Boundaries In Images," AI Tech Rpt 968/
 *
 * The histogram of the gradient magnitude should be a Rayleigh
 * distribution:  f(x) = (x/a)*(e^(-x^2/2a^2)) and F(x) = 1 - e^(-x^2/2a^2)
 * Our noise estimate is the peak of f(x), which occurs when x = a.
 */

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

template<class T>
float DCanny<T>::cannyNoiseEstimate(const _DPlane<float> &mag)     /* given gradient magnitude sqaured */
{					     /* (times 4) of the image, give an */
  std::vector<float> hist(NOISE_BUCKETS, 0.0), hists(NOISE_BUCKETS), histd(NOISE_BUCKETS);
  float tmp, scale, mxval;
  int index;

  int initial;
  float noise, sum, cdf;
  float sigma;
  int i, j;
  float minderiv, maxderiv, deriv;
  int estimate;

  mxval = -1.0;				     /* Get maximum gradient magnitude value */
  for (float *scan = mag[0], *end = scan + mag.total_pixel_count(); scan < end; scan++) {
    if ((tmp = *scan) > mxval)
      mxval = tmp;
  }

  mxval *= NOISE_RANGE * NOISE_RANGE;	     /* max mag-sqrd value to histogram */
  scale  = NOISE_BUCKETS / sqrt(mxval);	     /* scale for histograming mag vals */
  
					     /* now histogram the data */
  for (float *scan = mag[0], *end = scan + mag.total_pixel_count(); scan < end; scan++) {
    if ((tmp = *scan) <= mxval) {	     /* ok to histogram? */
      hist[((int) (sqrt(tmp)*scale))] += 1.0; /* yes .. scale & histogram */
      sum += 1.0;
    }
  }
  
  /* partly kill the peak there *//* maybe use init est to figure how many to kill */
  hist[0] = hist[1] = hist[2] = 0.0;
  
  for (sum = 0.0, initial = 0; initial < NOISE_BUCKETS; initial++) {
    sum += hist[initial];
    }
  
  sum *= 0.40;
  for (cdf = 0.0, initial = 0; initial < NOISE_BUCKETS && cdf < sum; initial++) {
    cdf += hist[initial];
  }


  noise = (initial + 0.5)/scale;
  noise = noise*noise;

  sigma = MIN(NOISE_MAX_SIGMA, MAX(NOISE_SIGMA*initial, NOISE_MIN_SIGMA));

  
  _DGaussianKernel<float> _gaussian(sigma, 5);

  hist[0] = 0.0;
  int len = _gaussian.cols();
  float *gaussian = _gaussian[0];

  for (i = 0; i < NOISE_BUCKETS; i++) {
    sum = gaussian[0] * hist[i];
    for (j = 1; j < len; j++) {
      sum += gaussian[j]*((i-j >= 0 ? hist[i-j] : -hist[j-i]) +
				  (i+j < NOISE_BUCKETS ? hist[i+j] : 0.0));
    }
    hists[i] = sum;
  }


  /* derivative based approach */
  minderiv = maxderiv = 0.0;		     /* compute histogram derivative */
  histd[0] = 0.0;
  for (index = 1; index < NOISE_BUCKETS; index++) {
    deriv = histd[index] = hists[index] - hists[index-1];
    
    if (deriv > maxderiv) maxderiv = deriv;
    if (deriv < minderiv) minderiv = deriv;
  }

  /* now use as estimate the beginning of the largest region of negative
     valued derivative. */
  
  { int nareas = 0, begin;
    typedef struct { int begin; float area; } AREA;
    AREA areas[2];
    float area;
    
    
    for (i = 0; i < NOISE_BUCKETS - 1; i++) {
      if (histd[i] >= 0.0  && histd[i+1] < 0.0) {
	begin = ++i;
	area  = histd[i];
	for ( ; i < NOISE_BUCKETS - 1 && histd[i] < 0.0; i++) {
	  area += histd[i];
	}
	area = -area; i--;
	areas[nareas >= 2
	      ? (areas[0].area < areas[1].area
		 ? 0
		 : 1)
	      : nareas++] = ((AREA) {begin: begin, area: area});
      }
    }

    estimate = (nareas == 2
		? (areas[0].area > areas[1].area
		   ? 0
		   : 1)
		: 0);
    estimate = areas[estimate].begin;
    noise = (estimate + 0.5)/scale;
    noise *= noise;

  }
      
  CannyNoiseEstimate = noise;		     /* save in global var */
  return noise;  
}


template<class T>
_DPlane<float> DCanny<T>::get_edgedir_plane()
{
  _DPlane<float> result(dx.rows(), dx.cols());
  result = 0;

  unsigned char *edge_cp = canny_edges[0];
  float *out_cp = result[0];
  for(int i=0; i<dx.rows(); i++)
    for(int j=0; j<dx.cols(); j++, edge_cp++, out_cp++)
      {
	if(*edge_cp)
	  {
	    float this_dx = dx[i][j], this_dy = dy[i][j];

	    if(this_dy == 0)
	      *out_cp = M_PI / 2.0;
	    else
	      *out_cp = atan(this_dx / this_dy);
	  }
      }

  return result;
}



/******************************* Local Functions *******************************/

template<class T>
int DCanny<T>::NMS(float *magn, const float mag, const float dx, const float dy)
{
  float theta, m1, m2;
  int n, o1, o2;
  float interp, interp2;
					     /* map [-pi/2 .. pi/2] to [0.0 .. 4.0] */
  if (dx == 0.0)			     /* Don't divide by 0.0 */
    theta = 0.0;			     /* 90 degrees */
  else
    theta = 4.0*(0.5 + atanpi(dy/dx));	     /* 4*(.5 + PI*theta) */

  n = (int) theta;			     /* region we're in */
  interp = theta - (float) n;		     /* interpolating factors */
  interp2 = 1.0 - interp;

  o1 = offset1[n];			     /* get offsets & interpolate */
  o2 = offset2[n];
  m1 = *(magn + o1)*interp2 + *(magn + o2)*interp;
  m2 = *(magn - o1)*interp2 + *(magn - o2)*interp;

  return (mag>=m1 && mag>=m2 && m1!=m2);     /* return 1 iff passes NMS */
}


template<class T>
int DCanny<T>::canny_edge(const _DPlane<T> &img, const _DPlane<float> &_dx, const _DPlane<float> &_dy, const _DPlane<float> &_magn,
		       const float thresh, _DPlane<unsigned char> &_output)
{
  int x_m = img.cols(), y_m = img.rows();
  offset1[0] = x_m;   offset2[0] = x_m-1;    /* 225..270 and 45..90 */
  offset1[1] = x_m-1; offset2[1] = -1;	     /* 180..225 and 0..45 */
  offset1[2] = 1;     offset2[2] = x_m+1;    /* 315..360 and 135..180 */
  offset1[3] = x_m+1; offset2[3] = x_m;	     /* 270..315 and 90..135 */ 
  offset1[4] = x_m;   offset2[4] = x_m-1;    /* 225..270 and 45..90 */

  unsigned char *output = _output[0] + x_m + 1;			     /* skip first row */
  float *dx     = _dx[0] + x_m + 1;			     /* leave zero=>no edge */
  float *dy     = _dy[0] + x_m + 1;
  
  float *magn_max, *magn = _magn[0];
  
  for (magn_max = magn + x_m*(y_m - 1) - 1, magn += x_m + 1; magn < magn_max; ) {
    for (float *magn_max_x = magn + x_m - 2; magn < magn_max_x;
	 output++, dx++, dy++, magn++)
    {
      float mag;

      

      if ((mag = *magn) < thresh)	     /* don't even need to do NMS */
	*output = 0;
      else				     /* interpolate gradient & do NMS */
	*output = NMS(magn, mag, *dx, *dy);
    }

    output += 2;			     /* skip over last pixel on this line */
    dx     += 2;			     /* and first pixel on next line */
    dy     += 2;
    magn   += 2;
  }

  return 0;
}


template<class T>
int DCanny<T>::canny_edge2(const _DPlane<T> &img, const _DPlane<float> &_dx, const _DPlane<float> &_dy, 
			   const _DPlane<float> &_magn, const float thresh, const float thresh2,
			   _DPlane<unsigned char> &output)
{
  int x_m = img.cols(), y_m = img.rows();
  offset1[0] = x_m;   offset2[0] = x_m-1;    /* 225..270 and 45..90 */
  offset1[1] = x_m-1; offset2[1] = -1;	     /* 180..225 and 0..45 */
  offset1[2] = 1;     offset2[2] = x_m+1;    /* 315..360 and 135..180 */
  offset1[3] = x_m+1; offset2[3] = x_m;	     /* 270..315 and 90..135 */ 
  offset1[4] = x_m;   offset2[4] = x_m-1;    /* 225..270 and 45..90 */

  unsigned char **stack = new unsigned char *[img.rows() * img.cols()];
  unsigned char **stktop = stack;

  unsigned char *out = output[0] + x_m + 1;			     /* skip over first row */
  float *dx  = _dx[0] + x_m + 1;
  float *dy  = _dy[0] + x_m + 1;
  float * magn = _magn[0];
  float *magn_max;

  for (magn_max = magn + x_m*(y_m - 1) - 1, magn += x_m + 1; magn < magn_max; ) {
    for (float *magn_max_x = magn + x_m - 2; magn < magn_max_x;
	 out++, dx++, dy++, magn++)
    {
      float mag;

      if ((mag = *magn) < thresh)	     /* no possible edge */
	*out = 0;
      else
      {
	if (NMS(magn, mag, *dx, *dy))	     /* check if passes NMS */
	  if (mag >= thresh2) {
	    *out = 1;			     /* definitely have an edge */

	    *(stktop++) = out;		     /* put edge pixel addr on stack */
	  }
	  else
	    *out = 2;			     /* maybe have an edge */
	else
	  *out = 0;			     /* no edge here */
      }
    }

    out  += 2;				     /* skip over last pixel on this line */
    dx   += 2;				     /* and first pixel on next line */
    dy   += 2;
    magn += 2;
  }

  while (stktop > stack) {
    out = *--stktop;

    /* look at neighbors. if a neighbor is 2 make it 1 and add to stack */
    /* continue until stack empty */

    if (*(out -= x_m + 1) == 2)              /* upper-left */
      *out = 1, *(stktop++) = out;
    if (*(++out) == 2)                       /* upper-middle */
      *out = 1, *(stktop++) = out;
    if (*(++out) == 2)                       /* upper-right */
      *out = 1, *(stktop++) = out;
    if (*(out += x_m) == 2)                  /* middle-right */
      *out = 1, *(stktop++) = out;
    if (*(out -= 2) == 2)                    /* middle-left */
      *out = 1, *(stktop++) = out;
    if (*(out += x_m) == 2)                  /* lower-left */
      *out = 1, *(stktop++) = out;
    if (*(++out) == 2)                       /* lower-middle */
      *out = 1, *(stktop++) = out;
    if (*(++out) == 2)                       /* lower-right */
      *out = 1, *(stktop++) = out;
  }

  { 
    unsigned char *outend = output[0] + x_m*y_m; /* get rid of any remaining 2's */
    
    for (out = output[0]; out < outend; out++ )
      if (*out == 2)
	*out = 0;
  }

  delete[] stack;

  return 0;
}


template<class T> 
int DCanny<T>::gradientFull(const _DPlane<float> &img, _DPlane<float> &dx, _DPlane<float> &dy, _DPlane<float> &mag)
{
  float *img_xendp;
  register float x, y;
  register float a, b, c, d;
  register float *img_endp;

  int width  = img.cols(), height = img.rows();

  dx = _DPlane<float>(height, width);
  dy = _DPlane<float>(height, width);
  mag = _DPlane<float>(height, width);

  float *imgp = img[0];
  float *dxp = dx[0], *dyp = dy[0], *magp = mag[0];

  for (img_endp = imgp + width*(height-1) - 1; imgp < img_endp; ) {
    a = *imgp;
    c = *(imgp+width);

    for (img_xendp = imgp + width - 1; imgp < img_xendp; ) { /* this really works! */
      imgp++;

      b = *(imgp);
      d = *(imgp+width);

      a = d - a;
      c = b - c;
      x = a + c;
      y = a - c;

      a = b;
      c = d;

      *(dxp++)  = x;
      *(dyp++)  = y;
      *(magp++) = x*x + y*y;
    }

    imgp++;                                  /* last column gets 0.0 */
    *(dxp++) = *(dyp++) = *(magp++) = 0.0;
  }

  for (img_endp += width; imgp <= img_endp; imgp++) /* last row gets 0.0 */
    *(dxp++) = *(dyp++) = *(magp++) = 0.0;

  return 0;
}


#define DECLARE(x) \
  template class DCanny<x>; 

DECLARE(double)
DECLARE(short)
DECLARE(int)
DECLARE(float)
DECLARE(char)
DECLARE(unsigned char)

