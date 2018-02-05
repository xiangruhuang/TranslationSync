//
//  DLib: A simple image processing library.
//
//  David Crandall, 2003-2005
//  crandall@cs.cornell.edu
//
//  Please do not redistribute this code.
//
//
//
//
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "DistanceTransform.h"
using namespace std;

/* Implementation of fast squared Euclidean distance transform algorithm
   using amortized algorithm for lower envelope of quadratics.  

   For a description see
   www.cs.cornell.edu/~dph/matchalgs/iccv2003-tutorial.pdf */

// hacked for C++ by crandall, 9/2003

template<class T>
DistanceTransform_2D<T>::DistanceTransform_2D(const _DMatrix<T> &sigma, bool _old_algo) : old_algo(_old_algo) 
{
  rotate = true;
  if(sigma.is_diagonal())
    rotate = false;
  
  if(rotate)
    {
      _DMatrix<T> eig_vectors;
      _DMatrix<T> eig_values = sigma.eigen(eig_vectors);
      
      if(eig_vectors[1][1] ==0)
	angle=0;
      else
	angle = -atan(eig_vectors[0][1]/eig_vectors[1][1]);
      
      scaling_x = 1/eig_values[0][1];
      scaling_y = 1/eig_values[0][0];
      
    }
  else
    {
      scaling_x = 1/sigma[1][1];
      scaling_y = 1/sigma[0][0];
    }      
}

template<class T>
inline T square(T n)
{
  return n*n;
}


/* dt helper function */
template<class T>
void dt(T *src, T *dst, int s1, int s2, int d1, int d2, T scale, T *lut, int col_count, int col_off) {
  //  if (d2 >= d1) 
  {
    int d = (d1+d2) >> 1;
    int s = s1;
    for (int p = s1; p <= s2; p++)
      //	if (src[s] + square(s-d) * scale> src[p] + square(p-d) * scale)
      if (src[s] + lut[s-d]> src[p] + lut[p-d])
	s = p;
    //      dst[d*col_count+col_off] = src[s] + square(s-d) * scale;
    dst[d*col_count+col_off] = src[s] + lut[s-d];
      
    if(d-1 >= d1)
      dt(src, dst, s1, s, d1, d-1, scale, lut, col_count, col_off);
    if(d2>=d+1)
      dt(src, dst, s, s2, d+1, d2, scale, lut, col_count, col_off);
  }
}

template<class T>
_DPlane<T> DistanceTransform_2D<T>::do_transform(const _DPlane<T> & in_im, bool need_locations, T scaling_x, T scaling_y)
{
  if(old_algo) 
    {
      _DPlane<T> im = in_im;

      // records closest pixel
      if(need_locations)
	row_locations = _DPlane< int >(in_im.rows(), in_im.cols()), col_locations = _DPlane< int >(in_im.rows(), in_im.cols());

      int width = im.cols();
      int height = im.rows();
      int k;
      float s;
      float sp;
      int x, y;
  
      int *z = (int *)malloc(sizeof(int)*(max(width, height)+1));
      int *v = (int *)malloc(sizeof(int)*(max(width, height)));
      float *vref = (float *)malloc(sizeof(float)*(max(width, height)));
      if ((z == NULL) || (v == NULL) || (vref == NULL)) {
	assert(0);
      }
  
      float *lut1 = new float[max(height,width)*2+1];

      float sx_inv = 1.0/scaling_x;
      for(int i=-width,j=0; i<width; i++, j++) 
	lut1[j] = 1.0/(2.0 * i) * sx_inv;
  
  
      /* do x transform */
      for (y = 0; y < height; y++) {
	k = 0;  /* Number of boundaries between parabolas */
	z[0] = 0;   /* Indexes of locations of boundaries,
		       order by increasing x */
	z[1] = width;    
	v[0] = 0;     /* Indexes of locations of visible parabola bases,
			 ordered by increasing x */

	T *im_y=im[y];
	int *closest_row_y, *closest_col_y;
	if(need_locations)
	  closest_row_y = row_locations[y], closest_col_y = col_locations[y];

	int v_k = v[k], v_k_2 = v[k]*v[k];

	for (x = 1; x < width; x++) 
	  {
	    // # of times here = 1.0
	    do {
	      /* compute Vornoi border: intersection of parabola at x
		 with rightmost currently visible parabola */
	      //	s = ((im[y][x] + x*x) - (im[y][v[k]] + v[k]*v[k])) /
	      //  (2 * (x - v[k]));
	      //	s = ((im_y[x] + scaling_x*x*x) - (im_y[v[k]] + scaling_x*v[k]*v[k])) /
	      //	  (2 * scaling_x * (x - v[k]));
	  
	  
	      // count 1.1
	      s = (im_y[x] - im_y[v_k] + scaling_x*(x*x - v_k_2)) *
		lut1[x-v_k+width]; // * sx_inv;
	      // / (2 * scaling_x * (x - v_k));
	  
	      sp = ceil(s); // floor(ceil(s));  // OPTME
	  
	      /* case one: intersection is to the right of the array, 
		 so this parabola is not visible (and nothing to do) */
	      if (sp >= width)
		{ 
		  // count 0.8
		  break;
		}
	  
	      /* case two: intersection is at larger x than rightmost current
		 intersection, so this parabola is visible on the right (add
		 it to the end) */
	      if (sp > z[k]) 
		{
		  z[k+1] = int(sp);
		  z[k+2] = width;
		  v[k+1] = x;
		  k++; 
		  v_k = v[k];
		  v_k_2 = v_k*v_k;
		  // count 0.2
		  break;
		}
	  
	      /* case three: intersection is at smaller x than the
		 rightmost current intersection, so this parabola hides the
		 rightmost parabola (remove the rightmost parabola, if there
		 are still remaining potentially visible parabolas iterate to
		 complete the addition of this parabola). */
	  
	      if (k == 0) 
		{
		  v[k] = x;
		  v_k = v[k];
		  v_k_2 = v_k*v_k;
		  // count 0.002
		  break;
		} 
	      else 
		{
		  z[k] = width;
		  k--;
		  v_k = v[k];
		  v_k_2 = v_k*v_k;
		  // count 0.09
		}
	    } while (1);
	
	  }
    
	/* compute transform values from visible parabolas */
    
	/* get value of input image at each parabola base */
	for (x = 0; x <= k; x++) 
	  {
	    vref[x] = im_y[v[x]];
	  }
	k = 0;
    
	/* iterate over pixels, calculating value for closest parabola */
	v_k=v[k];
	if(need_locations)
	  for (x = 0; x < width; x++) 
	    {
	      if (x == z[k+1]) k++, v_k=v[k];
	      im_y[x] = vref[k] + (v_k-x)*(v_k-x) * scaling_x;
	  
	      closest_row_y[x] = y, closest_col_y[x] = v_k; 
	    }
	else
	  for (x = 0; x < width; x++) 
	    {
	      if (x == z[k+1]) k++, v_k=v[k];
	      im_y[x] = vref[k] + (v_k-x)*(v_k-x) * scaling_x;
	    }
      }


      /* do y transform - analogous computation in y-direction applied to
	 result of x-transform */

      float sy_inv = 1.0/scaling_y;
      for(int i=-height,j=0; i<height; i++, j++)
	lut1[j] = 1.0/(2.0 * i) * sy_inv;

      _DPlane< int > closest2_row(row_locations), closest2_col(col_locations);

      for (x = 0; x < width; x++) 
	{
	  k = 0;
	  z[0] = 0;
	  z[1] = height;
	  v[0] = 0;
      
	  int v_k = v[k];
	  int v_k_2 = v[k]*v[k];
	  T im_vk_x = im[v_k][x];
	  T *im_x_cp=im[1]+x;
      
	  for (y = 1; y < height; y++, im_x_cp+=width) 
	    {

	      do {
		/* compute vornoi border */
		float s1 = *im_x_cp;
		s1+=(y*y-v_k_2)*scaling_y;
		s1 -= im_vk_x;            // OPTME 
		//	float s2 = (2 * (y - v_k) * scaling_y);
	    
		float s=s1*lut1[y-v_k+height]; //*sy_inv;
		sp = ceil(s);               // OPTME get rid of ceil
	    
		/* case one */
		if (sp >= height)
		  break;
	    
		/* case two */
		if (sp > z[k]) {
		  z[k+1] = int(sp);
		  z[k+2] = height;
		  v[k+1] = y;
		  k++;
		  v_k = v[k];
		  v_k_2 = v_k*v_k;
		  im_vk_x = im[v_k][x];
		  break;
		}
	    
		/* case three */
		if (k == 0) {
		  v[0] = y;
		  v_k = v[k];
		  v_k_2 = v_k*v_k;
		  im_vk_x = im[v_k][x];
		  break;
		} else {
		  z[k] = height;
		  k--;
		  v_k = v[k];
		  v_k_2 = v_k*v_k;
		  im_vk_x = im[v_k][x];
		}
	      } while (1);
	  
	    }
      
	  for (y = 0; y <= k; y++) 
	    vref[y] = im[v[y]][x];

	  k = 0;

	  v_k = v[k];

	  if(need_locations)
	    {
	      int *closest2_row_x = closest2_row[0]+x, *closest2_col_x = closest2_col[0]+x;
	      int *closest_in_row_x = row_locations[v_k]+x, *closest_in_col_x = col_locations[v_k]+x;
	
	      for (y = 0; y < height; y++, closest2_row_x+=width, closest2_col_x+=width) 
	  
		{
		  if (y == z[k+1])
		    {
		      k++;
		      closest_in_row_x += width*(v[k]-v_k);
		      closest_in_col_x += width*(v[k]-v_k);
		      v_k = v[k];
		    }
	    
		  im[y][x] = vref[k] + (v_k-y)*(v_k-y)*scaling_y;
	    
		  *closest2_row_x = *closest_in_row_x;
		  *closest2_col_x = *closest_in_col_x;
		}
	
	    }
	  else
	    {
	      T *im_x_ptr = im[0]+x;

	      for (y = 0; y < height; y++, im_x_ptr += width)
		{
		  if (y == z[k+1])
		    {
		      k++;
		      v_k = v[k];
		    }
	    
		  *im_x_ptr = vref[k] + (v_k-y)*(v_k-y)*scaling_y;   // OPTME ::: get rid of im[y]
		}
	    }
	}

      row_locations = closest2_row, col_locations = closest2_col;
  
      free(z);
      free(v);
      free(vref);

      delete[] lut1;

      return im;
    }
  else
    {
      assert(!need_locations);

      int lut_pivot = max(in_im.rows(), in_im.cols())+5;
      int lut_size = lut_pivot * 2;

      //      static T *lut=0;
      //      if(!lut)
	T *lut = new T[lut_size];
  
      for(int i=0; i<lut_size; i++)
	lut[i] = square(T(i)-lut_pivot) * scaling_x;
  
      _DPlane<T> result(in_im.cols(), in_im.rows());
  
      for(int i=0; i<in_im.rows(); i++)
	dt(in_im[i], result[0], 0, in_im.cols()-1,0, in_im.cols()-1, scaling_x, lut+lut_pivot, in_im.rows(), i);

  
      for(int i=0; i<lut_size; i++)
	lut[i] = square(T(i)-lut_pivot) * scaling_y;

      //  result=result.transpose();
      _DMatrix<T> result2(in_im.rows(), in_im.cols());

      for(int i=0; i<in_im.cols(); i++)
	dt(result[i], result2[0], 0, in_im.rows()-1,0, in_im.rows()-1, scaling_y, lut + lut_pivot, in_im.cols(), i);

      delete[] lut;
      return result2;//.transpose();

    }
}


// distance transform, with arbitrary (i.e. possibly non-diagonal) covariance matrix
template<class T>
_DPlane<T> DistanceTransform_2D<T>::do_transform(const _DPlane<T> & orig_im, bool need_locations)
{
  _DPlane<T> in_im;

  if(rotate)
    in_im = orig_im.rotate_image(-angle);
  else
    in_im = orig_im;


  _DPlane<T> dt_result = do_transform(in_im, need_locations, scaling_x, scaling_y);
 
  if(rotate)
    {
      _DPlane<T> I4 = dt_result.rotate_image(angle);

      DPoint p1 = orig_im.size();
      DPoint p2 = I4.size();
      DPoint new_half = in_im.size() / 2;
      DPoint pp1 = (p2 - p1) / 2;
      DPoint pp2 = (in_im.size() - p1) / 2;

      dt_result = I4.extract(DRect(pp1, pp1+p1-DPoint(1,1)));

      // fix closest_rows, closest_cols here (to compensate for rotation)
      if(need_locations)
	{
	  row_locations = row_locations.rotate_image(angle);
	  col_locations = col_locations.rotate_image(angle);
	  row_locations = row_locations.extract(DRect(pp1, pp1+p1-DPoint(1,1)));
	  col_locations = col_locations.extract(DRect(pp1, pp1+p1-DPoint(1,1)));

	  double cos_ = cos(-angle);
	  double sin_ = sin(-angle);

	  for(int i=0; i<p1.row(); i++)
	    {
	      int* closest_row = row_locations[i], *closest_col = col_locations[i];
	      
	      int c1 = p1.col();
	      for(int j=0; j<c1; j++)
		{
		  DPoint old = DPoint(closest_row[j], closest_col[j]) - new_half;
		  
		  DPoint new_pt(int(cos_ * old.row() + sin_ * old.col()), int(-sin_ * old.row() + cos_ * old.col()));

		  DPoint pt = new_pt + new_half - pp2;
		  closest_row[j] = pt.row(), closest_col[j] = pt.col();
		}
	    }
	}
    }

  return dt_result;
}


#define DECLARE(x)			\
  template class DistanceTransform_2D<x>; 

DECLARE(double);
DECLARE(float);
