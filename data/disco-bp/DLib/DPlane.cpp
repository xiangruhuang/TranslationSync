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
#include "DPlane.h"
#include "thin_lookup.h"
#include "DGaussianKernel.h"
//
#include <FrequencyTransforms.h>


using namespace std;

template<class T>
void _DPlane<T>::draw(const DRect &rect, T color)
{
  for(int i=rect.top(); i <= rect.bottom(); i++)
    {
      if(i >= 0 && i < _DMatrix<T>::rows() && rect.left() >= 0 && rect.left() < _DMatrix<T>::cols())
	 (*this)[i][rect.left()] = color;
      if(i >= 0 && i < _DMatrix<T>::rows() && rect.right() >= 0 && rect.right() < _DMatrix<T>::cols())
	(*this)[i][rect.right()] = color;
    }

  for(int j=rect.left(); j <= rect.right(); j++)
    {
      if(rect.top() >= 0 && rect.top() < _DMatrix<T>::rows() && j >= 0 && j < _DMatrix<T>::cols())
	(*this)[rect.top()][j] = color;
      if(rect.bottom() >= 0 && rect.bottom() < _DMatrix<T>::rows() && j >= 0 && j < _DMatrix<T>::cols())
	(*this)[rect.bottom()][j] = color;
    }
}


template<class T>
void _DPlane<T>::draw(const DPoint &p1, const DPoint &p2, T color)
{
  float xi=1, yi=1;
  int steps;

  int xdiff = p2.col() - p1.col(), ydiff = p2.row() - p1.row();

  if(abs(xdiff)>abs(ydiff))
    {
      yi=ydiff/(float)abs(xdiff);
      steps=abs(xdiff);
      if(xdiff<0)
	xi=-1;
    }
  else
    {
      xi=xdiff/(float)abs(ydiff);
      steps=abs(ydiff);
      if(ydiff<0)
	yi=-1;
    }

  float cur_x = static_cast< float >(p1.col()), cur_y = static_cast< float >(p1.row());
  for(int i=0; i<steps; i++, cur_x += xi, cur_y += yi)
    {
      DPoint p(int(rint(cur_y)), int(rint(cur_x)));

      if(_DMatrix<T>::is_valid_index(p))
	(*this)[p.row()][p.col()] = color;
    }
}



template<class T>
_DPlane<T> _DPlane<T>::get_x_gradient() const
{
  _DPlane<T> result(_DMatrix<T>::rows()-1, _DMatrix<T>::cols()-1);

  result=0;

  for(int i=0; i<_DMatrix<T>::rows()-1; i++)
    for(int j=0; j<_DMatrix<T>::cols()-1; j++)
      result[i][j] = (*this)[i+1][j+1] - (*this)[i][j];

  return result;
}

template<class T>
_DPlane<T> _DPlane<T>::get_y_gradient() const
{
  _DPlane<T> result(_DMatrix<T>::rows()-1, _DMatrix<T>::cols()-1);

  result=0;

  for(int i=0; i<_DMatrix<T>::rows()-1; i++)
    for(int j=0; j<_DMatrix<T>::cols()-1; j++)
      result[i][j] = (*this)[i+1][j] - (*this)[i][j+1];

  return result;
}

/*
// FIXME:
// Warning : cross_correlate and cross_correlate_fft give different results for
// elements near matrix borders!
#define CROSS_CORR(T) \
template<>  \
_DPlane<T> _DPlane<T>::cross_correlate_fft(const _DPlane<T> &kernel, bool normalize_flag) const \
{ \
  \
  int new_rows = kernel.rows() + rows() - 1; \
  int new_cols = kernel.cols() + cols() - 1; \
\
  _DPlane<T> padded_kernel(new_rows, new_cols), padded_input(new_rows, new_cols);\
\
  DiscreteRealFourierTransform<T> fft(DPoint(new_rows, new_cols), true);\
  DiscreteInverseRealFourierTransform<T> ifft(DPoint(new_rows, new_cols), true);\
\
  padded_kernel = 0;\
  padded_kernel.set_submatrix(DPoint(0,0), kernel.rotate_180());\
\
  padded_input = 0;\
  padded_input.set_submatrix(DPoint(0,0), (*this));\
\
  _DComplexMatrix<T> c1 = fft.do_transform(padded_input);\
  _DComplexMatrix<T> c2 = fft.do_transform(padded_kernel);\
\
  _DComplexMatrix<T> product = pointwise_multiply(c1, c2);\
\
  return (ifft.do_transform(product)).extract(DRect(kernel.size()/2, kernel.size()/2 + _DMatrix<T>::size() - DPoint(1,1)));\
}

CROSS_CORR(float)
CROSS_CORR(double)
*/


template<class T>  
_DPlane<T> _DPlane<T>::cross_correlate_fft(const _DPlane<T> &kernel, bool normalize_flag) const
{
    throw std::string("cross_correlate_fft only implemented for double and float types");
}


template<class T>
_DPlane<T> _DPlane<T>::cross_correlate_1d_rows(const _DPlane<T> &kernel, bool normalize) const
{
  _DPlane<T> result(_DMatrix<T>::rows(), _DMatrix<T>::cols());

  assert(_DMatrix<T>::cols() > kernel.cols());
  assert(kernel.rows() == 1);
  assert(!normalize);

  T *result_cp = result[0];
  T *kernel_cp = kernel[0];
  int half_width = kernel.cols() / 2;

  for(int i=0; i<_DMatrix<T>::rows(); i++)
    {
      T *img_cp = (*this)[i];

      for(int j=0; j < half_width; j++)
	{
	  T sum = 0;
	  for(int l=-half_width; l<=half_width; l++)
	    {
	      T kernel_val = kernel_cp[half_width+l];
	      
	      if(j+l < 0)
		sum += kernel_val * img_cp[0];
	      else
		sum += kernel_val * img_cp[j+l];
	    }

	  *(result_cp++) = sum;
	}
      
      for(int j=half_width; j<_DMatrix<T>::cols() - half_width; j++)
	{
	  T sum = 0;
	  //	  for(int l=-half_width; l<=half_width; l++)
	  T *img_cp2=img_cp+j-half_width;
	  for(int l=0; l<kernel.cols(); l++)
	    {
	      sum += kernel_cp[l] * img_cp2[l];
	    }

	  *(result_cp++) = sum;
	}

      for(int j=_DMatrix<T>::cols() - half_width; j < _DMatrix<T>::cols(); j++)
	{
	  T sum = 0;
	  for(int l=-half_width; l<=half_width; l++)
	    {
	      T kernel_val = kernel_cp[half_width+l];
	      
	      if(j+l >= _DMatrix<T>::cols()) 
		sum += kernel_val * img_cp[_DMatrix<T>::cols()-1];
	      else
		sum += kernel_val * img_cp[j+l];
	    }

	  *(result_cp++) = sum;
	}
    }

  return result;
}

template<class T>
_DPlane<T> _DPlane<T>::cross_correlate_1d_rows_sep(const _DPlane<T> &kernel, bool normalize) const
{
  _DPlane<T> result(_DMatrix<T>::cols(), _DMatrix<T>::rows());

  assert(kernel.rows() == 1);
  assert(!normalize);

  T *kernel_cp = kernel[0];
  int half_width = kernel.cols() / 2;

  int kernel_cols = kernel.cols();

  // optimized version assumes image is wider than kernel
  if(_DMatrix<T>::cols() > kernel.cols())
    {
      //      profiler->begin(10);
      for(int i=0; i<_DMatrix<T>::rows(); i++)
	{
	  T *result_cp = &(result[0][i]);
	  
	  T *img_cp = (*this)[i];

	  //	  profiler->begin(12);
	  for(int j=0; j < half_width; j++)
	    {
	      T sum = 0;
	      for(int l=-half_width; l<=half_width; l++)
		{
		  T kernel_val = kernel_cp[half_width+l];
		  
		  if(j+l < 0)
		    sum += kernel_val * img_cp[0];
		  else
		    sum += kernel_val * img_cp[j+l];
		}
	      
	      *(result_cp) = sum;
	      result_cp += result.cols();
	    }
	  //	  profiler->end(12);
	  
	  //	  profiler->begin(13);
	  int sz = _DMatrix<T>::cols() - half_width;
	  for(int j=half_width; j<sz; j++)
	    {
	      T sum = 0;
	      //	  for(int l=-half_width; l<=half_width; l++)
	      T *img_cp2=img_cp+j-half_width;
	      for(int l=0; l<kernel_cols; l++)
		{
		  sum += kernel_cp[l] * img_cp2[l];
		}
	      
	      *(result_cp) = sum;
	      result_cp += result.cols();
	    }
	  //	  profiler->end(13);
	  
	  //	  profiler->begin(14);
	  sz = _DMatrix<T>::cols();
	  for(int j=_DMatrix<T>::cols() - half_width; j < sz; j++)
	    {
	      T sum = 0;
	      for(int l=-half_width; l<=half_width; l++)
		{
		  T kernel_val = kernel_cp[half_width+l];
		  
		  if(j+l >= _DMatrix<T>::cols()) 
		    sum += kernel_val * img_cp[_DMatrix<T>::cols()-1];
		  else
		    sum += kernel_val * img_cp[j+l];
		}
	      
	      *(result_cp) = sum;
	      result_cp += result.cols();
	    }
	  //	  profiler->end(14);
	}
      //      profiler->end(10);
    }
  else // unoptimized version
    {
      //      profiler->begin(11);
      for(int i=0; i<_DMatrix<T>::rows(); i++)
	{
	  T *result_cp = &(result[0][i]);
	  T *img_cp = (*this)[i];

	  for(int j=0; j < _DMatrix<T>::cols(); j++)
	    {
	      T sum = 0;
	      
	      for(int l=-half_width; l<=half_width; l++)
		{
		  T kernel_val = kernel_cp[half_width+l];
		  
		  if(j+l >= _DMatrix<T>::cols()) 
		    sum += kernel_val * img_cp[_DMatrix<T>::cols()-1];
		  else if(j+l < 0)
		    sum += kernel_val * img_cp[0];
		  else
		    sum += kernel_val * img_cp[j+l];
		}
	      
	      *(result_cp) = sum;
	      result_cp += result.cols();
	    }
	}
      //      profiler->end(11);
    }


  return result;
}

template<class T>
_DPlane<T> _DPlane<T>::flip_horizontally() const
{
  _DPlane<T> result(_DMatrix<T>::rows(), _DMatrix<T>::cols());

  for(int i=0; i<result.rows(); i++)
    {
      const T *in_ptr = (*this)[i];
      T *out_ptr = result[i] + result.cols() -1;
      
      for(int j=0; j<result.cols(); j++)
	*(out_ptr--) = *(in_ptr++);
    }

  return result;
}

template<class T>
_DPlane<T> _DPlane<T>::cross_correlate_1d_cols(const _DPlane<T> &kernel, bool normalize) const
{
  _DPlane<T> result(_DMatrix<T>::rows(), _DMatrix<T>::cols());

  assert(_DPlane<T>::rows() > kernel.rows());
  assert(kernel.cols() == 1);
  assert(!normalize);

  T *kernel_cp = kernel[0];
  int half_width = kernel.rows() / 2;

  for(int j=0; j<_DMatrix<T>::cols(); j++)
    {
      T *img_cp = &((*this)[0][j]);
      T *result_cp = &(result[0][j]);

      for(int i=0; i < half_width; i++)
	{
	  T sum = 0;
	  for(int l=-half_width; l<=half_width; l++)
	    {
	      T kernel_val = kernel_cp[half_width+l];
	      
	      if(i+l < 0)
		sum += kernel_val * img_cp[0];
	      else
		sum += kernel_val * img_cp[(i+l)*_DMatrix<T>::cols()];
	    }

	  *result_cp = sum;
	  result_cp += _DMatrix<T>::cols();
	}
      
      for(int i=half_width; i<_DMatrix<T>::rows() - half_width; i++)
	{
	  T sum = 0;
	  T *img_cp_cp = img_cp + _DMatrix<T>::cols() * (i-half_width);
	  for(int l=0; l<kernel.rows(); l++)
	    {
	      sum += kernel_cp[l] * (*img_cp_cp);
	      img_cp_cp += _DMatrix<T>::cols();
	    }

	  *result_cp = sum;
	  result_cp += _DMatrix<T>::cols();
	}

      for(int i=_DMatrix<T>::rows() - half_width; i < _DMatrix<T>::rows(); i++)
	{
	  T sum = 0;
	  for(int l=-half_width; l<=half_width; l++)
	    {
	      T kernel_val = kernel_cp[half_width+l];
	      
	      if(i+l >= _DMatrix<T>::rows()) 
		sum += kernel_val * img_cp[(_DMatrix<T>::rows()-1)*_DMatrix<T>::cols()];
	      else
		sum += kernel_val * img_cp[(i+l)*_DMatrix<T>::cols()];
	    }

	  *result_cp = sum;
	  result_cp += _DMatrix<T>::cols();
	}
    }

  return result;
}

template<class T>
_DPlane<T> _DPlane<T>::cross_correlate_separable(const _DPlane<T> &row_kernel, const _DPlane<T> &col_kernel, bool normalize) const
{
  return (cross_correlate_1d_rows_sep(row_kernel, normalize)).cross_correlate_1d_rows_sep(col_kernel.transpose(), normalize);
}


// cross_correlate with a gaussian
template<class T>
_DPlane<T> _DPlane<T>::convolve_gaussian(const _DMatrix<T> &sigma, float sigma_width, bool box_approximate) const
{
  bool rotate = true;
  double angle=0;
  _DPlane<T> in_im(*this);
  T row_sigma = sigma[1][1], col_sigma = sigma[0][0];

  if(sigma[0][1] == 0 && sigma[1][0] == 0)
    rotate=false;

  if(rotate)
    {
      _DMatrix<T> eig_vector;
      _DMatrix<T> res = sigma.eigen(eig_vector);
      
      if(eig_vector[0][1] ==0)
	angle=0;
      else
	angle = -atan((double) eig_vector[1][1]/eig_vector[0][1]);
      
      in_im = in_im.rotate_image(angle, 0);

      col_sigma = res[0][1];
      row_sigma = res[0][0];
    }

  int f_width = int(ceil(sqrt((double) row_sigma)*sigma_width))*2+1;
  int f_height = int(ceil(sqrt((double) col_sigma)*sigma_width))*2+1;
  if(f_width >= in_im.cols()-1)
    f_width = in_im.cols()-1 - (in_im.cols()) % 2;
  if(f_height >= in_im.rows()-1)
    f_height = in_im.rows()-1 - (in_im.rows()) % 2;

  //  profiler->begin(9); 
  //  _DPlane<T> result;
  //  if(!box_approximate)
  _DPlane<T>  result = in_im.cross_correlate_separable(_DGaussianKernel<T>(sqrt((double) row_sigma), f_width),
					     _DGaussianKernel<T>(sqrt((double) col_sigma), f_height).transpose(), false);
  //  else
  //    result = in_im.convolve_gaussian_approximate(sqrt(row_sigma), sqrt(col_sigma));

  //  profiler->end(9);
  if(rotate)
    {
      _DPlane<T> I4 = result.rotate_image(-angle, 0);
      
      int r1 = this->rows(), c1 = this->cols();
      int r2 = I4.rows(), c2 = I4.cols();
      
      int rr1 = ((r2-r1)/2);
      int cc1 = ((c2-c1)/2);
      
      result = I4.extract(DRect(rr1, cc1, rr1+r1-1, cc1+c1-1));
    }

  return result;
}


template<class T>
_DPlane<T> _DPlane<T>::box_filter_rows_1D(int kernel_cols, bool normalize) const
{
  _DPlane<T> result(cols(), rows());

  int half_width = kernel_cols / 2;

  // optimized version assumes image is wider than kernel
  if(_DMatrix<T>::cols() > kernel_cols)
    {
      for(int i=0; i<_DMatrix<T>::rows(); i++)
	{
	  T *result_cp = &(result[0][i]);
	  
	  T *img_cp = (*this)[i];
	  
	  for(int j=0; j < half_width; j++)
	    {
	      T sum = 0;
	      for(int l=-half_width; l<=half_width; l++)
		{
		  if(j+l < 0)
		    sum += img_cp[0];
		  else
		    sum += img_cp[j+l];
		}
	      
	      *(result_cp) = sum;
	      result_cp += result.cols();
	    }
	  
	  T cum_sum=0;
	  for(int l=0; l<kernel_cols; l++)
	    cum_sum += img_cp[l];
	  *(result_cp) = cum_sum;
	  result_cp += result.cols();

	  int sz = _DMatrix<T>::cols() - half_width;
	  for(int j=half_width+1; j<sz; j++)
	    {
	      cum_sum = *(result_cp) = cum_sum - img_cp[j-half_width-1] + img_cp[j-half_width-1+kernel_cols];

	      result_cp += result.cols();
	    }
	  
	  sz = _DMatrix<T>::cols();
	  for(int j=_DMatrix<T>::cols() - half_width; j < sz; j++)
	    {
	      T sum = 0;
	      for(int l=-half_width; l<=half_width; l++)
		{
		  if(j+l >= _DMatrix<T>::cols()) 
		    sum += img_cp[_DMatrix<T>::cols()-1];
		  else
		    sum += img_cp[j+l];
		}
	      
	      *(result_cp) = sum;
	      result_cp += result.cols();
	    }
	}
    }
  else // unoptimized version
    {
      for(int i=0; i<_DMatrix<T>::rows(); i++)
	{
	  T *result_cp = &(result[0][i]);
	  T *img_cp = (*this)[i];

	  for(int j=0; j < _DMatrix<T>::cols(); j++)
	    {
	      T sum = 0;
	      
	      for(int l=-half_width; l<=half_width; l++)
		{
		  if(j+l >= _DMatrix<T>::cols()) 
		    sum += img_cp[_DMatrix<T>::cols()-1];
		  else if(j+l < 0)
		    sum += img_cp[0];
		  else
		    sum += img_cp[j+l];
		}
	      
	      *(result_cp) = sum;
	      result_cp += result.cols();
	    }
	}
    }

 
  return normalize ? result * (1.0 / T(kernel_cols)) : result;
}


template<class T>
_DPlane<T> _DPlane<T>::box_filter(int kernel_rows, int kernel_cols, bool normalize) const
{
  return box_filter_rows_1D(kernel_cols, normalize).box_filter_rows_1D(kernel_rows, normalize);
}


static int factorial(int in)
{
  if(in == 0) return 1;
  else return in*factorial(in-1);
}



template<class T>
_DPlane<T> _DPlane<T>::convolve_gaussian_approximate(T row_sigma, T col_sigma)
{
  // compute # of columns necessary for box filter;
  const int reps = 4;
  float kernel_cols = float( 0.005 ), kernel_rows = float( 0.005 );
  
  for(int i=0; i<=(reps-1)/2; i++)
    kernel_cols += (float) (pow(-1.0,i) * reps / float(factorial(i) * factorial(reps-i)) * pow(float(reps/2 - i), reps-1));
  kernel_cols *= float( row_sigma * sqrt(2 * M_PI) );

  for(int i=0; i<=(reps-1)/2; i++)
    kernel_rows += (float) (pow(-1.0,i) * reps / float(factorial(i) * factorial(reps-i)) * pow(float(reps/2 - i), reps-1));
  kernel_rows *= float( col_sigma * sqrt(2 * M_PI) );

  kernel_rows = ceil(kernel_rows);
  kernel_cols = ceil(kernel_cols);

  std::cout << kernel_cols << " " << kernel_rows << std::endl;

  _DPlane<T> result = (*this);

  for(int i=0; i<reps; i++)
    result = result.box_filter_rows_1D((int)kernel_cols, true);

  result = result.transpose();

  for(int i=0; i<reps; i++)
    result = result.box_filter_rows_1D((int)kernel_rows, true);
  
  return result.transpose();
}



template<class T>
_DPlane<T> _DPlane<T>::cross_correlate(const _DPlane<T> &kernel, bool normalize) const
{
  _DPlane<T> result(_DMatrix<T>::rows(), _DMatrix<T>::cols());

  result = 0;

  int half_height = kernel._DMatrix<T>::rows()/2;
  int half_width = kernel._DMatrix<T>::cols()/2;

  int low_half_height = -half_height, low_half_width = -half_width;
  int high_half_height = half_height, high_half_width = half_width;
  if(!(half_height % 2))
    high_half_height--;
  if(!(half_width % 2))
    high_half_width--;

  T kernel_norm;

  if(normalize)
    kernel_norm = kernel.sum();
  else
    kernel_norm = T(1.0);

  for(int i=0; i<_DMatrix<T>::rows(); i++)
    for(int j=0; j<_DMatrix<T>::cols(); j++)
    {
      T sum=0;
      for(int k=low_half_height; k<=high_half_height; k++)
	{
	  T *kernel_cp, *img_cp;
	  kernel_cp = kernel[k+half_height];

	  if(i+k < 0)
	    img_cp = (*this)[0];
	  else if(i+k >= _DMatrix<T>::rows())
	    img_cp = (*this)[_DMatrix<T>::rows()-1];
	  else
	    img_cp = (*this)[i+k];

	  if(j-half_width < 0 || j + half_width >= _DMatrix<T>::cols())
	    {
	      for(int l=low_half_width; l<=high_half_width; l++)
		{
		  T kernel_val = kernel_cp[half_width+l];

		  if(j+l < 0)
		    sum += kernel_val * img_cp[0];
		  else if(j+l >= _DMatrix<T>::cols()) 
		    sum += kernel_val * img_cp[_DMatrix<T>::cols()-1];
		  else
		    sum += kernel_val * img_cp[j+l];
		}
	    }
	  else
	    {
	      for(int l=low_half_width; l<=high_half_width; l++)
		{
		  sum += kernel_cp[l+half_width] * img_cp[j+l];
		}
	    }
	}

      result[i][j] = sum/kernel_norm;
    }


  return result;
}

template<class T>
_DPlane<T> _DPlane<T>::subsample(int row_factor, int col_factor)
{
  int r = (int) ceil((double) (_DMatrix<T>::rows()) / row_factor);
  int c = (int) ceil((double) (_DMatrix<T>::cols()) / col_factor);
  _DPlane<T> result(r, c);

  T *out_cp = result[0];
  int n=0;
  for(int i=0; i<_DMatrix<T>::rows(); i+=row_factor)
    for(int j=0; j<_DMatrix<T>::cols(); j+=col_factor, out_cp++,n++)
      *out_cp = (*this)[i][j];

  return result;
}

// sub-samples, taking the *maximum* value within each factor x factor block
// 
// currently only works for factor == 2
//
// just chops off odd row or column
template<class T>
_DPlane<T> _DPlane<T>::downsample_max(int factor) const
{
  //  assert(factor == 2 || factor == 3);

  int r = (int) floor(float(_DMatrix<T>::rows()) / factor);
  int c = (int) floor(float(_DMatrix<T>::cols()) / factor);
  _DPlane<T> result(r, c);
  
  if(factor == 2)
    {
      T *out_cp = result[0];
      for(int i=0; i < r; i++)
	{
	  const T *in1_cp = (*this)[i*2];
	  const T *in2_cp = (*this)[i*2 + 1];
	  
	  for(int j=0; j < c; j++, in1_cp+=2, in2_cp+=2, out_cp++)
	    *out_cp = std::max( std::max(*in1_cp, *(in1_cp+1)), std::max(*in2_cp, *(in2_cp+1)) );
	}
    }
  else if(factor==3) 
    {
      T *out_cp = result[0];
      for(int i=0; i < r; i++)
	{
	  const T *in1_cp = (*this)[i*3];
	  const T *in2_cp = (*this)[i*3 + 1];
	  const T *in3_cp = (*this)[i*3 + 2];
	  
	  for(int j=0; j < c; j++, in1_cp+=3, in2_cp+=3, in3_cp+=3, out_cp++)
	    *out_cp = std::max( std::max( std::max( std::max(*in1_cp, *(in1_cp+1)), std::max(*(in1_cp+2), *(in2_cp))),
					  std::max( std::max(*(in2_cp+1), *(in2_cp+2)), std::max(*in3_cp, *(in3_cp+1)))),
				*(in3_cp+2));
					  
	}
    }
  else 
    {
      int in_rows = int(_DMatrix<T>::rows() / factor)*factor, in_cols = int(_DMatrix<T>::cols() / factor)*factor;
      _DPlane<T> tmp(in_rows, c);

      for(int i=0; i<in_rows; i++)
	{
	  T *out_cp = tmp[i];
	  T *in_ptr = (*this)[i];

	  for(int j=0; j<in_cols; ++out_cp)
	    {
	      T local_max = in_ptr[j];
	      j++;

	      for(int k=1; k<factor; k++, j++)
		{
		  if(in_ptr[j] > local_max)
		    local_max = in_ptr[j];
		}
		
	      *out_cp = local_max;
	    }
	}

      for(int j=0; j<c; j++)
	{
	  T *out_cp = result[0]+j;
	  T *in_ptr = tmp[0]+j;

	  for(int i=0; i<in_rows; out_cp += c)
	    {
	      T local_max = *in_ptr;
	      i++;
	      in_ptr+=c;

	      for(int k=1; k<factor; k++, i++, in_ptr+=c)
		{
		  if(*in_ptr > local_max)
		    local_max = *in_ptr;
		}
		
	      *out_cp = local_max;
	    }
	}
    }

  return result;
}


inline
double float_part(double i)
{
  return i - floor(i);
}

template<class T>
_DPlane<T> _DPlane<T>::bilinear_interpolate(T row_offset, T col_offset) const
{
  assert(row_offset >= 0 && col_offset >= 0);
  assert(row_offset < 1 && col_offset < 1);

  int new_rows = _DMatrix<T>::rows() - (int)ceil((double) row_offset);
  int new_cols = _DMatrix<T>::cols() - (int)ceil((double) col_offset);

  _DPlane<T> result(new_rows, new_cols);

  T row_pos, col_pos;

  row_pos = row_offset;
  for(int i=0; i<new_rows; i++, row_pos++)
    {
      col_pos = col_offset;
      for(int j=0; j<new_cols; j++, col_pos++)
	{
	  double row_fraction = float_part(row_pos);
	  double col_fraction = float_part(col_pos);
	  int row_int = (int)floor((double) row_pos);
	  int col_int = (int)floor((double) col_pos);
	  
	  result[i][j] = 
	    T((1.0-row_fraction) * (1.0-col_fraction) * (*this)[row_int][col_int] +
	    (1.0-row_fraction) * (    col_fraction) * (*this)[row_int][col_int+1] +
	    (    row_fraction) * (1.0-col_fraction) * (*this)[row_int+1][col_int] +
	      (    row_fraction) * (    col_fraction) * (*this)[row_int+1][col_int+1]);
	    
	}
    }

  return result;
}


template<class T>
void _DPlane<T>::get_subsample_vector(float start_row, float row_count, int &sz, float *&ptr) const
{
  if(floor(start_row + row_count) == floor(start_row)) // start and end in same pixel
    {
      sz = 1;
      ptr = new float[sz];
      ptr[0] = row_count;
    }
  else 
    {
      int i_part = int(start_row);
      double f_part = float_part(start_row);

      int ii=0;
      sz = int(ceil(start_row + row_count) - floor(start_row));   
      ptr = new float[sz];

      if(f_part > 0)
	ptr[ii++] = static_cast< float >(1.0 - f_part); // first pixel

      for(int i=0; i<int(floor(start_row + row_count) - ceil(start_row)); i++) 
	ptr[ii++] = 1.0;

      i_part = int(start_row + row_count);
      f_part = float_part(start_row + row_count);

      if(f_part > 0 && ii < sz)
	ptr[ii++] = static_cast< float >(f_part); // last pixel

      if(ii != sz)
	{
	  cout << start_row << " " << row_count << " " << sz << endl;
	  cout << "--- " << i_part << " " << f_part << " " << sz << endl;
	  assert(0);
	}
    }
}


template<class T>
_DPlane<T> _DPlane<T>::bilinear_rescale(float row_factor, float col_factor) const
{
  _DPlane<T> result1(int(this->rows() * row_factor), this->cols());
  _DPlane<T> result2(int(this->rows() * row_factor), int(this->cols() * col_factor));
  float row_inc = float( 1.0 ) / row_factor;

  assert(!_DMatrix<T>::is_any_nan());

  if(row_factor == 1.0)
    result1 = *this;
  else if(row_factor < 1.0)  // subsample
    {
      for(int i=0; i<result1.rows(); i++)
	{
	  int sz;
	  float *vec;
	  get_subsample_vector(i * row_inc, row_inc, sz, vec);
	  int i_beg = int(i * row_inc);
	  T *result1_cp = result1[i];
	  
	  for(int j=0; j<result1.cols(); j++)
	    {
	      float sum = 0;
	      float norm=0;

	      int i_max = i_beg + sz;
	      for(int ii = i_beg, cp=0; ii < i_max; ii++, cp++)
		if(ii < this->rows())
		  sum += (*this)[ii][j] * vec[cp], norm += vec[cp];
	      
	      result1_cp[j] = T(sum / norm);
	    }

	  delete[] vec;
	}
    }
  else
    {
      for(int i=0; i<result1.rows(); i++)
	{
	  double new_i = i * row_inc;
	  int i_part = int(new_i);
	  double f_part = float_part(new_i);
	  T *in_ptr_i = (*this)[i_part], *in_ptr_i1 = (*this)[i_part+1], *result1_ptr = result1[i];

	  for(int j=0; j<result1.cols(); j++)
	    {
	      if(i_part+1 < this->rows())
		result1_ptr[j] = T(in_ptr_i[j] * (1.0 - f_part) + in_ptr_i1[j] * f_part);
	      else
		result1_ptr[j] = T(in_ptr_i[j]);
	    }
	}
    }

  assert(!result1.is_any_nan());
  float col_inc = float( 1.0 ) / col_factor;

  if(col_factor == 1.0)
    result2 = result1;
  else if(col_factor < 1.0)  // subsample
    {
      for(int j=0; j<result2.cols(); j++)
	{
	  float *vec = NULL;
	  int sz = 0;
	  get_subsample_vector(j * col_inc, col_inc, sz, vec);
	  int j_beg = int(j * col_inc);
	  
	  for(int i=0; i<result1.rows(); i++)
	    {
	      float sum = 0;
	      float norm = 0;

	      T *result1_i = result1[i];

	      for(int jj = j_beg, cp=0; jj < j_beg + sz; jj++, cp++)
		if(jj < this->cols())
		  {
		    sum += result1_i[jj] * vec[cp], norm += vec[cp];
		  }
	      
	      result2[i][j] = T(sum / norm);
	    }

	  delete[] vec;
	}
    }
  else
    {
      for(int i=0; i<result2.rows(); i++)
	{
	  T *result2_ptr = result2[i], *result1_ptr = result1[i];

	  for(int j=0; j<result2.cols(); j++)
	    {
	      double new_j = j * col_inc;
	      int i_part = int(new_j);
	      double f_part = float_part(new_j);
	      
	      if(i_part+1 < this->cols())
		result2_ptr[j] = T(result1_ptr[i_part] * (1.0 - f_part) + result1_ptr[i_part+1] * f_part);
	      else
		result2_ptr[j] = T(result1_ptr[i_part]);
	    }
	}
    }
  assert(!result2.is_any_nan());

  return result2;
}


template<class T>
_DPlane<T> _DPlane<T>::rotate_image(double angle, T bg_val) const
{
  //  cerr << "in rotate image" << endl;
  DRect current_size = DRect(DPoint(0,0), _DMatrix<T>::size());
  DPoint new_corner1 = current_size.top_left().rotate(static_cast< float >(angle));
  DPoint new_corner2 = current_size.top_right().rotate(static_cast< float >(angle));
  DPoint new_corner3 = current_size.bottom_left().rotate(static_cast< float >(angle));
  DPoint new_corner4 = current_size.bottom_right().rotate(static_cast< float >(angle));

  //  cout << " -- corners " << new_corner1 << " " << new_corner2 << " " << new_corner3 << " " << new_corner4 << endl;

  DPoint new_bottomright = elementwise_max( elementwise_max( new_corner1, new_corner2 ), elementwise_max( new_corner3, new_corner4 ) );
  DPoint new_topleft = elementwise_min( elementwise_min( new_corner1, new_corner2 ), elementwise_min( new_corner3, new_corner4 ) );
  DPoint new_size = new_bottomright - new_topleft + DPoint(5,5);

  //  cout << "  - sizes " << new_bottomright << " " << new_size << " " << new_topleft << endl;

  int new_rows = new_size.row() + 1-(new_size.row() % 2);
  int new_cols = new_size.col() + 1-(new_size.col() % 2);

  //  cerr << new_rows << " " << new_cols << endl;
  _DPlane<T> new_I(new_rows, new_cols);

  int row_half = (_DMatrix<T>::rows()/2);
  int col_half = (_DMatrix<T>::cols()/2);

  int new_row_half = (new_rows/2);
  int new_col_half = (new_cols/2);

  float c = static_cast< float >(cos(angle));
  float s = static_cast< float >(sin(angle));

  T *cp = new_I[0];
  for(int i=0; i<new_rows; ++i)
    {
      float i2 = (c*(i-new_row_half)+s*(-1-new_col_half)) + row_half;
      float j2 = (-s*(i-new_row_half)+c*(-1-new_col_half)) + col_half;

      for(int j=0; j<new_cols; ++j, ++cp)
	{
      	  i2 += s; j2 += c;

	  if(i2 >= 0 && j2 >= 0 && i2 < _DMatrix<T>::rows()-1 && j2 < _DMatrix<T>::cols()-1)
	    {
	      // compute upper-left pixel
	      int low_i = int(i2);
	      int low_j = int(j2);
	      float i_diff = i2-low_i;
	      float j_diff = j2-low_j;

	      T *cp1 = (*this)[low_i];
	      T *cp2 = (*this)[low_i+1];
	      
	      *cp = T(((1-i_diff)*(1-j_diff)*cp1[low_j]) +
		      ((1-i_diff)*j_diff*cp1[low_j+1]) + 
		      (i_diff*(1-j_diff)*cp2[low_j]) +
		      (i_diff * j_diff * cp2[low_j+1]));
	    }
	  else
	    {
	      *cp = T(bg_val);
	    }
	}
    }
  
  //  cerr << "returning new I" << endl;
  return new_I;
}


template<class T>
_DPlane<T> _DPlane<T>::rotate_image_nn(double angle) const
{
  int new_rows = int(ceil(sqrt((double) _DMatrix<T>::rows()*_DMatrix<T>::rows()+_DMatrix<T>::cols()*_DMatrix<T>::cols())));
  int new_cols = new_rows;

  _DPlane<T> new_I(new_rows, new_cols);
  new_I = T(1e100);

  int row_half = (_DMatrix<T>::rows()/2);
  int col_half = (_DMatrix<T>::cols()/2);

  int new_row_half = (new_rows/2);
  int new_col_half = (new_cols/2);

  float c = static_cast< float >(cos(angle));
  float s = static_cast< float >(sin(angle));


  T *cp = new_I[0];
  for(int i=0; i<new_rows; ++i)
    {
      float i2 = (c*(i-new_row_half)+s*(-1-new_col_half)) + row_half;
      float j2 = (-s*(i-new_row_half)+c*(-1-new_col_half)) + col_half;

      for(int j=0; j<new_cols; ++j, ++cp)
	{
      	  i2 += s; j2 += c;

	  if(i2 < 0 || j2 < 0 || i2>=_DMatrix<T>::rows()-1 || j2>=_DMatrix<T>::cols()-1)
	    continue;

	  int rnd_i = (int) round(i2);
	  int rnd_j = (int) round(j2);

	  *cp = T((*this)[rnd_i][rnd_j]);
	  
	}
    }

  return new_I;
}

template<class T>
_DPlane<T> _DPlane<T>::binary_thin(void) const
{
  static bool lookup_table[256];

  make_thin_lookup_table(lookup_table);

  _DPlane<T> result = (*this);

  for(int i=1; i<_DMatrix<T>::rows()-1; i++)
    for(int j=1; j<_DMatrix<T>::cols()-1; j++)
      {
	int val=0, n=1;

	for(int i2=-1; i2<=1; i2++)
	  {
	    T *result_row_cp = result[i+i2];
	    for(int j2=-1; j2<=1; j2++)
	      {
		if(i2 == 0 && j2 == 0)
		  continue;
		
		if(result_row_cp[j+j2])
		  val += n; 
		
		n=n << 1;
	      }
	  }

	if(!lookup_table[val])
	  result[i][j] = 0;
      }

  return result;
}

template<class T>
_DPlane<T> _DPlane<T>::dilate(int square_size) const
{
  _DPlane<T> kernel(square_size, square_size);
  kernel = 1;

  return dilate(kernel);
}

template<class T>
_DPlane<T> _DPlane<T>::dilate(const _DPlane<T> &kernel) const
{
  return dilate(kernel, DPoint(kernel.rows()/2, kernel.cols()/2));
}

template<class T>
_DPlane<T> _DPlane<T>::dilate(const _DPlane<T> &kernel, const DPoint &se_origin) const
{
  int seorigin = kernel.cols() * se_origin.row() + se_origin.col();
  int *se2 = new int[kernel.total_pixel_count()];
  T *se = kernel[0];

  int numpoints = 0;
  for (int i=0; i < kernel.total_pixel_count(); i++) 
    {
      if(se[i])
	se2[numpoints++] = ( (i / kernel.cols()) - seorigin / kernel.cols()) * this->cols() +
	  ((i % kernel.cols())-(seorigin % kernel.cols()));
    }
  
  _DPlane<T> out(this->rows(), this->cols());
  out = 0;

  T *outbuf = out[0];
  T *sbuf = (*this)[0];
  int i_end = this->total_pixel_count() - kernel.rows() * this->cols() + kernel.cols();
  for (int i = kernel.rows() * this->cols() + kernel.cols(); i < i_end; i++)
    {
      if (sbuf[i]) 
	{
	  outbuf[i] = 1;

	  for (int j=0; j < numpoints; j++) 
	    outbuf[i+se2[j]] = 1;
	}
    }

  delete[] se2;

  return out;
}




#define DECLARE(x) \
  template class _DPlane<x>;

DECLARE(double)
DECLARE(short)
DECLARE(int)
DECLARE(float)
DECLARE(char)
DECLARE(unsigned char)






