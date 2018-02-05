/*#include <FrequencyTransforms.h>

template<class T>
DiscreteCosineTransform<T>::DiscreteCosineTransform(const DPoint &size, bool _normalize) : normalize(_normalize)
{    
  input = _DMatrix<double>(size.row(), size.col());
  output = _DMatrix<double>(size.row(), size.col());
  
  p = fftw_plan_r2r_2d(size.col(), size.row(), input[0], output[0], FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
}

template<class T>
_DMatrix<T> DiscreteCosineTransform<T>::do_transform(const _DMatrix<T> &_input)
{
  change_type(_input, input);
  
  fftw_execute(p);
  
  if(normalize)
    output /= 2.0 * output.rows() * 2.0 * output.cols();
  
  _DMatrix<T> result;
  change_type(output, result);
  
  return result;
}

template<class T>
DiscreteCosineTransform<T>::~DiscreteCosineTransform()
{
  fftw_destroy_plan(p);
}

template<class T>
DiscreteRealFourierTransform<T>::DiscreteRealFourierTransform(const DPoint &_size, bool _normalize)  : size(_size), normalize(_normalize)
{    
  input = _DMatrix<double>(size.row(), size.col());
  output = _DComplexMatrix<double>(size.row(), size.col() / 2 + 1);
  
  p = fftw_plan_dft_r2c_2d(size.row(), size.col(), input[0], reinterpret_cast<fftw_complex*>(output[0]), FFTW_ESTIMATE);
}

template<class T>
_DComplexMatrix<T> DiscreteRealFourierTransform<T>::do_transform(const _DMatrix<T> &_input)
{
  change_type(_input, input);
  
  fftw_execute(p);
  
  _DComplexMatrix<T> result;
  change_type(output, result);
  
  return result;
}

template<class T>
DiscreteRealFourierTransform<T>::~DiscreteRealFourierTransform()
{
  fftw_destroy_plan(p);
}


template<class T>
DiscreteInverseRealFourierTransform<T>::DiscreteInverseRealFourierTransform(const DPoint &_size, bool _normalize)   : size(_size), normalize(_normalize)
{    
  input = _DComplexMatrix<double>(size.row(), size.col() / 2 + 1);
  output = _DMatrix<double>(size.row(), size.col());
  
  p = fftw_plan_dft_c2r_2d(size.row(), size.col(),   reinterpret_cast<fftw_complex*>(input[0]), output[0], FFTW_ESTIMATE);
}

template<class T>
_DMatrix<T> DiscreteInverseRealFourierTransform<T>::do_transform(const _DComplexMatrix<T> &_input)
{
  change_type(_input, input);
  
  fftw_execute(p);
  
  if(normalize)
    output /= output.rows() * output.cols();
  
  _DMatrix<T> result;
  change_type(output, result);
  
  return result;
}

template<class T>
DiscreteInverseRealFourierTransform<T>::~DiscreteInverseRealFourierTransform()
{
  fftw_destroy_plan(p);
}

template<class T>
FFT_Convolution<T>::FFT_Convolution(const DPoint &_image_size, const DPoint &_kernel_size) : kernel_size(_kernel_size), image_size(_image_size)
{
  new_size = kernel_size + image_size - DPoint(1,1);
  
  fft = DiscreteRealFourierTransform<T>(new_size, true);
  ifft = DiscreteInverseRealFourierTransform<T>(new_size, true);
}

template<class T>
void FFT_Convolution<T>::compute_image_fft(const _DPlane<T> &in_image)
{
  assert(in_image.size() == image_size);
  
  _DPlane<T> padded_input(new_size.row(), new_size.col());
  padded_input = 0;
  padded_input.set_submatrix(DPoint(0,0), in_image);
  
  _input_fft = fft.do_transform(padded_input);
}

template<class T>
void FFT_Convolution<T>::compute_kernel_fft(const _DPlane<T> &kernel)
{
  assert(kernel.size() == kernel_size);
  
  _DPlane<T> padded_kernel(new_size.row(), new_size.col());
  padded_kernel = 0;
  padded_kernel.set_submatrix(DPoint(0,0), kernel.rotate_180());
  
  _kernel_fft = fft.do_transform(padded_kernel);
}

template<class T>
_DPlane<T> FFT_Convolution<T>::do_transform()
{
  _DComplexMatrix<T> product = pointwise_multiply(_input_fft, _kernel_fft);
  
  _DPlane<T> p = (ifft.do_transform(product)).extract(DRect(kernel_size/2, kernel_size/2 + image_size - DPoint(1,1)));

  return p;
}

template<class T>
_DPlane<T> FFT_Convolution<T>::do_transform(const _DPlane<T> &in_image, const _DPlane<T> &kernel)
{
  compute_image_fft(in_image);
  compute_kernel_fft(kernel);
  
  return do_transform();
}



#define DECLARE(x) \
  template class FFT_Convolution<x>; \
  template class DiscreteInverseRealFourierTransform<x>; \
  template class DiscreteRealFourierTransform<x>; \
  template class DiscreteCosineTransform<x>;

DECLARE(double);
DECLARE(float);

*/
