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
#include<DMultiDMatrix.h>

using namespace std;

/*
template<class T>
_DMultiDMatrix<T> _DMultiDMatrix<T>::pointwise_min(const _DMultiDMatrix<T> &m1) const
{
  const _DMultiDMatrix &m2 = *this;
  assert(same_size(m1, m2));
  
  _DMultiDMatrix<T> result(m1);
  for(int i=0; i<m1.planes(); i++)
    result.get(i) = _DMatrix<T>::pointwise_min<T>(m1.get(i), m2.get(i));
  
  return result;
}





template<class T>
_DMultiDMatrix<T> _DMultiDMatrix<T>::pointwise_max(const _DMultiDMatrix<T> &m1) const
{
  const _DMultiDMatrix &m2 = *this;
  assert(same_size(m1, m2));
  
  _DMultiDMatrix<T> result(m1);
  for(int i=0; i<m1.planes(); i++)
    result.get(i) = ::pointwise_max(m1.get(i), m2.get(i));
  
  return result;
}

template<class T>
_DMatrixArray<T> _DMatrixArray<T>::pointwise_min(const _DMatrixArray<T> &m1) const
{
  const _DMatrixArray &m2 = *this;
  assert(same_size(m1, m2));
  
  _DMatrixArray<T> result(m1);
  for(int i=0; i<m1.planes(); i++)
    result.get(i) = ::pointwise_min(m1.get(i), m2.get(i));
  
  return result;
}

template<class T>
_DMatrixArray<T> _DMatrixArray<T>::pointwise_max(const _DMatrixArray<T> &m1) const
{
  const _DMatrixArray &m2 = *this;
  assert(same_size(m1, m2));
  
  _DMatrixArray<T> result(m1);
  for(int i=0; i<m1.planes(); i++)
    result.get(i) = ::pointwise_max(m1.get(i), m2.get(i));
  
  return result;
}
*/

template<class T>
_DMatrix<T> pointwise_min(const _DMultiDMatrix<T> &m1, _DMatrix<int> &index) 
{
  index = _DMatrix<int>(m1.rows(), m1.cols());
  index = 0;

  _DMatrix<T> result(m1.get(0));

  int sz = m1.rows() * m1.cols();

  for(int p=1; p<m1.planes(); p++)
    {
      T *res_ptr = result[0];
      int *ind_ptr = index[0];
      T *this_ptr = m1.get(p)[0];

      for(int i=0; i<sz; i++)
	if(this_ptr[i] < res_ptr[i])
	  res_ptr[i] = this_ptr[i], ind_ptr[i] = p;
    }

  return result;
}


template<class T>
_DMatrix<T> pointwise_max(const _DMultiDMatrix<T> &m1, _DMatrix<int> &index)
{
  index = _DMatrix<int>(m1.rows(), m1.cols());
  index = 0;

  _DMatrix<T> result(m1.get(0));

  int sz = m1.rows() * m1.cols();

  for(int p=1; p<m1.planes(); p++)
    {
      if(!same_size(result, m1.get(p)))
	throw string("in pointwise max, not all planes same size");

      T *res_ptr = result[0];
      int *ind_ptr = index[0];
      T *this_ptr = m1.get(p)[0];

      for(int i=0; i<sz; i++)
	if(this_ptr[i] > res_ptr[i])
	  res_ptr[i] = this_ptr[i], ind_ptr[i] = p;
    }

  return result;
}


template<class T2>
bool same_size(const _DMatrixArray<T2> &m1, const _DMatrixArray<T2> &m2)
{
  if(m1.planes() != m2.planes())
    return false;

  for(int i=0; i<m1.planes(); i++)
    if(!same_size(m1.get(i), m2.get(i)))
      return false;

  return true;
}


template<class T2>
bool same_size(const _DMultiDMatrix<T2> &m1, const _DMultiDMatrix<T2> &m2)
{
  return m1.rows() == m2.rows() && m1.cols() == m2.cols() && m1.planes() == m2.planes();
}



#define DECLARE(x) \
  template bool same_size<x>(const _DMultiDMatrix<x> &m1, const _DMultiDMatrix<x> &m2); \
  template bool same_size<x>(const _DMatrixArray<x> &m1, const _DMatrixArray<x> &m2); \
  template _DMatrix<x> pointwise_min(const _DMultiDMatrix<x> &m1, _DMatrix<int> &index); \
  template _DMatrix<x> pointwise_max(const _DMultiDMatrix<x> &m1, _DMatrix<int> &index);


/*  template _DMultiDMatrix<x> pointwise_min<x>(const _DMultiDMatrix<x> &m1, const _DMultiDMatrix<x> &m2); \
    template _DMultiDMatrix<x> pointwise_max<x>(const _DMultiDMatrix<x> &m1, const _DMultiDMatrix<x> &m2); \ */
  /*  template _DMatrixArray<x> pointwise_min<x>(const _DMatrixArray<x> &m1, const _DMatrixArray<x> &m2); \
      template _DMatrixArray<x> pointwise_max<x>(const _DMatrixArray<x> &m1, const _DMatrixArray<x> &m2); \ */



DECLARE(double)
DECLARE(short)
DECLARE(int)
DECLARE(float)
DECLARE(char)
DECLARE(unsigned char)


template class _DMatrixArray<double>;
template class _DMatrixArray<float>;
template class _DMatrixArray<int>;
template class _DMatrixArray<short>;
template class _DMatrixArray<char>;

template class _DMultiDMatrix<double>;
template class _DMultiDMatrix<float>;
template class _DMultiDMatrix<int>;
template class _DMultiDMatrix<short>;
template class _DMultiDMatrix<char>;








