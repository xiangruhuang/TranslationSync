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
#include "DMatrix.h"
#include <string>
#include <vector>
#include <complex>
#include <fstream>

using namespace std;

template<class T>
void _DTwoDimArray<T>::deallocate_storage()
{
  if(data)
    {
      delete[] data;
      delete[] data_area;

      data = 0;
      data_area = 0;
    }
}

template<class T>
void _DTwoDimArray<T>::initialize_storage()
{
  //  profiler->begin(2);
  if(data)
    deallocate_storage();

  if(_rows > 0)
    {
      data = new T *[_rows];
      data_area = new T[_rows * _cols];
    }
  else
    {
      data = 0;
      data_area = 0;
    }
  
  T *cp = data_area;
  for(int i=0; i<_rows; i++, cp+=_cols)
    data[i] = cp;
  //  profiler->end(2);
}

template<class T>
_DTwoDimArray<T>::_DTwoDimArray()
{
  data = 0;
  data_area = 0;
  _rows = _cols = 0;
}

template<class T>
_DTwoDimArray<T>::_DTwoDimArray(int __rows, int __cols)
{
    _rows = __rows;
    _cols = __cols;

    data = 0;
    data_area = 0;

    initialize_storage();
}

template<class T>
_DTwoDimArray<T>::_DTwoDimArray(int __rows, int __cols, const T *array)
{
    _rows = __rows;
    _cols = __cols;

    data = 0;
    data_area = 0;

    initialize_storage();

    memcpy(data_area, array, _rows * _cols * sizeof(T));
}

template<class T>
_DMatrix<T>::_DMatrix(const std::vector<T> &vec)
{
  _rows = 1;
  _cols = (int) vec.size();

  data = 0;
  data_area = 0;
  
  initialize_storage();

  set_row(0, vec);
}


template<class T>
_DMatrix<T>::_DMatrix(int __rows, int __cols, T val)
{
  _rows = __rows;
  _cols = __cols;
  
  data = 0;
  data_area = 0;
  
  initialize_storage();
  
  (*this) = val;
}

template<class T>
_DMatrix<T>::_DMatrix(int __rows, int __cols, matrix_init_type type)
{
  _rows = __rows;
  _cols = __cols;
  
  data = 0;
  data_area = 0;
  
  initialize_storage();
  
  if(type == random)
    {
      int sz = _rows * _cols;
      for(int i=0; i<sz; i++)
	data_area[i] = T(drand48());
    }
  else if(type == identity)
    {
      for(int i=0; i<_rows; i++)
	for(int j=0; j<_cols; j++)
	  if(i == j)
	    (*this)[i][j] = 1;
	  else
	    (*this)[i][j] = 0;
    }
  else if(type == zeros)
    {
      (*this) = 0;
    }
  else if(type == ones)
    {
      (*this) = 1;
    }
  else
    throw(std::string("unknown matrix initializer type"));

}

template<class T>
bool same_size(const _DMatrix<T> &m1, const _DMatrix<T> &m2)
{
    return(m1.rows() == m2.rows() && m1.cols() == m2.cols());
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator+(const _DMatrix<T> &other) const
{
    _DMatrix<T> result(_rows, _cols);

    if(!same_size(*this, other))
        throw string("Size mismatch in DMatrix operator +");

    T *cp1 = data_area, *cp2 = other.data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; j++, cp1++, cp2++, cp_out++)
            *cp_out = *cp1 + * cp2;

    return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator+(T value) const
{
    _DMatrix<T> result(_rows, _cols);

    T *cp1 = data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; j++, cp1++, cp_out++)
            *cp_out = *cp1 + value;

    return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator-(const _DMatrix<T> &other) const
{
    _DMatrix<T> result(_rows, _cols);

    if(!same_size(*this, other))
        throw string("Size mismatch in DMatrix operator -");

    T *cp1 = data_area, *cp2 = other.data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; ++j, ++cp1, ++cp2, ++cp_out)
            *cp_out = *cp1 - * cp2;

    return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator-(T value) const
{
    _DMatrix<T> result(_rows, _cols);

    T *cp1 = data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; ++j, ++cp1, ++cp_out)
            *cp_out = *cp1 - value;

    return result;
}

template<class T>
_DMatrix<T> operator-(const _DMatrix<T> &m)
{
  _DMatrix<T> result(m.rows(), m.cols());

  T *cp_in = m.data_ptr(), *cp_out = result.data_ptr();

  for(int i=0; i<m.rows(); ++i)
    for(int j=0; j<m.cols(); ++j, ++cp_in, ++cp_out)
      *cp_out = -*cp_in;

  return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator*(const _DMatrix<T> &other) const
{
    if(_cols != other._rows)
        throw string("Size mismatch in DMatrix operator *");

    _DMatrix<T> result(_rows, other._cols);

    for(int i=0; i<_rows; i++)
      {
	T *data_i = data[i];

        for(int j=0; j<other._cols; j++)
	  {
            T res=0;
	    
            for(int k=0; k<other._rows; k++)
	      res += data_i[k] * other.data[k][j];
	    
            result[i][j] = res;
	  }
      }

    return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator*(T value) const
{
    _DMatrix<T> result(_rows, _cols);

    T *cp1 = data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; j++, cp1++, cp_out++)
            *cp_out = T(*cp1 * value);

    return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator/(double value) const
{
    _DMatrix<T> result(_rows, _cols);

    T *cp1 = data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
      for(int j=0; j<_cols; j++, cp1++, cp_out++)
	*cp_out = T(*cp1 / value);
    
    return result;
}

template<class T>
_DMatrix<T> operator*(T value, const _DMatrix<T> &other)
{
    return(other * value);
}


template<class T>
void _DMatrix<T>::operator+=(const _DMatrix<T> &other)
{
  assert(same_size(other, *this));

  T *in = other.data_ptr();
  T *out = data_ptr();

  int count = rows() * cols();
  for(int i=0; i<count; i++)
    out[i] += in[i];
}

template<class T>
void _DMatrix<T>::operator-=(const _DMatrix<T> &other)
{
  assert(same_size(other, *this));

  T *in = other.data_ptr();
  T *out = data_ptr();

  int count = rows() * cols();
  for(int i=0; i<count; i++)
    out[i] -= in[i];
}

template<class T>
void _DMatrix<T>::operator+=(T value)
{
  T *out = data_ptr();

  int count = rows() * cols();
  for(int i=0; i<count; i++)
    out[i] += value;
}
  
template<class T>
void _DMatrix<T>::operator-=(T value)
{
  T *out = data_ptr();

  int count = rows() * cols();
  for(int i=0; i<count; i++)
    out[i] -= value;
}

template<class T>
void _DMatrix<T>::operator*=(double value)
{
  T *out = data_ptr();

  int count = rows() * cols();
  for(int i=0; i<count; i++)
    out[i] = T(out[i] * value);
}

template<class T>
void _DMatrix<T>::operator/=(T value)
{
  T *out = data_ptr();

  int count = rows() * cols();
  for(int i=0; i<count; i++)
    out[i] /= value;
}
  


template<class T>
T& _DMatrix<T>::operator[](const DPoint &point) const
{
    return data[point.row()][point.col()];
}


template<class T>
_DTwoDimArray<T>::~_DTwoDimArray()
{
    deallocate_storage();
}

template<class T>
_DTwoDimArray<T> &_DTwoDimArray<T>::operator=(const _DTwoDimArray<T> &other)
{
  if(this == &other)
    return *this;

  // profiler->begin(4);
  if(!data || _rows != other.rows() || _cols != other.cols())
    {
      _rows = other.rows();
      _cols = other.cols();
      
      initialize_storage();
    }

  memcpy(data_area, other.data_area, _rows * _cols * sizeof(T));
  
// profiler->end(4);
  return *this;
}


template<class T>
_DMatrix<T> &_DMatrix<T>::operator=(T other)
{
    T *cp = data_ptr();
    int sz = rows() * cols();
    for(int i=0; i<sz; i++)
      *(cp++) = other;

    return *this;
}

template<class T>
_DTwoDimArray<T>::_DTwoDimArray(const _DTwoDimArray<T> &other)
{
  assert(this != &other);

  data = 0;
  data_area = 0;

  *this = other;
}

template<class T>
_DMatrix<T> _DMatrix<T>::transpose() const
{
    _DMatrix<T> result(_cols, _rows);

    T *in_cp = data_ptr();
    for(int i=0; i<_rows; i++)
      {
	//	T *out_cp = result.data[0]+i;
	//
	for(int j=0; j<_cols; j++) //, out_cp += _rows)
	  result[j][i] = *(in_cp++);
      }
    
    return result;
}


template<class T>
_DMatrix<T> _DMatrix<T>::rotate_180() const
{
    _DMatrix<T> result(_rows, _cols);

    int sz = _rows * _cols;
    const T *in_cp = data_ptr() + sz - 1;
    T *out_cp = result.data_ptr();

    for(int i=0; i<sz; i++)
      *(out_cp++) = *(in_cp--);
    
    return result;
}


template<class T>
_DMatrix<T> _DMatrix<T>::flip_horizontal() const
{
    _DMatrix<T> result(_rows, _cols);

    T *out_cp = result.data_ptr();

    for(int i=0; i<_rows; i++)
      {
	const T *in_cp = data[i] + _cols - 1;

	for(int j=0; j<_cols; j++)
	  *(out_cp++) = *(in_cp--);
      }
    
    return result;
}


template<class T>
_DMatrix<T> _DMatrix<T>::flip_vertical() const
{
    _DMatrix<T> result(_rows, _cols);

    T *out_cp = result.data_ptr();

    for(int i=0; i<_rows; i++)
      {
	const T *in_cp = data[_rows-i-1];

	for(int j=0; j<_cols; j++)
	  *(out_cp++) = *(in_cp++);
      }
    
    return result;
}



template<class T>
_DMatrix<T> operator+(T value, const _DMatrix<T> &other)
{
    return(other + value);
}

template<class T>
_DMatrix<T> operator-(T value, const _DMatrix<T> &other)
{
    _DMatrix<T> result(other.rows(), other.cols());

    T *cp1 = other.data_area, *cp_out = result.data_area;
    for(int i=0; i<other.rows(); i++)
        for(int j=0; j<other.cols(); j++, cp1++, cp_out++)
            *cp_out = -*cp1 + value;

    return result;
}

template<class T>
istream &operator>>(istream &is, _DMatrix<T> &matrix)
{
  int _rows=0, _cols=0;
  is >> _rows >> _cols;

  matrix = _DMatrix<T>(_rows, _cols);

  for(int i=0; i<matrix.rows(); i++)
    for(int j=0; j<matrix.cols(); j++)
      is >> matrix[i][j];
  
  return is;
}

template<class T>
ostream &operator<<(ostream &os, const _DMatrix<T> &matrix)
{
  ofstream os2;
  os2.copyfmt(os);

  //  ios::fmtflags fflags = os.flags();

  os << matrix.rows() << " " << matrix.cols() << endl;

  for(int i=0; i<matrix.rows(); i++)
    {
      for(int j=0; j<matrix.cols(); j++)
	{
	  os.copyfmt(os2);
	  os << matrix[i][j] << " ";
	}
      os << endl;
    }
  
  return os;
}

template<class T2>
void fread(_DMatrix<T2> &matrix, FILE *fp, bool enforce_size)
{
  int rows, cols, type;

  fread(&rows, 1, sizeof(int), fp);
  fread(&cols, 1, sizeof(int), fp);
  fread(&type, 1, sizeof(int), fp);

  assert(type == sizeof(T2));
  if(enforce_size && (matrix.rows() != rows || matrix.cols() != cols))
    throw(std::string("unexpected matrix size encountered during fread"));

  if(!enforce_size)
    matrix = _DMatrix<T2>(rows, cols);

  fread(matrix.data_area, rows*cols, sizeof(T2), fp);
}

template<class T2>
void fread(_DMatrix<T2> &matrix, FILE *fp)
{
  return fread(matrix, fp, false);
}

template<class T2>
void fwrite(const _DMatrix<T2> &matrix, FILE *fp)
{
  int rows = matrix.rows(), cols = matrix.cols(), type = sizeof(T2);

  fwrite(&rows, 1, sizeof(int), fp);
  fwrite(&cols, 1, sizeof(int), fp);
  fwrite(&type, 1, sizeof(int), fp);

  fwrite(matrix.data_area, rows*cols, sizeof(T2), fp);
}


template<class T>
_DMatrix<T> _DMatrix<T>::LU_factor()
{
    _DMatrix<T> result((*this));

    for(int j=0; j<cols(); j++)
        for(int i=j+1; i<rows(); i++)
        {
            T alpha = result[j][j] / result[j][i];
            _DMatrix<T> this_row = extract_row(i);
            result.set_row(i, this_row-this_row*alpha);
        }
    
    return result; 
}

template<class T>
_DMatrix<T> _DMatrix<T>::extract_row(int row) const
{
    return extract(DRect(row, 0, row, cols()-1));
}

template<class T>
_DMatrix<T> _DMatrix<T>::extract_col(int col) const
{
    return extract(DRect(0, col, rows()-1, col));
}

template<class T>
_DMatrix<T> _DMatrix<T>::extract_rows(const std::vector<int> &row_list) const
{
  _DMatrix<T> result((int) row_list.size(), cols());

  vector<int>::const_iterator iter;
  int out_i=0;
  for(iter = row_list.begin(); iter != row_list.end(); ++iter, ++out_i)
    {
      result.copy_row_from(*this, *iter, out_i);
    }

  return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::extract_cols(const std::vector<int> &col_list) const
{
  _DMatrix<T> result(rows(),(int) col_list.size());

  vector<int>::const_iterator iter;
  int out_j=0;
  for(iter = col_list.begin(); iter != col_list.end(); ++iter, ++out_j)
    {
      result.copy_col_from(*this, *iter, out_j);
    }

  return result;
}




template<class T>
_DMatrix<T> _DMatrix<T>::extract(const DRect &rect) const
{
  if(rect.top() < 0 || rect.left() < 0 || rect.right() >= cols() || rect.bottom() >= rows())
    {
      ostringstream oss;
      oss << "in _DMatrix::extract() : couldn't extract " << rect << " from " << DRect(DPoint(0,0), size()) << endl;
      throw oss.str();
    }

  _DMatrix<T> result(rect.height(), rect.width());
		      
  if(result.rows() == 0 || result.cols() == 0)
    return result;
  
  T *cp = result.data_ptr();
  
  int rect_bottom = rect.bottom(), rect_right = rect.right();
  for(int i=rect.top(); i <= rect_bottom; i++)
    {
      const T *in_cp = (*this)[i];

      for(int j=rect.left(); j <= rect_right; j++, cp++)
	*cp = in_cp[j];
    }

  return result;
}

template<class T>
void _DMatrix<T>::set_submatrix(const DPoint &pt, const _DMatrix<T> &in)
{
  if(!(pt.row() + in.rows() <= rows() && pt.row() >= 0) ||
     !(pt.col() + in.cols() <= cols() && pt.col() >= 0))
    {
      ostringstream oss;
      oss << "in _DMatrix::set_submatrix() : couldn't do set_submatrix at " << pt << " on matrix of size " << size() << " with matrix of size " << in.size() << endl;
      throw oss.str();
    }
  
  const T *in_ptr = in.data_ptr();
  
  for(int i=0; i<in.rows(); i++)
    {
      T *out_ptr = (*this)[i+pt.row()] + pt.col();
      
      for(int j=0; j<in.cols(); j++)
	out_ptr[j]=*(in_ptr++);
    }
}

template<class T>
void _DMatrix<T>::set_submatrix(const DRect &rect, T val)
{
  assert(rect.top() >= 0 && rect.left() >= 0 && rect.bottom() < rows() && rect.right() < cols());
  
  for(int i=rect.top(); i <= rect.bottom(); i++)
    {
      T *out_ptr = (*this)[i];
      
      for(int j=rect.left(); j <= rect.right(); j++)
	out_ptr[j] = val;
    }
}


template<class T>
void _DMatrix<T>::set_row(int row, const _DMatrix<T> &in)
{
  assert(in.cols() == cols());
  assert(row >= 0 && row < rows());

  copy_row_from(in, 0, row);
}

template<class T>
void _DMatrix<T>::set_row(int row, const std::vector<T> &in)
{
  assert(int(in.size()) == cols());
  assert(row >= 0 && row < rows());

  T *out_cp = (*this)[row];
  typename std::vector<T>::const_iterator iter;
  for(iter = in.begin(); iter != in.end(); ++iter, ++out_cp)
    *out_cp = *iter;
}

template<class T>
void _DMatrix<T>::set_col(int col, const _DMatrix<T> &in)
{
  assert(in.rows() == rows());
  assert(col >= 0 && col < cols());

  copy_col_from(in, 0, col);
}

template<class T>
_DMatrix<T> _DMatrix<T>::reshape(int new_rows, int new_cols)
{
    assert(new_rows * new_cols == rows() * cols());

    _DMatrix<T> result(new_rows, new_cols);

    T *in_cp = data_ptr();
    T *out_cp = result.data_ptr();

    for(int i=0; i<new_rows*new_cols; i++, out_cp++, in_cp++)
       *out_cp = *in_cp;

    return result;
}

template<class T>
_DMatrix<T> sqrt(const _DMatrix<T> &m1)
{
  _DMatrix<T> result(m1.rows(), m1.cols());

  T *in_cp = m1.data_ptr(), *out_cp = result.data_ptr();
  for(int i=0; i<m1.rows(); i++)
    for(int j=0; j<m1.cols(); j++, in_cp++, out_cp++)
      *out_cp = T(sqrt((double) *in_cp));

  return result;
}

template<class T>
_DMatrix<T> sqr(const _DMatrix<T> &m1)
{
  _DMatrix<T> result(m1.rows(), m1.cols());

  T *in_cp = m1.data_ptr(), *out_cp = result.data_ptr();
  for(int i=0; i<m1.rows(); i++)
    for(int j=0; j<m1.cols(); j++, in_cp++, out_cp++)
      *out_cp = (*in_cp) * (*in_cp);

  return result;
}

template<class T>
void _DMatrix<T>::sqr_ip()
{
  T *in_cp = data_ptr();
  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++, in_cp++)
      *in_cp = (*in_cp) * (*in_cp);

  return;
}

template<class T>
_DMatrix<T> pointwise_multiply(const _DMatrix<T> &m1, const _DMatrix<T> &m2)
{
  assert(same_size(m1, m2));

  _DMatrix<T> result(m1.rows(), m2.cols());

  const T *cp1 = m1.data_ptr(), *cp2 = m2.data_ptr();
  T *out_cp = result.data_ptr();

  int count = m1.rows() * m1.cols();
  for(int i=0; i<count; i++)
    out_cp[i] = cp1[i] * cp2[i];

  return result;
}

template<class T>
_DMatrix<T> pointwise_divide(const _DMatrix<T> &m1, const _DMatrix<T> &m2)
{
  assert(same_size(m1, m2));

  _DMatrix<T> result(m1.rows(), m2.cols());

  for(int i=0; i<m1.rows(); i++)
    for(int j=0; j<m1.cols(); j++)
      result[i][j] = m1[i][j] / m2[i][j];

  return result;
}

template<class T>
_DMatrix<T> pointwise_min(const _DMatrix<T> &m1, const _DMatrix<T> &m2)
{
  assert(same_size(m1, m2));

  _DMatrix<T> result(m1.rows(), m2.cols());

  T *cp1 = m1.data_ptr(), *cp2 = m2.data_ptr(), *cpr = result.data_ptr();
  for(int i=0; i<m1.rows(); i++)
    for(int j=0; j<m1.cols(); j++, cp1++, cp2++, cpr++)
      if(*cp1 < *cp2)
	*cpr = *cp1;
      else
	*cpr = *cp2;

  return result;
}

template<class T>
_DMatrix<T> pointwise_max(const _DMatrix<T> &m1, const _DMatrix<T> &m2)
{
  assert(same_size(m1, m2));

  _DMatrix<T> result(m1.rows(), m2.cols());

  T *cp1 = m1.data_ptr(), *cp2 = m2.data_ptr(), *cpr = result.data_ptr();
  for(int i=0; i<m1.rows(); i++)
    for(int j=0; j<m1.cols(); j++, cp1++, cp2++, cpr++)
      if(*cp1 > *cp2)
	*cpr = *cp1;
      else
	*cpr = *cp2;

  return result;
}

template<class T>
_DMatrix<T> pointwise_min(const _DMatrix<T> &m1, T val)
{
  _DMatrix<T> result(m1);

  T *cpr = result.data_ptr();
  int sz  = m1.rows() * m1.cols();
  for(int i=0; i<sz; i++)
    if(val < cpr[i])
      cpr[i] = val;

  return result;
}

template<class T>
_DMatrix<T> pointwise_max(const _DMatrix<T> &m1, T val)
{
  _DMatrix<T> result(m1);

  T *cpr = result.data_ptr();
  int sz  = m1.rows() * m1.cols();
  for(int i=0; i<sz; i++)
    if(val > cpr[i])
      cpr[i] = val;

  return result;
}

template<class T>
_DMatrix<T> fabs(const _DMatrix<T> &m)
{
  _DMatrix<T> result(m.rows(), m.cols());

  const T *in = m.data_ptr();
  T *out = result.data_ptr();
  int sz = m.rows() * m.cols();

  for(int i=0; i<sz; i++)
    out[i] = T(fabs(in[i]));

  return result;
}

template<class T>
_DMatrix<T> log(const _DMatrix<T> &m)
{
  _DMatrix<T> result(m.rows(), m.cols());

  const T *in = m.data_ptr();
  T *out = result.data_ptr();
  int sz = m.rows() * m.cols();

  for(int i=0; i<sz; i++)
    out[i] = T(log((double) in[i]));

  return result;
}

template<class T>
_DMatrix<T> exp(const _DMatrix<T> &m, bool fast)
{
  _DMatrix<T> result(m.rows(), m.cols());

  const T *in = m.data_ptr();
  T *out = result.data_ptr();
  int sz = m.rows() * m.cols();

  if(!fast)
    for(int i=0; i<sz; i++)
      out[i] = T(exp((double) in[i]));
  else
    for(int i=0; i<sz; i++)
      out[i] = T(fast_exp(in[i]));

  return result;
}

template<class T>
_DMatrix<T> exp(const _DMatrix<T> &m)
{
  return exp(m, false);
}

// compute mean of all entries in matrix
template<class T>
T _DMatrix<T>::mean() const
{
  return sum() / (rows() * cols());
}


// compute sum of all entries in matrix
template<class T>
T _DMatrix<T>::sum() const
{
  T sum = 0;

  const T *in_ptr = data_ptr();
  int sz = rows() * cols();

  for(int i=0; i<sz; i++)
    sum += in_ptr[i];

  return sum;
}

// compute sum of all entries in each matrix row
template<class T>
_DMatrix<T> _DMatrix<T>::sum_rows() const
{
  _DMatrix<T> result(rows(), 1);

  const T *in_ptr = data_ptr();
  for(int i=0; i<rows(); i++)
    {
      T row_sum = 0;
      for(int j=0; j<cols(); j++, in_ptr++)
	row_sum += *in_ptr;

      result[i][0] = row_sum;
    }
      
  return result;
}


// compute sum of all entries in each matrix column
template<class T>
_DMatrix<T> _DMatrix<T>::sum_cols() const
{
  _DMatrix<T> result(1, cols());
  result = 0;

  T *res_ptr = result[0];
  for(int i=0; i<rows(); i++)
    {
      const T *row_ptr = (*this)[i];

      for(int j=0; j<cols(); j++)
	res_ptr[j] += row_ptr[j];
    }
      
  return result;
}



// compute median of all entries in matrix
template<class T>
T _DMatrix<T>::median() const
{
  vector<T> vect;

  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++)
      vect.push_back((*this)[i][j]);

  sort(vect.begin(), vect.end());

  return vect[vect.size()/2];
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator==(T value)
{
  _DMatrix<T> result(rows(), cols());

  int sz = rows() * cols();

  const T *in_ptr = data_ptr();
  T *out_ptr = result.data_ptr();

  for(int i=0; i<sz; i++)
    out_ptr[i] = (in_ptr[i] == value);

  return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator<(T value)
{
  _DMatrix<T> result(rows(), cols());

  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++)
      result[i][j] = (*this)[i][j] < value;

  return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator>(T value)
{
  _DMatrix<T> result(rows(), cols());

  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++)
      result[i][j] = (*this)[i][j] > value;

  return result;
}


// computes the mean of *each column* of the input matrix.
//
template<class T>
_DMatrix<T> _DMatrix<T>::means() const
{
  _DMatrix<T> result(1, cols());
  T *res_ptr = result.data_ptr();
  const T *in_ptr = data_ptr();

  for(int j=0; j<cols(); j++)
  {
    const T *in_cp = in_ptr+j;
    T mean=0;

    for(int i=0; i<rows(); i++, in_cp += cols())
      mean += *(in_cp);

    res_ptr[j] = mean / (rows());
  }

  return result;
}


template<class T>
_DMatrix<T> _DMatrix<T>::trimmed_means(int leaveout_count) const
{
  _DMatrix<T> sorted = (*this);
  sorted.sort_cols();

  assert(sorted.rows() > leaveout_count * 2);

  return (sorted.extract(DRect(leaveout_count, 0, rows()-1-leaveout_count, cols()-1))).means();
}


template<class T>
void _DMatrix<T>::sort_cols()
{
  assert(0);

  //  (*this) = transpose();
  //  sort_rows();
  //  (*this) = transpose();
}



template<class T>
int compar_sortrows(const void *t1, const void *t2)
{
  T i1 = *((T *) t1), i2 = *((T *) t2);

  if(i1 > i2)
    return 1;
  else if(i1 < i2)
    return -1;
  else
    return 0;
}

// FIXME: get rid of global variable
int __sort_colcount = 0;

template<class T>
int compar_sortwholerows(const void *t1, const void *t2)
{
  T *i1 = (T *) t1, *i2 = (T *) t2;

  for(int i=0; i<__sort_colcount; i++)
    if(i1[i] < i2[i])
      return -1;
    else if(i1[i] > i2[i])
      return 1;

  return 0;
}


template<class T>
void _DMatrix<T>::sort_rows()
{
  for(int i=0; i<rows(); i++)
    qsort((*this)[i], cols(), sizeof(T), compar_sortrows<T>);
}

template<class T>
void _DMatrix<T>::sort_wholerows()
{
  __sort_colcount = cols();
  qsort(data_ptr(), rows(), sizeof(T) * cols(), compar_sortwholerows<T>);
}


template<class T>
_DMatrix<T> _DMatrix<T>::unique_rows(std::vector<int> &row_counts) const
{
  row_counts.clear();

  _DMatrix<T> tmp(*this);
  tmp.sort_wholerows();

  int cur_count = 1;
  for(int i=1; i<rows(); i++)
    if(compar_sortwholerows<T>((void *)((*this)[i-1]), (void *)((*this)[i])) != 0)
      {
	row_counts.push_back(cur_count);
	cur_count = 1;
      }
    else
      cur_count++;

  row_counts.push_back(cur_count);

  _DMatrix<T> result((int) row_counts.size(), cols());
  int cum_row = 0;
  for(int i=0; i<(int)row_counts.size(); i++)
    {
      result.copy_row_from((*this), cum_row, i);
      cum_row += row_counts[i];
    }

  return result;
}


template<class T>
T _DMatrix<T>::max(int &out_row, int &out_col) const
{
  const T *in_ptr = data_ptr();

  int max_i = 0;
  T max_element = in_ptr[0];
  int sz = rows() * cols();

  for(int i=0; i < sz; i++)
    if(max_element < in_ptr[i])
      {
	max_element = in_ptr[i];
	max_i = i;
      }

  out_row = max_i / cols(), out_col = max_i % cols();

  return max_element;
}


template<class T>
T _DMatrix<T>::min(int &out_row, int &out_col) const
{
  const T *cp = data_ptr();

  out_row = out_col = 0;
  T min_element = *cp;
  
  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++, cp++)
      if(min_element > *cp)
	{
	  min_element = *cp;
	  out_row = i, out_col = j;
	}

  return min_element;
}

// Returns the same value as max(). If there is a unique global maximum,
// also returns the same point as max(); otherwise returns a randomly-chosen
// maximum.
template<class T>
T _DMatrix<T>::max_random(DPoint &out_location) const
{
  const T *cp = data_ptr();

  int out_row = 0, out_col = 0;
  T max_element = *cp;
  int max_count = 0;
  
  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++, cp++)
      if(max_element < *cp)
	  max_element = *cp, out_row = i, out_col = j, max_count = 1;
      else if(max_element == *cp)
	{
	  max_count++;

	  // now possibly switch max_element to the new maximum,
	  // with probability 1/max_count. This is the same as choosing among
	  // the maximums with a uniform distribution.
	  if(drand48() * max_count < 1.0)
	    out_row = i, out_col = j;
	}

  out_location = DPoint(out_row, out_col);
  return max_element;
}

template<class T>
bool _DMatrix<T>::is_any_nan() const
{
  int sz = rows() * cols();
  for(int i=0; i<sz; i++)
    if(isnan(data_area[i]))
      return true;

  return false;
}


template<class T>
bool _DMatrix<T>::is_any_inf() const
{
  int sz = rows() * cols();
  for(int i=0; i<sz; i++)
    if(isinf(data_area[i]))
      return true;

  return false;
}


// assumes each row is an observation, each column is
// a feature.
template<class T>
_DMatrix<T> _DMatrix<T>::covariance() const
{
  _DMatrix<T> mean_vector = means();

  return covariance(mean_vector);
}


template<class T>
_DMatrix<T> _DMatrix<T>::covariance(const _DMatrix<T> &mean_vector) const
{
  _DMatrix<T> cov(cols(), cols());
  const T *mean_vec = mean_vector.data_ptr();

  // handle special case of cols=2 separately, for speed
  if(cols() == 2)
    {
      T var00=0, var01=0, var11=0;

      const T *row_cp = data_ptr();
      for(int k=0; k<rows(); k++, row_cp += cols())
	{
	  T a0 = (row_cp[0] - mean_vec[0]);
	  T a1 = (row_cp[1] - mean_vec[1]);

	  var00 += a0 * a0;
	  var11 += a1 * a1;
	  var01 += a0 * a1;
	}

      cov[0][0] = var00 / (rows()-1);
      cov[1][1] = var11 / (rows()-1);
      cov[0][1] = cov[1][0] = var01 / (rows()-1);
    }
  else
    {
      for(int i=0; i<cols(); i++)
	for(int j=0; j<cols(); j++)
	  {
	    T var = 0;
	    
	    const T *row_cp = data_ptr();
	    for(int k=0; k<rows(); k++, row_cp += cols())
	      var += (row_cp[i]-mean_vec[i]) * (row_cp[j]-mean_vec[j]); // was mean_vec[i]
	    
	    cov[i][j] = var / (rows()-1);
	  }
    }
  return cov;
}


template<class T>
void _DMatrix<T>::swap_rows(int row1, int row2)
{
  // ideally, we'd just have to swap the pointers, except that
  // some routines assume that the matrix is arranged
  // contiguously in row-major order. So we have to move
  // memory around also (unfortunately).

  if(row1 == row2)
    return;

  _DMatrix<T> temp(1, cols());
  temp.copy_row_from(*this, row1, 0);  
  copy_row_from(*this, row2, row1);    
  copy_row_from(temp, 0, row2);
}

template<class T>
void _DMatrix<T>::copy_row_from(const _DMatrix<T> &other, int src_row, int dest_row)
{
  assert(other.cols() == cols());
  assert(src_row >= 0 && src_row < other.rows());
  assert(dest_row >= 0 && dest_row < rows());

  memcpy((*this)[dest_row], other[src_row], sizeof(T) * other.cols());
}

template<class T>
void _DMatrix<T>::copy_row_from(const DPoint &other, int dest_row)
{
  assert(cols() == 2);

  (*this)[dest_row][0] = other.row();
  (*this)[dest_row][1] = other.col();
}

template<class T>
void _DMatrix<T>::copy_col_from(const _DMatrix<T> &other, int src_col, int dest_col)
{
  assert(other.rows() == rows());
  assert(src_col >= 0 && src_col < other.cols());
  assert(dest_col >= 0 && dest_col < cols());

  T *in_ptr = &(other.data_ptr()[src_col]);
  T *out_ptr = &(data_ptr()[dest_col]);

  for(int i=0; i<rows(); i++, in_ptr += other.cols(), out_ptr += cols())
    *out_ptr = *in_ptr;
}



template<class T>
vector<DPoint> _DMatrix<T>::sample_probabilities(int sample_count) const
{
  int s=0;
  vector<T> samples(sample_count);
  vector<DPoint> result(sample_count);

  T _sum = sum();
  for(s=0; s<sample_count; s++)
    samples[s] = T(drand48()) * _sum;
  
  sort(samples.begin(), samples.end());
  
  _sum = 0;
  s=0;

  T *prob_map_cp = data_ptr();
  T this_sample = samples[s];

  for(int i=0, cp=0; i<rows(); i++)
    for(int j=0; j<cols(); j++, cp++)
      {
	_sum += prob_map_cp[cp];
	while(_sum > this_sample)
	  {
	    result[s].row(i), result[s].col(j); 
	    
	    s++;
	    this_sample = samples[s];
	    if(s >= sample_count)
	      goto done;
	  }
      }
  
 done:
  
  static int count = 0;
  count++;
  
  return result;
}


template<class T>
vector<DPoint> _DMatrix<T>::sample_likelihoods(int sample_count, bool fast) const
{
  _DMatrix<T> prob_map = exp((*this) - max(), fast);
  
  //  prob_map = pointwise_max<T>(prob_map, T(1e-20));

  return prob_map.sample_probabilities(sample_count);
}


#ifdef GSL_SUPPORT
template<class T>
gsl_matrix *DMatrix_to_gsl(const _DMatrix<T> &dm)
{
  gsl_matrix *gm = gsl_matrix_alloc(dm.rows(), dm.cols());

  T *cp = dm.data_ptr();
  for(int i=0; i<dm.rows(); i++)
    for(int j=0; j<dm.cols(); j++, cp++)
      gsl_matrix_set(gm, i, j, *cp);

  return gm;
}


template<class T>
_DMatrix<T> gsl_to_DMatrix(gsl_matrix *gm)
{
  _DMatrix<T> dm((int) gm->size1, (int) gm->size2);

  T *cp = dm.data_ptr();
  for(int i=0; i<dm.rows(); i++)
    for(int j=0; j<dm.cols(); j++, cp++)
      *cp = T(gsl_matrix_get(gm, i, j));

  return dm;
}

template<class T>
_DMatrix<T> gsl_to_DMatrix(gsl_vector *gm)
{
  _DMatrix<T> dm(1, (int) gm->size);

  T *cp = dm.data_ptr();
  for(int i=0; i<dm.cols(); i++, cp++)
    *cp = T(gsl_vector_get(gm, i));

  return dm;
}

#endif

template<class T>
_DMatrix<T> _DMatrix<T>::inverse() const
{
#ifdef GSL_SUPPORT
  assert(rows() == cols());
  if(rows() == 0)
    return _DMatrix<T>();

  gsl_matrix *gm = DMatrix_to_gsl(*this);
  gsl_permutation *p = gsl_permutation_alloc(rows());;
  int signum;

  gsl_linalg_LU_decomp(gm, p, &signum);

  gsl_matrix *inverse = gsl_matrix_alloc(rows(), cols());
  gsl_linalg_LU_invert(gm, p, inverse);

  _DMatrix<T> dm = gsl_to_DMatrix<T>(inverse);
  gsl_permutation_free(p);
  gsl_matrix_free(gm);
  gsl_matrix_free(inverse);

  return dm;
#else
  throw string("no gsl support");
#endif
}

template<class T>
T _DMatrix<T>::determinant() const
{
#ifdef GSL_SUPPORT
  assert(rows() == cols());

  if(rows() == 0)
    return T(1);
  else if(rows() == 1)
    return (*data_area);
  else if(rows() == 2)
    return (data_area[0] * data_area[3] - data_area[1] * data_area[2]);

  gsl_matrix *gm = DMatrix_to_gsl(*this);
  gsl_permutation *p = gsl_permutation_alloc(rows());;
  int signum;

  gsl_linalg_LU_decomp(gm, p, &signum);

  T det = T(gsl_linalg_LU_det(gm, signum));

  _DMatrix<T> dm = gsl_to_DMatrix<T>(gm);
  gsl_permutation_free(p);
  gsl_matrix_free(gm);

  return det;
#else
  throw string("no gsl support");
#endif
}

template<class T>
_DMatrix<T> _DMatrix<T>::eigen(_DMatrix<T> &d_eigvec) const
{
#ifdef GSL_SUPPORT
  assert(rows() == cols());

  gsl_matrix *gm = DMatrix_to_gsl(*this);
  gsl_vector *eigval = gsl_vector_alloc(rows());
  gsl_matrix *eigvec = gsl_matrix_alloc(rows(), rows());
  gsl_eigen_symmv_workspace *worksp = gsl_eigen_symmv_alloc(rows());

  gsl_eigen_symmv(gm, eigval, eigvec, worksp);

  gsl_eigen_symmv_sort (eigval, eigvec, GSL_EIGEN_SORT_ABS_ASC);

  d_eigvec = gsl_to_DMatrix<T>(eigvec);
  _DMatrix<T> d_eigval = gsl_to_DMatrix<T>(eigval);

  gsl_matrix_free(gm);
  gsl_matrix_free(eigvec);
  gsl_vector_free(eigval);
  gsl_eigen_symmv_free(worksp);

  return d_eigval;
#else
  throw string("no gsl support");
#endif
}



template<class T>
void _DMatrix<T>::svd(_DMatrix<T> &u, _DMatrix<T> &s, _DMatrix<T> &v) const
{
#ifdef GSL_SUPPORT
  gsl_matrix *gm = DMatrix_to_gsl(*this);
  gsl_vector *S = gsl_vector_alloc(cols());
  gsl_matrix *V = gsl_matrix_alloc(cols(), cols());
  gsl_vector *work = gsl_vector_alloc(cols());

  gsl_linalg_SV_decomp(gm, V, S, work);

  u = gsl_to_DMatrix<T>(gm);
  v = gsl_to_DMatrix<T>(V);
  s = gsl_to_DMatrix<T>(S);

  gsl_matrix_free(gm);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);
  return;
#else
  throw string("no gsl support");
#endif
}


template<class T2>
_DMatrix<T2> horiz_concat(const _DMatrix<T2> &m1, const _DMatrix<T2> &m2)
{
  assert(m1.rows() == m2.rows());

  _DMatrix<T2> result(m1.rows(), m1.cols() + m2.cols());

  result.set_submatrix(DPoint(0, 0), m1);
  result.set_submatrix(DPoint(0, m1.cols()), m2);

  return result;
}


template<class T2>
_DMatrix<T2> vert_concat(const _DMatrix<T2> &m1, const _DMatrix<T2> &m2)
{
  assert(m1.cols() == m2.cols());

  _DMatrix<T2> result(m1.rows() + m2.rows(), m1.cols());

  result.set_submatrix(DPoint(0, 0), m1);
  result.set_submatrix(DPoint(m1.rows(), 0), m2);

  return result;
}

template<class T>
std::vector<T> _DMatrix<T>::extract_row_as_vector(int row) const
{
  assert(row >= 0 && row < rows());

  std::vector<T> result(cols());

  const T *row_ptr = (*this)[row];
  for(int j=0; j<cols(); j++)
    result[j] = row_ptr[j];

  return result;
}





template<class T>
void _DMatrix<T>::search_and_replace_ip(T src, T dest)
{
  T *cp = data_ptr();
  int sz = rows() * cols();
  for(int i=0; i<sz; i++)
    if(cp[i] == src) cp[i] = dest;

  return;
}

template<class T2> 
_DMatrix<T2> operator>(const _DMatrix<T2> &m1, T2 val)
{
  _DMatrix<T2> result(m1);

  T2 *cp = result.data_ptr();
  int sz = m1.rows() * m1.cols();
  for(int i=0; i<sz; i++)
    if(cp[i] > val)
      cp[i] = 1;
    else
      cp[i] = 0;

  return result;
}


template<class T2> 
_DMatrix<T2> operator<(const _DMatrix<T2> &m1, T2 val)
{
  _DMatrix<T2> result(m1);

  T2 *cp = result.data_ptr();
  int sz = m1.rows() * m1.cols();
  for(int i=0; i<sz; i++)
    if(cp[i] < val)
      cp[i] = 1;
    else
      cp[i] = 0;

  return result;
}


template<class T2, class T3>
void change_type(const _DMatrix<T2> &in, _DMatrix<T3> &result)
{
  if(!(in.rows() == result.rows() && in.cols() == result.cols()))
    result = _DMatrix<T3>(in.rows(), in.cols());

  int sz = in.rows() * in.cols();
  T3 *res_ptr = result.data_ptr();
  const T2 *in_ptr = in.data_ptr();

  for(int i=0; i<sz; i++)
    res_ptr[i] = (T3) in_ptr[i];
}



#define DECLARE(x) \
  template class _DTwoDimArray<x>; \
  template class _DTwoDimArray<std::complex<x> >; \
  template class _DMatrix<x>; \
  template _DMatrix<x> operator-(x value, const _DMatrix<x> &other); \
  template _DMatrix<x> pointwise_multiply(const _DMatrix<x> &m1, const _DMatrix<x> &m2); \
  template _DMatrix<x> pointwise_divide(const _DMatrix<x> &m1, const _DMatrix<x> &m2); \
  template ostream &operator<<(ostream &os, const _DMatrix<x> &matrix); \
  template istream &operator>>(istream &is, _DMatrix<x> &matrix); \
  template _DMatrix<x> operator-(const _DMatrix<x> &m); \
  template _DMatrix<x> pointwise_min(const _DMatrix<x> &, const _DMatrix<x> &); \
  template _DMatrix<x> pointwise_max(const _DMatrix<x> &, const _DMatrix<x> &); \
  template _DMatrix<x> pointwise_min(const _DMatrix<x> &, x val); \
  template _DMatrix<x> pointwise_max(const _DMatrix<x> &, x val); \
  template _DMatrix<x> sqr(const _DMatrix<x> &); \
  template _DMatrix<x> log(const _DMatrix<x> &); \
  template _DMatrix<x> exp(const _DMatrix<x> &); \
  template _DMatrix<x> exp(const _DMatrix<x> &, bool); \
  template _DMatrix<x> sqrt(const _DMatrix<x> &); \
  template _DMatrix<x> fabs(const _DMatrix<x> &); \
  template _DMatrix<x> horiz_concat(const _DMatrix<x> &m1, const _DMatrix<x> &m2); \
  template _DMatrix<x> vert_concat(const _DMatrix<x> &m1, const _DMatrix<x> &m2); \
  template _DMatrix<x> operator<(const _DMatrix<x> &m1, x val); \
  template _DMatrix<x> operator>(const _DMatrix<x> &m1, x val); \
  template _DMatrix<x> operator*(x value, const _DMatrix<x> &other); \
  template void fread(_DMatrix<x> &, FILE *fp, bool enforce_size); \
  template void fread(_DMatrix<x> &, FILE *fp); \
  template void fwrite(const _DMatrix<x> &, FILE *fp);  \
  template bool same_size(const _DMatrix<x> &m1, const _DMatrix<x> &m2);

#define DECLARE2(x,y) \
  template void change_type(const _DMatrix<x> &m1, _DMatrix<y> &m2); 

DECLARE(double)
DECLARE(short)
DECLARE(int)
DECLARE(float)
DECLARE(char)
DECLARE(unsigned char)


  DECLARE2(double, double)
  DECLARE2(float, float)
  DECLARE2(double, unsigned char)
  DECLARE2(double, int)
  DECLARE2(double, short)
  DECLARE2(double, float)
  DECLARE2(float, unsigned char)
  DECLARE2(float, int)
  DECLARE2(float, short)
  DECLARE2(float, double)
  DECLARE2(int, unsigned char)
  DECLARE2(int, float)
  DECLARE2(int, short)
  DECLARE2(int, double)
  DECLARE2(short, float)
  DECLARE2(short, double)
  DECLARE2(short, int)
  DECLARE2(short, unsigned char)
  DECLARE2(unsigned char, int)
  DECLARE2(unsigned char, float)
  DECLARE2(unsigned char, short)
  DECLARE2(unsigned char, double)

