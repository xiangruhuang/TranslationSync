#include <vector>
#include <numeric>
#include <DGaussianModel.h>

using namespace std;

template<class T>
DGaussianModel<T>::DGaussianModel(const _DMatrix<T> &__covariance, const _DMatrix<T> &__mean, bool reverse_dims) :
  _covariance(__covariance), _mean(__mean)
{
  parameter_sanity_check();
  
  if(reverse_dims && _covariance.cols() == 2)
    {
      T tmp = _covariance[0][0];
      _covariance[0][0] = _covariance[1][1];
      _covariance[1][1] = tmp;
      
      tmp = _mean[0][0];
      _mean[0][0] = _mean[1][0];
      _mean[1][0] = tmp;
    }
  
  constants_clean = false;
}

template<class T>
DGaussianModel<T>::DGaussianModel(const DGaussianModel<T> &other)
{
  *this = other;
}

template<class T>
DGaussianModel<T>::DGaussianModel(const _DMatrix<T> &data)
{
  learn(data);
}

template<class T>
DGaussianModel<T>::DGaussianModel() 
{
  constants_clean = false;
}

template<class T>
void DGaussianModel<T>::learn(const _DMatrix<T> &data)
{
  constants_clean = false;
  _mean = data.means().transpose();
  _covariance = data.covariance(_mean);
}

template<class T>
T DGaussianModel<T>::get_data_likelihood(const _DMatrix<T> &data, std::vector<T> &likelihoods)
{
  if(data.rows() == 0 ||  data.cols() == 0)
    return 0;

  compute_constants();

  if(fabs(determinant_of_covariance()) < 1e-10)
    throw std::string("nearly singular");

  likelihoods = vector<T>(data.rows());

  const T *data_cp = data[0];

  if(data.cols() == 2)
    {
      T i_cov_11 = _inverse_covariance[0][0] / 2.0, i_cov_22 = _inverse_covariance[1][1] / 2.0;
      T i_cov_12 = _inverse_covariance[0][1];
      T mean_row = _mean[0][0], mean_col = _mean[1][0];

      typename vector<T>::iterator likelihood_iter = likelihoods.begin();
      for(; likelihood_iter != likelihoods.end(); ++likelihood_iter)
	{
	  T a = *(data_cp++) - mean_row;
	  T b = *(data_cp++) - mean_col;
	      
	  *likelihood_iter = log_c - ( a * a * i_cov_11 + a * b * i_cov_12 + b * b * i_cov_22);
	}
    }
  else
    {
      for(int i=0; i<data.rows(); i++)
	{
	  _DMatrix<T> _data = (data.extract_row(i)).transpose();
	  likelihoods[i] = log_c - ((_data - _mean).transpose() * _inverse_covariance * (_data - _mean) / 2.0)[0][0];
	}
    }

  return accumulate(likelihoods.begin(), likelihoods.end(), T(0.0));
}


template<class T>
T DGaussianModel<T>::get_data_likelihood(T data)
{
  compute_constants();
  assert(dimensionality() == 1);

  if(fabs(determinant_of_covariance()) < 1e-10)
    throw std::string("nearly singular");
  
  double __mean = _mean[0][0];
  T likelihood = log_c - ((data - __mean) * _inverse_covariance[0][0] * (data - __mean) / 2.0);

  return likelihood;
}



template<class T>
_DMatrix<T> DGaussianModel<T>::get_likelihood_pointwise(const _DMatrix<T> &data)
{
  _DMatrix<T> result(data.rows(), data.cols());

  compute_constants();
  assert(dimensionality() == 1);

  if(fabs(determinant_of_covariance()) < 1e-10)
    throw std::string("nearly singular");
  
  double __mean = _mean[0][0];
  double __invcov = _inverse_covariance[0][0];
  
  int sz = data.rows() * data.cols();
  const T *in_ptr = data[0];
  T *out_ptr = result[0];

  for(int i=0; i<sz; i++)
    out_ptr[i] = log_c - ((in_ptr[i] - __mean) * __invcov * (in_ptr[i] - __mean) / 2.0);

  return result;
}


template<class T>
T DGaussianModel<T>::get_data_likelihood(const _DMatrix<T> &data)
{
  vector<T> likelihoods;

  return get_data_likelihood(data, likelihoods);
}

template<class T>
class like_item
{
 public:
  T likelihood;
  int index;
};

template<class T>
inline bool operator<(const like_item<T> &i1, const like_item<T> &i2)
{
  return i1.likelihood < i2.likelihood;
}



// trim percent is the fraction **you want to trim off**
template<class T>
T DGaussianModel<T>::get_data_likelihood_trimmed(const _DMatrix<T> &data, float trim_percent, std::vector<T> &result_likelihoods, _DMatrix<T> &result_data)
{
  vector<T> likelihood_list;

  get_data_likelihood(data, likelihood_list);

  // FIXME
  vector<like_item<T> > likelihoods(data.rows());
  typename vector<like_item<T> >::iterator out_iter = likelihoods.begin();
  typename vector<T>::iterator in_iter = likelihood_list.begin();

  
  for(int i=0 ; out_iter != likelihoods.end(); ++in_iter, ++out_iter, i++)
    {
      out_iter->likelihood = *in_iter;
      out_iter->index = i;
    }

  int worst_count = int(likelihoods.size() * trim_percent);
  int good_count = data.rows() - worst_count;

  nth_element(likelihoods.begin(), likelihoods.begin() + worst_count, likelihoods.end());

  result_data = _DMatrix<T>(good_count, data.cols());
  result_likelihoods = std::vector<T>(good_count);

  T total_likelihood = 0;
  for(int i=worst_count, new_row=0; i<data.rows(); i++, new_row++)
    {
      result_data.copy_row_from(data, likelihoods[i].index, new_row);
      result_likelihoods[new_row] = likelihoods[i].likelihood;
      total_likelihood += likelihoods[i].likelihood;
    }

  return total_likelihood;
}


template<class T>
_DMatrix<T> DGaussianModel<T>::get_likelihood_plane(int rows, int cols)
{
  compute_constants();

  if(fabs(determinant_of_covariance()) < 1e-10)
    throw std::string("nearly singular");

  T i_cov_11 = _inverse_covariance[0][0] / 2.0, i_cov_22 = _inverse_covariance[1][1] / 2.0;
  T i_cov_12 = _inverse_covariance[0][1];
  T mean_row = _mean[0][0], mean_col = _mean[1][0];
      
  _DMatrix<T> result(rows, cols);

  for(int i=0; i<rows; i++)
    {
      T a = i - mean_row;
      T row_pre = log_c - (a * a * i_cov_11);
      T b = 0 - mean_col;

      T *out_cp = result[i];

      for(int j=0; j<cols; j++)
	{
	  out_cp[j] = row_pre - ( a * b * i_cov_12 + b * b * i_cov_22 );
	  b++;
	}
    }

  return result;
}

template<class T>
void DGaussianModel<T>::compute_constants()
{
  if(constants_clean) 
    return;

  // assume one data point
  const T N=1.0;

  T det = _covariance.determinant();
  log_c = -0.5 * _covariance.rows() * N * log(2.0 * M_PI) - 0.5 * N * log(det);
  _inverse_covariance = _covariance.inverse();

  constants_clean = true;
}

template<class T>
void DGaussianModel<T>::parameter_sanity_check()
{
  assert(_covariance.rows() == _covariance.cols());
  // assert(covariance.is_positive_definite();
  // assert(covariance.is_symmetric());
  assert(_mean.rows() == _covariance.rows() && _mean.cols() == 1);
}

template<class T>
DGaussianModel<T> &DGaussianModel<T>::operator=(const DGaussianModel<T> &other)
{
  _covariance = other.covariance();
  _mean = other.mean();
  _inverse_covariance = other._inverse_covariance;
  constants_clean = other.constants_clean;
  log_c = other.log_c;

  return *this;
}

template<class T>
T DGaussianModel<T>::log_constant() 
{
  compute_constants();

  return log_c;
}


template<class T>
_DMatrix<T> DGaussianModel<T>::inverse_covariance() 
{
  compute_constants();
  return _inverse_covariance; 
}


template<class T>
void DGaussianModel<T>::set_mean(const _DMatrix<T> &new_mean)
{
  _mean = new_mean;

  parameter_sanity_check();
}

template<class T>
void DGaussianModel<T>::set_mean(const DPoint &new_mean)
{
  _mean = _DMatrix<T>(2,1);
  _mean[0][0] = new_mean.row();
  _mean[0][1] = new_mean.col();

  parameter_sanity_check();
}

template<class T>
void DGaussianModel<T>::set_covariance(const _DMatrix<T> &new_covariance)
{
  _covariance = new_covariance;
  constants_clean = false;

  parameter_sanity_check();
}

template<class T>
void DGaussianModel<T>::force_diagonal_covariance()
{
  constants_clean = false;

  set_covariance(pointwise_multiply(covariance(), _DMatrix<T>(_covariance.rows(), _covariance.cols(), _DMatrix<T>::identity)));
}


template<class T>
ostream & operator<<(ostream &os, const DGaussianModel<T> &model)
{
  cout << "mean: " << model.mean() << endl;
  cout << "covariance: " << model.covariance() << endl;

  return os;
}


// check whether the covariance in at least one dimension (when
// rotated to be axis oriented) is less than some threshold.
// This function checks whether the smallest eigenvalue is less than
// the threshold.
//
template<class T>
inline bool DGaussianModel<T>::is_covariance_bad(double threshold) const
{
  _DMatrix<T> cov = covariance();

  T a = cov[0][0], b = cov[0][1], d = cov[1][1];

  T low_eigenvalue = 0.5 * ((a+d) - sqrt(4*b*b + (a-d)*(a-d)));

  return low_eigenvalue < threshold;
}


template<class T>
void DGaussianModel<T>::rescale(float scale_factor)
{
  set_mean(mean() * scale_factor);
  set_covariance(covariance() * scale_factor * scale_factor);
}


#define DECLARE(x) \
  template class DGaussianModel<x>; \
  template ostream &operator<<(ostream &os, const DGaussianModel<x> &matrix); 

DECLARE(double)
DECLARE(float)
