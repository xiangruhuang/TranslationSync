#include <DGaussianMixtureModel.h>
#include <DCluster.h>
#include <DMultiDMatrix.h>
#include <float.h>
#include <DMath.h>

using namespace std;

template<class T>
void DGaussianMixtureModel<T>::learn(const _DMatrix<T> &data, int mixture_count, 
				     float delta_thresh, int max_iters, int replicate_count)
{
  KMeansCluster<T> km(mixture_count, replicate_count);

  _DMatrix<int> assignments = km.do_clustering(data);

  //  cout << assignments << endl;

  // There are two special results of do_clustering that must be handled here:
  //
  // (a) it could be that do_clustering produced fewer than mixture_count
  //     clusters, due to the make-up of the data. 
  // (b) it could be that a cluster contains only one (unique) point, in
  //     which case the estimate of the covariance matrix will be singular.
  //

  // handle case (a). This works because do_clustering() gives empty clusters
  //  the highest cluster id's. case(b) is handled in other learn().
  mixture_count = km.get_actual_cluster_count();


  vector<DGaussianModel<T> > initial_models(mixture_count);
  vector<int> counts(mixture_count, 0);
  vector<float> initial_weights(mixture_count, 0);
  
  for(int i=0; i<data.rows(); i++)
    counts[assignments[0][i]]++;

  _DMatrixArray<T> datas(3,mixture_count);
  for(int i=0; i<mixture_count; i++)
    {
      datas.get(i) = _DMatrix<T>(counts[i], 2);
      initial_weights[i] = counts[i] / float(data.rows());
    }

  vector<int> indices(mixture_count, 0);
  for(int i=0; i<data.rows(); i++)
    {
      int assign = assignments[0][i];

      (datas.get(assign))[indices[assign]][0] = data[i][0];
      (datas.get(assign))[indices[assign]][1] = data[i][1];

      indices[assign]++;
    }

  for(int i=0; i<mixture_count; i++)
    {
      //      cout << "-------> " << i << datas.get(i) << endl;

      
            // cout<< "kkkkk " << i << endl;
            // cout<< datas.get(i) << endl;

      // handle case (b) above.
      /*      if(counts[i] == 1)
	{
	  counts[i]++;
	  datas.get(i) = vert_concat(datas.get(i), datas.get(i));
	}
      */
      initial_models[i] = DGaussianModel<T>(datas.get(i));

      if(initial_models[i].is_covariance_bad() || isnan(initial_models[i].log_constant()))
	{
	  _DMatrix<T> cov = initial_models[i].covariance();
	  cov[0][1] = cov[1][0] = 0;
	  // cout << "here" << endl;
	  initial_models[i].set_covariance(cov);

	  if(initial_models[i].is_covariance_bad()  || isnan(initial_models[i].log_constant()))
	    initial_models[i].set_covariance(_DMatrix<T>(2,2,_DMatrix<T>::identity) * 10.0);
	}

      //	initial_models[i] = DGaussianModel<T>(datas.get(i) + _DMatrix<T>(counts[i], 2, _DMatrix<T>::random) * 5.0 - 2.5);

      //             cout<< initial_models[i] << endl;

    }

  learn(data, initial_models, initial_weights, delta_thresh, max_iters);

  return;
}




template<class T>
void DGaussianMixtureModel<T>::learn(const _DMatrix<T> &data, const std::vector<DGaussianModel<T> > &initial_model,
				     const std::vector<float> &initial_weights, float delta_thresh,
				     int max_iters)
{
  int mixture_count = initial_model.size();
  assert((int)initial_weights.size() == mixture_count);

  weights = initial_weights;
  models = initial_model;

  _DMatrix<float> G(data.rows(), mixture_count);
  float last_like = -1e50, this_like = 0;

  _DMatrix<T> s1(1,data.cols());
  T *s1_ptr = s1[0];
  
  _DMatrix<float> P(data.rows(), mixture_count);
  
  _DMatrix<T> new_covars(2,2);
  T *new_covars_ptr = new_covars[0];

  int q=0;
  for(q=0; q<max_iters; q++)
    {
      // cout << "ITERATION " << q << endl;
      this_like = 0;

      // dp->begin(4);
      for(int j=0; j<mixture_count; j++)
	{
	  std::vector<T>  likelihoods;

	  // cout<< "mixture: " << j << endl << "+++++++++++ cov" << endl;
	  // cout<< models[j].covariance() << endl << "+++++++ inv cov " << models[j].inverse_covariance() << " +++++ log c " << models[j].log_constant() << " +++++ det " << models[j].determinant_of_covariance() << " ++++++ mean " << models[j].mean() << endl;

	  models[j].get_data_likelihood(data, likelihoods);

	  typename vector<T>::const_iterator likelihood_iter;
	  float *G_ptr = G[0]+j;

	  // FIXME
	  int i=0;
	  for(likelihood_iter = likelihoods.begin(); likelihood_iter != likelihoods.end(); ++likelihood_iter, G_ptr += G.cols(), i++)
	    {
	      *G_ptr = fast_exp(*likelihood_iter);
	      // cout << "yyy " << i << " " << j << " " << *likelihood_iter << " " << fast_exp(*likelihood_iter) << " " << exp(*likelihood_iter) << endl;
	    }
	}
      // dp->end(4);

      // dp->begin(5);

      float *G_ptr = G[0];
      for(int i=0; i<data.rows(); i++)
	{
	  float l = 0;
	  for(int j=0; j<mixture_count; j++)
	    {
	      l += *(G_ptr++) * weights[j];
	    }

	  //	  cout << " === " << l << " " << log(l) << " " << fast_log(l) << endl;

	  this_like += log(l);
	}

      //      for(float iii = 1e-10; iii < 1e10; iii *= 2.0)
      // cout << " yyy " << iii << " " << log(iii) << " " << fast_log(iii) << endl;

      //      cout << "ITERATION " << q << " " << this_like << " " << fabs(this_like - last_like) << " " << delta_thresh << endl;


      // if likelihood hasn't changed much this iteration, then break
      if(fabs(this_like - last_like) < delta_thresh)
	break;



      last_like = this_like;

      // dp->end(5);
      // dp->begin(6);
      // the e-step

      for(int i=0; i<data.rows(); i++)
	{
	  float sum_j = 0;
	  G_ptr = G[i];
	  float *P_ptr = P[i];

	  for(int j=0; j<mixture_count; j++)
	    {
	      float tmp = weights[j] * G_ptr[j];
	      sum_j += tmp;
	      P_ptr[j] = tmp;
	    }
	  
	  for(int j=0; j<mixture_count; j++)
	    {
	      if(sum_j != 0)
		P_ptr[j] /= sum_j;
	      else  if(P_ptr[j])
		P_ptr[j] = 1e100;
	      else
		P_ptr[j] = -1e100;
	      
	      // cout<< "xxx " << i << " " << j << " " << P_ptr[j] << " " << weights[j] << " " << G_ptr[j] << " " << weights[j] * G_ptr[j] << endl;
	    }
	}
      // dp->end(6);
      // dp->begin(7);
      // the m-step

      for(int j=0; j<mixture_count; j++)
	{
	  s1=0;
	  float s=0;
	  T *data_cp = data[0];
	  float *P_ptr = P[0]+j;
	  for(int i=0; i<data.rows(); i++, P_ptr += P.cols())
	    {
	      float p = *P_ptr;
	      for(int l=0; l<data.cols(); l++)
		s1_ptr[l] += p * *(data_cp++);

	      s += p;
	    }
	  
	  // cout << "s = " << s << endl;

	  weights[j] = s / data.rows();

	  // not really clear what to do when s is 0. It should never happen, but
	  // presumably could. 
	  models[j].set_mean(s1.transpose() / s);

	  if(isnan(models[j].mean()[0][0]))
	     models[j].set_mean(_DMatrix<T>(2,1,0.0));
		
	  // cout<< "MEAN " << j << " " << models[j] << " " << s << endl;
		 

	  new_covars = 0;

	  data_cp = data[0];

	  float x_mean = (models[j].mean())[0][0];
	  float y_mean = (models[j].mean())[0][1];
	  P_ptr = P[0]+j;
	  for(int i=0; i<data.rows(); i++, P_ptr += P.cols())
	    {
	      float x_diff = *(data_cp++) - x_mean;
	      float y_diff = *(data_cp++) - y_mean;

	      float p = *P_ptr;

	      new_covars_ptr[0] += p * x_diff * x_diff;
	      new_covars_ptr[1] += p * x_diff * y_diff;
	      new_covars_ptr[3] += p * y_diff * y_diff;
	    }
	  
	  new_covars[1][0] = new_covars[0][1];
	  models[j].set_covariance(new_covars / s);


	  // a negative determinant here is bad, since the covariance matrix should be
	  // positive definite.
	  //
	  // We want to make sure that the smallest eigenvalue is not too small (e.g. < k).
	  // In this context the smallest eigenvalue corresponds to the covariance of
	  // the dimension with the lowest covariance when rotated to be an axis-oriented
	  // gaussian. So k=1 is probably good.
	  
	  
	  if(models[j].is_covariance_bad() || isnan(models[j].log_constant()))
	    {
	      _DMatrix<T> cov = models[j].covariance();
	      cov[0][1] = cov[1][0] = 0;
	      models[j].set_covariance(cov);
	      
	      if(models[j].is_covariance_bad() || isnan(models[j].log_constant()))
		models[j].set_covariance(_DMatrix<T>(2,2,_DMatrix<T>::identity) * 10.0);
	    }
	  


	}
      // dp->end(7);
    }

  // cout << "FINAL MODEL: " << endl;
  //  for(int i=0; i<mixture_count; i++)
  //    cout << i << " " << models[i] << " peak " << exp(models[i].log_constant()) * weights[i] << " weight " << weights[i] << endl;

  _model_likelihood = this_like;
}

template<class T>
T DGaussianMixtureModel<T>::get_data_likelihood(const _DMatrix<T> &data)
{
  std::vector<T> tmp;
  return get_data_likelihood(data, tmp);
}

template<class T>
T DGaussianMixtureModel<T>::get_data_likelihood(const _DMatrix<T> &data, std::vector<T> &likelihoods)
{
  //  std::cout << "--> in gdl data " << data << " weights " << weights[0] << " " << weights[1] << " " << weights[2] << std::endl;
  //  std::cout << models[0] << " " << models[1] << " " << models[2] << std::endl;

  likelihoods = vector<T>(data.rows(), 0);

  for(int i=0; i<models.size(); i++)
    {
      vector<T> likes;
      models[i].get_data_likelihood(data, likes);
      //      cout << " l of model " << i << " " << likes[0] << endl;

      for(int j=0; j<data.rows(); j++)
	likelihoods[j] += weights[i] * exp(likes[j]);
    }

  for(int j=0; j<data.rows(); j++)
    likelihoods[j] = log(likelihoods[j]);

  return accumulate(likelihoods.begin(), likelihoods.end(), T(0.0));
}


#define DECLARE(x) \
  template class DGaussianMixtureModel<x>; 

DECLARE(double)
DECLARE(float)
