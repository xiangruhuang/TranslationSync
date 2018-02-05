#include <DCluster.h>
#include <set>

using namespace std;

#define SWAP(a, b) { int tmp = b; b = a; a = tmp; }


template<class T>
_DMatrix<int> KMeansCluster<T>::do_clustering(const _DMatrix<T> &in_matrix)
{
  _DMatrix<int> min_assignments;
  _DMatrix<T> min_centers;
  double min_error;
  int min_cluster_count;

  min_assignments = do_clustering_singlereplicate(in_matrix);
  min_centers = get_cluster_centers();
  min_error = get_total_error();
  min_cluster_count = get_actual_cluster_count();

  for(int i=1; i<replicate_count; i++)
    {
      do_clustering_singlereplicate(in_matrix);

      if(get_total_error() < min_error)
	min_error = get_total_error(), min_centers = get_cluster_centers(), 
	  min_assignments = get_assignments(), min_cluster_count = get_actual_cluster_count();
    }

  assignments = min_assignments;
  cluster_centers = min_centers;
  total_error = min_error;
  actual_cluster_count = min_cluster_count;

  return assignments;
}

template<class T>
bool KMeansCluster<T>::find_assignments(int beg, int end, const _DMatrix<T> &in_matrix, int *assignments_ptr, double &err, const _DMatrix<T> &cluster_centers, int actual_cluster_count)
{

  err = 0;
  bool done = true;
  int dim_count = in_matrix.cols();
  //  cerr << "--> " << beg << " " << end << " " << actual_cluster_count << " " << dim_count <<  endl;
  for(int i=beg; i <= end; i++)
    {
      //      cerr << "== " << i << " " << actual_cluster_count << " " << dim_count << endl;
      const T *in_matrix_cp = in_matrix[i];
      T _min = 1e20;
      int _min_i = 0;
      if(!(i % 10000))
	cerr << ".";
	      
      const T *cluster_centers_cp = cluster_centers.data_ptr();
      for(int j=0; j < actual_cluster_count; j++)
	{
	  float dist = 0;
	  for(int k=0; k<dim_count; k++)
	    {
	      float a = in_matrix_cp[k] - *(cluster_centers_cp++);
		      
	      dist += a * a;
	    }
		  
	  if(dist < _min)
	    {
	      _min = dist;
	      _min_i = j;
	    }
	}
	  
      // a point being reassigned means that convergence
      // hasn't yet been reached.
      if(assignments_ptr[i] != _min_i)
	{
	  assignments_ptr[i] = _min_i;
	  done = false;
	}
	      
      err += double((_min));
    }

  return done;
}

template<class T>
void *KMeansCluster<T>::cluster_thread(void *_p)
{
  ClusterThreadParams<T> *p = (ClusterThreadParams<T> *)_p;
  double err;
  //  cerr << "CT " << p->actual_cluster_count << " " << p->in_matrix->cols() << endl;
  bool done = find_assignments(p->beg, p->end, *(p->in_matrix), p->assignments_ptr, err, *p->cluster_centers, p->actual_cluster_count);
  ClusterThreadResult *res = new ClusterThreadResult(done, err);
  pthread_exit((void *)res);
  return res;
}


template<class T>
bool KMeansCluster<T>::run_cluster_thread(int thread_count, int beg, int end, const _DMatrix<T> &in_matrix, int *assignments, double &total_error)
{
  if(thread_count == 0)
    return true;
  
  pthread_t p;
  ClusterThreadParams<T> params(beg, (end-beg)/thread_count + beg, &in_matrix, &cluster_centers, assignments, actual_cluster_count);

  //  cerr << "CT33 " << params.actual_cluster_count << endl;
  pthread_create(&p, NULL, &cluster_thread, (void *) &params);
  
  double err2=0;
  bool done2 = run_cluster_thread(thread_count-1, params.end+1, end, in_matrix, assignments, err2);
  
  ClusterThreadResult *result;
  pthread_join(p, (void **)&result);
  //  cerr << "joined " << (long int) result << endl;
  //  cerr << "joined2 " <<  result->error<< endl;
  
  total_error = err2 + result->error;
  bool done = result->done;
  delete result;

  //  cerr << "returning" << endl;
  return done && done2;
}



// columns are features, rows are observations
//
// If a cluster is empty, then nothing in assignments will point to it and 
// the corresponding entries of cluster_centers are nan
//
template<class T>
_DMatrix<int> KMeansCluster<T>::do_clustering_singlereplicate(const _DMatrix<T> &in_matrix)
{
  cerr << " clustering: setting up... " << endl;
  int dim_count = in_matrix.cols();
  int pt_count = in_matrix.rows();
  actual_cluster_count = min(requested_cluster_count, pt_count);

  cerr << " clustering:  " << requested_cluster_count << " " << dim_count << endl;
  cluster_centers = _DMatrix<T>(requested_cluster_count, dim_count);
  set<int> used_pts;

  cerr << " clustering: assigning initial clusters..." << endl;
  // randomly select points for initial cluster centers
  // (also make sure no 2 clusters have the same center)
  // FIXME: this only checks that clusters have different
  // data index ids, not actually that the pts are different
  //
  for(int i=0; i<actual_cluster_count; i++)
    {
      int pt;

      do {
	pt = random() % pt_count;
      } while (used_pts.find(pt) != used_pts.end());

      cluster_centers.set_row(i, in_matrix.extract_row(pt));
      used_pts.insert(pt);
    }

  cerr << " clustering..." << endl;
  bool done=false;
  assignments = _DMatrix<int>(1, pt_count);
  assignments = -1;

  int *assignments_ptr = assignments.data_ptr();
  int iter_count = 0;

  cerr << " starting..." << endl;
  while(!done && iter_count < max_iters)
    {
      // first: assign points to closest clusters
      total_error=0;

      cerr << "  pt count " << pt_count;

#ifdef HAVE_PTHREADS
      if(thread_count > 1)
        done = run_cluster_thread(thread_count, 0, pt_count-1, in_matrix, assignments_ptr, total_error);
      else
#endif
	done = find_assignments(0, pt_count-1, in_matrix, assignments_ptr, total_error, cluster_centers, actual_cluster_count);

      cerr << endl;
      cerr << "  total error " << total_error << endl;

      // second: recompute cluster centers
      cluster_centers = 0;
      
      int counts[actual_cluster_count];
      memset(counts, 0, sizeof(int) * actual_cluster_count);
      T *in_cp = in_matrix.data_ptr();
      
      for(int i=0; i < pt_count; i++)
	{
	  int which_cluster = assignments_ptr[i];
	  counts[which_cluster]++;

	  T *center_cp = cluster_centers[which_cluster];

	  for(int j=0; j < dim_count; j++)
	    center_cp[j] += *(in_cp++);
	}

      T *cluster_centers_cp = cluster_centers.data_ptr();
      for(int i=0; i < actual_cluster_count; i++)
	for(int j=0; j < dim_count; j++)
	  *(cluster_centers_cp++) /= counts[i];

      cerr << "  second " << endl;

      // check for the case of an empty cluster. Handle this by
      // moving cluster to the end and decreasing actual_cluster_count.
      for(int i=0; i<actual_cluster_count; i++)
	{
	  if(counts[i] == 0)
	    {
	      assignments.search_and_replace_ip(actual_cluster_count-1, i);
	      cluster_centers.swap_rows(actual_cluster_count-1, i);
	      swap(counts[actual_cluster_count-1], counts[i]);
	      actual_cluster_count--; i--;
	    }
	}

      cerr << "iteration " << iter_count << endl;
      iter_count++;
    }

  // mark empty cluster centers as "nan"
  _DMatrix<T> no_cluster(1, dim_count, T(nan("")));
  for(int i=actual_cluster_count; i < requested_cluster_count; i++)
    cluster_centers.set_row(i, no_cluster);

  return assignments;
}

template<class T>
pair<int, T> KMeansCluster<T>::get_assignment(const vector<T> &vec)
{
  assert(vec.size() == cluster_centers.cols());
  double _min = 1e100;
  int _min_i;

  const T *cluster_centers_cp = cluster_centers.data_ptr();
  for(int j=0; j < actual_cluster_count; j++)
    {
      T dist = 0;
      for(int k=0; k<cluster_centers.cols(); k++)
	{
	  T a = vec[k] - *(cluster_centers_cp++);
	  
	  dist += a * a;
	}

      if(dist < _min)
	{
	  _min = dist;
	  _min_i = j;
	}
    }

  return make_pair(_min_i, _min);
}
//#define DECLARE(x) \
//  template class KMeansCluster<x>; 


//DECLARE(double)
//DECLARE(float)

template class KMeansCluster<float>;
template class KMeansCluster<double>;
