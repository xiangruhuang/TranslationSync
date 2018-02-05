#include "DImage.h"
#include <DMultiDMatrix.h>
#include <vector>
#include <istream>
#include <fstream>
#include <ext/hash_map>
#include <set>
using namespace std;
#include <io.h>
#include <algorithm>
#include <numeric>
#include <map>
#include <math.h>
#include "bp-common.h"
#include <stack>

using namespace __gnu_cxx;
using namespace std;

// input files
char *pairwise_file_g = 0; //, *unary_file_g = 0;
// # of BP iterations to run
int max_ii_g = 100;
// 
double pair_trunc_g = 1.0;
// older file format doesn't have confidence field
bool no_confidence_in_file_g = false, two_confidence_in_file_g = false;

// ground truth file
char *gt_g = 0;

// for debugging...
double stretch_g = 4;
int bufsz_g = 4000;
int verbose_g=0; 
int node_g = -1;
char *rescore_g = 0;

//
bool linear_cost_g = false;
bool discrete_g = true;
// print to stdout info about geotag-based edge priors then exit
bool write_geo_priors_g = true;

bool scale_geotag_unary_g = false;
double scaled_geotag_weight_g = 0.22;
double unary_geotag_weight_g = 0.05;
double unary_geotag_trunc_g = 0.7;
double geotag_stdev_g = 0; // 53 
int num_unary_samples_g = 1;
bool geotag_pan_g = false;

/// Rotations

//int TWIST_BIN_COUNT = 5;
//double TWIST_BIN_SPACING = 0.15;
//double TWIST_BIN_MIN = -0.3;
int TWIST_BIN_COUNT = 1;
double TWIST_BIN_SPACING = 0;
double TWIST_BIN_MIN = 0;
int GRID_COUNT = 21; // 21; // 11; // 21;
int LABEL_COUNT = GRID_COUNT * GRID_COUNT * GRID_COUNT * TWIST_BIN_COUNT;
double DIST_GRAN = 10.0; // 10.0; // 5.0; // 10.0;
double DIST_STEP = 0.1;
int ref_node_g = -1;
int thread_count_g = 1;
//
bool random_ref_g = false;
bool no_remove_cut_g = false;
bool ccref_highest_degree_g = false;
bool no_labels_g = false;
bool write_edges_g = false;
bool prepropagate_g = true;
bool param_est_g = false;

char *vanish_g = 0;
char *perfect_pairs_g = 0;
char *geoplanar_g = 0;

double unary_scale_g = 1.0, unary_trunc_g = 1.0;

DistributionMap unary_dists;

typedef _DMultiDMatrix<double> Msgtmp_type ;
//typedef hash_map<int, hash_map<int, PairwiseDiff > > Pairs;

template<class T>
inline double vector_norm(const _DMatrix<T> &matrix)
{
  return sqrt(pointwise_multiply(matrix, matrix).sum()); 
}

inline int unpack_tw(int packed) { return packed / (GRID_COUNT * GRID_COUNT * GRID_COUNT); }
inline int unpack_t0(int packed) { return (packed / (GRID_COUNT * GRID_COUNT)) % GRID_COUNT; }
inline int unpack_t1(int packed) { return (packed / GRID_COUNT) % GRID_COUNT; }
inline int unpack_t2(int packed) { return packed % GRID_COUNT; }

inline int from_dist(double s) { return (int)round(s * double(DIST_GRAN) + GRID_COUNT/2); }
inline double to_dist(double s) { return (s - GRID_COUNT/2) / double(DIST_GRAN); }
inline int pack(int t0, int t1, int t2, int tz=0) { return tz * GRID_COUNT * GRID_COUNT * GRID_COUNT + t0 * GRID_COUNT * GRID_COUNT + t1 * GRID_COUNT + t2; }

inline int from_twist(double s) { int s2 = (int) round((s-TWIST_BIN_MIN) / TWIST_BIN_SPACING); if(s2 < 0) s2 = 0; if(s2 >= TWIST_BIN_COUNT) s2 = TWIST_BIN_COUNT-1; return s2; }
inline double to_twist(double s) { return TWIST_BIN_MIN + s * TWIST_BIN_SPACING; }

template<class T>
inline T sqr(T a) { return a * a; }

// FIXME: this function is called by inner_loop. Make it faster! 
// it would be good to get rid of the conditional here. That's possible with integers using bit twiddling,
// but I'm not sure how to do it with floats. We could possibly move to fixed-precision arithmetic.
inline double dist(double s)
{
  /*
    if(linear_cost_g)
    {
    if(s < -pair_trunc_g) return pair_trunc_g;
    else if(s > pair_trunc_g) return pair_trunc_g;
    else if(s < 0) return -s;
    else return s;
    }
  */

  //  if(s < -pair_trunc_g) return pair_trunc_g*pair_trunc_g;
  //  else if(s > pair_trunc_g) return pair_trunc_g*pair_trunc_g;
  //  else return s*s;
  return min(s*s, pair_trunc_g * pair_trunc_g);
}

/*
  double distance2_correct(int my_state, int other_state, const PairwiseDiff &diff,  const vector< _DMatrix<TYPE> > &rotation_matrices)
  {
  _DMatrix<TYPE> diff2;
  change_type(-diff.rot.extract_col(2), diff2);

  _DMatrix<TYPE> predicted_dir = rotation_matrices[my_state] * diff2;
  _DMatrix<TYPE> other = -rotation_matrices[other_state].extract_col(2);

  return(dist(predicted_dir[0][0]-other[0][0]) + dist(predicted_dir[0][1]-other[0][1]) + dist(predicted_dir[0][2] - other[0][2]));
  }
*/

double distance2_correct(const PairwiseDiff &diff, DMatrix &mystate_rot, const DMatrix &other)
{
  DMatrix predicted_dir = mystate_rot * -diff.rot.extract_col(2);
  //  double predicted_dir[3];
  //  for(int i=0; i<3; i++)
  //    for(int k=0; k<3; k++)
  //      predicted_dir[i] = mystate_rot[i][k] * -diff.rot[2+k*3]; 

  return(dist(predicted_dir[0][0]-other[0][0]) + dist(predicted_dir[0][1]-other[0][1]) + dist(predicted_dir[0][2] - other[0][2]));
  //  double *other_ptr = other[0];
  //  return dist(predicted_dir[0] - other_ptr[0]) + dist(predicted_dir[1] - other_ptr[1]) + dist(predicted_dir[2] - other_ptr[2]);
}


int inside_sphere(double x, double y, double z, double r = 1.0)
{
  return x*x + y*y + z*z <= r*r;
}

bool is_state_on_sphere(int state)
{
  static bool first = true;
  static bool *intersect = 0;

  if(first)
    {
      intersect = new bool[LABEL_COUNT];
      for(int s=0; s<LABEL_COUNT; s++)
        {
          double x = to_dist(unpack_t0(s));
          double y = to_dist(unpack_t1(s));
          double z = to_dist(unpack_t2(s));

          // there are eight corners of this cube. if at least 1 is inside the sphere, and
          // at least 1 is outside the sphere, then the surface of the sphere intersects this cube.
          double half_bin = 1.0/(2.0 * DIST_GRAN);
          int in_count=0;
          for(int xx = -1; xx <= 1; xx+=2)
            for(int yy = -1; yy <= 1; yy+=2)
              for(int zz = -1; zz <= 1; zz+=2)
                in_count += inside_sphere(x + xx*half_bin, y + yy*half_bin, z + zz*half_bin);

          intersect[s] = in_count > 0 && in_count < 8;
        }

      int c=0;
      for(int i=0; i<LABEL_COUNT; i++)
        if(intersect[i]) c++;
      first = false;
      cout << "(" << c << " " << LABEL_COUNT << " " << DIST_GRAN << ") " << endl;
    }

  return intersect[state];
}

// *all* geotags projected onto a plane
void read_2d_locations_file(const string &fname, hash_map<int, DMatrix> &locations) {
  ifstream f(fname.c_str());
  if (!f.good()) {
    cerr << "cannot find " << fname << endl;
  }
  assert(f.good());
  while (f.good()) {
    int id;
    DMatrix xyz(3,1);
    f >> id >> xyz[0][0] >> xyz[2][0];
    xyz[1][0] = 0;
    locations[id] = xyz;
    //cerr << unary[id].position << endl;
  }
  f.close();
}

void parse_opts(int argc, char *argv[])
{
  int ii=1;
  for(ii=1; ii<argc; ii++)
    {
      if(!strcmp(argv[ii], "--iters"))
        max_ii_g = atoi(argv[++ii]);
      else if(!strcmp(argv[ii], "--gt"))
        gt_g = argv[++ii];
      else if(!strcmp(argv[ii], "--threads"))
        thread_count_g = atoi(argv[++ii]);
      else if(!strcmp(argv[ii], "--verbose"))
        verbose_g++;
      else if(!strcmp(argv[ii], "--writeedges"))
        write_edges_g = true;
      else if(!strcmp(argv[ii], "--twist"))
        {
          TWIST_BIN_COUNT=5;
          TWIST_BIN_SPACING = 0.1;
          TWIST_BIN_MIN = -0.2;
        }
      else if(!strcmp(argv[ii], "--noconf"))
        no_confidence_in_file_g = true;
      else if(!strcmp(argv[ii], "--hideg"))
        ccref_highest_degree_g = true;
      else if(!strcmp(argv[ii], "--nocut"))
        no_remove_cut_g = true;
      else if(!strcmp(argv[ii], "--randomref"))
        random_ref_g = true;
      else if(!strcmp(argv[ii], "--twoconf"))
        two_confidence_in_file_g = true;
      else if(!strcmp(argv[ii], "--rescore"))
        rescore_g = argv[++ii];
      else if(!strcmp(argv[ii], "--pairtrunc"))
        pair_trunc_g = atof(argv[++ii]);
      else if(!strcmp(argv[ii], "--linear"))
        linear_cost_g = true;
      else if(!strcmp(argv[ii], "--refnode"))
        ref_node_g = atoi(argv[++ii]);
      else if(!strcmp(argv[ii], "--nolabels"))
        no_labels_g = true;
      else if(!strcmp(argv[ii], "--nopreprop"))
        prepropagate_g = false;
      else if(!strcmp(argv[ii], "--makeperfectpairs"))
        perfect_pairs_g = argv[++ii];
      else if(!strcmp(argv[ii], "--vanish"))
        vanish_g = argv[++ii];
      else if(!strcmp(argv[ii], "--paramest"))
        param_est_g = true;
      else if(!strcmp(argv[ii], "--labelspace"))
        {
          GRID_COUNT =  atoi(argv[++ii]); 
          LABEL_COUNT = GRID_COUNT * GRID_COUNT * GRID_COUNT * TWIST_BIN_COUNT;
          DIST_GRAN = atof(argv[++ii]); 
        }
      else if (!strcmp(argv[ii], "--geoplanar"))
        geoplanar_g = argv[++ii];
      else if (!strcmp(argv[ii], "--unarygeotagweight"))
        unary_geotag_weight_g = atof(argv[++ii]);
      else if (!strcmp(argv[ii], "--unarygeotagtrunc"))
        unary_geotag_trunc_g = atof(argv[++ii]);
      else if (!strcmp(argv[ii], "--geotagstdev"))
        geotag_stdev_g = atof(argv[++ii]);
      else if (!strcmp(argv[ii], "--numunarysamples"))
        num_unary_samples_g = atoi(argv[++ii]);
      else if (!strcmp(argv[ii], "--scalegeotagunary"))
        scale_geotag_unary_g = atoi(argv[++ii]);
      else if (!strcmp(argv[ii], "--scaledgeotagweight"))
        scaled_geotag_weight_g = atof(argv[++ii]);
      else if (!strcmp(argv[ii], "--geotagpan"))
        geotag_pan_g = atoi(argv[++ii]);
      else
        break;
    }

  // linear pairwise costs currently disabled
  assert(!linear_cost_g);

  DIST_STEP = 1.0 / DIST_GRAN;

  if(ii < argc)
    pairwise_file_g = argv[ii++];

  cerr << pairwise_file_g  << endl;

  if(argc != ii)
    {
      cerr << "bad commandline options." << endl;
      exit(1);
    }
}


inline TYPE inner_loop(double *scratch, double mu)
{
  double _min=INFINITY; 
  
  for(int z2=0; z2<GRID_COUNT; z2++) //"mine"
    {
      double this_dist = scratch[z2] + dist(to_dist(z2) - mu);
      
      if(this_dist < _min) 
        _min = this_dist;
    }
 
  return _min;
}  

// this version gives slightly different results than the above function (due to numerical issues),
//  but results should be of similar quality and this version is twice as fast
//
// FIXME: replace call to dist() with a lookup table?
inline TYPE inner_loop_opt(double *scratch, double mu)
{
  double _min=INFINITY; 
  
  double inc = 1.0/DIST_GRAN; 
  double z3= (-GRID_COUNT / 2) / double(DIST_GRAN) - mu;
  for(int z2=0; z2<GRID_COUNT; z2++, z3 += inc) //"mine"
    {
      double this_dist = scratch[z2] + dist(z3);
      
      if(this_dist < _min) 
        _min = this_dist;
    }

  /*
  // unrolled version -- doesn't seem to help much. The conditionals in dist() are probably the bottleneck.

  double td[12];
  switch(GRID_COUNT)
  {
  case 11: td[10] = scratch[10] + dist(z3+inc*10); 
  case 10: td[9] = scratch[9] + dist(z3+inc*9); if(td[10] < td[9]) td[9] = td[10];
  case 9:  td[8] = scratch[8] + dist(z3+inc*8); if(td[9] < td[8]) td[8] = td[9];
  case 8:  td[7] = scratch[7] + dist(z3+inc*7); if(td[8] < td[7]) td[7] = td[8];
  case 7:  td[6] = scratch[6] + dist(z3+inc*6); if(td[7] < td[6]) td[6] = td[7];
  case 6:  td[5] = scratch[5] + dist(z3+inc*5); if(td[6] < td[5]) td[5] = td[6];
  case 5:  td[4] = scratch[4] + dist(z3+inc*4); if(td[5] < td[4]) td[4] = td[5];
  case 4:  td[3] = scratch[3] + dist(z3+inc*3); if(td[4] < td[3]) td[3] = td[4];
  case 3:  td[2] = scratch[2] + dist(z3+inc*2); if(td[3] < td[2]) td[2] = td[3];
  case 2:  td[1] = scratch[1] + dist(z3+inc*1); if(td[2] < td[1]) td[1] = td[2];
  case 1:  td[0] = scratch[0] + dist(z3+inc*0); if(td[1] < td[0]) td[0] = td[1]; break;
  default: assert(0);
  }
  _min = td[0];*/

  return _min;
}  

inline double multiply_row_with_vector(TYPE *row, TYPE *vec)
{
  return row[0] * vec[0] + row[1] * vec[1] + row[2] * vec[2];
}

// have rotation matrix C0 estimated based on v0, so multiply C0*(P01*r).
// FIXME: there used to be a heuristic here that ignored columns that are all INF on the first dimension,
//  leading to minor (few %) speedup
//
// d1, d2, d3: some permutation of the variable names x, y, z. DT is done over dimension d3 (i.e. for each value of d1,d2)
// D1, D2, D3: same as d1, d2, d3, but with capital variable names X, Y, Z
// row1, row2, row3: some permutation of 0, 1, 2: same order as d1, d2, d3, where 0=x, 1=y, 2=z
// P01r : final column of expected rotation matrix between the two cameras.
//
#define LOOP_CORRECT(d1,d2,d3, D1, D2, D3, row1, row2, row3, P01r, __msg_tmp, rot_mats) \
  {                                                                     \
    for(int d1=0; d1<GRID_COUNT; d1++)                                  \
      for(int d2=0; d2<GRID_COUNT; d2++)                                \
        {                                                               \
          for(int d3=0; d3<GRID_COUNT; d3++)                            \
            scratch[d3] = __msg_tmp.get(x)[y][z];                       \
                                                                        \
          for(int d3=0; d3<GRID_COUNT; d3++)                            \
            {                                                           \
              double mu = multiply_row_with_vector(rot_mats[pack(x,y,z,twist)][row3], P01r); \
              __msg_tmp.get(x)[y][z] = inner_loop_opt(scratch, mu);     \
            }                                                           \
        }                                                               \
  }
// FIXME: remove twist in above call to pack to speed up.


void subtract_dc(Message &this_msg)
{
  // subtract off dc component
  double sum_min_cost = 0;
  int count=0;
  for(int other_s = 0; other_s < LABEL_COUNT; other_s++)
    if(!isinf(this_msg[other_s] ))
      sum_min_cost += this_msg[other_s], count++;
  
  if(count == 0) count=1;
  sum_min_cost /= double(count);
  for(int i=0; i< LABEL_COUNT; i++)
    this_msg[i] -= sum_min_cost;
  
  // This step unnecessary, but makes messages more readable (for debugging)
  for(int i=0; i< LABEL_COUNT; i++)
    if(!is_state_on_sphere(i))
      this_msg[i] = INFINITY;
}


void make_H(Distribution &H, int me, int cc_ref)
{
  H = unary_dists[me];
}

void read_messages(Distribution &H, int me, int him, const Pairs &deltas, const MessageSet &last_ms)
{
  //  NeighborPairs &my_neighbors = deltas[me];
  Pairs::const_iterator deltas_iter_me = deltas.find(me);
  assert(deltas_iter_me != deltas.end());
  const NeighborPairs &my_neighbors = deltas_iter_me->second;

  //  hash_map<int, Message> &my_incoming_messages = last_ms[me];
  MessageSet::const_iterator lastms_me_iter = last_ms.find(me);
  assert(lastms_me_iter != last_ms.end());
  const hash_map<int, Message> &my_incoming_messages = lastms_me_iter->second;

  for(NeighborPairs::const_iterator iter = my_neighbors.begin(); iter != my_neighbors.end(); ++iter)
    {
      if(iter->first != him) 
        {
          //	  Message &this_msg = my_incoming_messages[iter->first];
          hash_map<int, Message>::const_iterator this_msg_iter = my_incoming_messages.find(iter->first);
          assert(this_msg_iter != my_incoming_messages.end());
          const Message &this_msg = this_msg_iter->second;

          for(int my_s = 0; my_s < LABEL_COUNT; my_s++)
            H[my_s] += this_msg[my_s];
        }
    }

}

DistributionMap compute_final_distribution(MessageSet &current_ms, const vector<Edge> &edge_list, int cc_ref)
{
  DistributionMap result;

  for(int i=0; i<(int)edge_list.size(); i++)
    {
      Edge edge = edge_list[i];
      int me = edge.first, him = edge.second;
    
      if(result.find(him) == result.end())
        result[him] = vector<TYPE>(LABEL_COUNT);
    
      const Message &this_msg = current_ms[him][me];
      vector<TYPE> &this_final = result[him];

      for(int j=0; j<LABEL_COUNT; j++)
        this_final[j] += this_msg[j];
    }

  for(DistributionMap::iterator iter = result.begin(); iter != result.end(); ++iter)
    {
      Distribution H;
      make_H(H, iter->first, cc_ref);

      for(int j=0; j<LABEL_COUNT; j++)
        iter->second[j] += H[j];
    }

  return result;
}

DMatrix norm_dir_of_label(int my_state, bool include_twist)
{
  DMatrix vec(4,1);
  vec[0][0] = to_dist(unpack_t0(my_state)); vec[0][1] = to_dist(unpack_t1(my_state)); vec[0][2] = to_dist(unpack_t2(my_state)); vec[0][3] = 0;
  double norm = (sqrt(pointwise_multiply(vec, vec).sum()));
  if(norm > 0) vec /= norm;
  vec[0][3] = to_twist(unpack_tw(my_state));
  
  return vec;
}

vector<int> find_map_estimates(DistributionMap &final_distribution, int photo_count, map<int, DMatrix> &estimated_dirs, int ii)
{
  vector<int> labels(photo_count);

  for(int him=0; him<photo_count; him++)
    {
      if(!no_labels_g)
        cout << "labels " << ii << " " << him << " ";

      if(final_distribution.find(him) == final_distribution.end())
        {
          cout << -1 << " " << 0 << " " << 0 << " " << 1 << " " << 0 << " " << 0 << " " << 0 << " " << 1 << endl;
          continue;
        }
      
      double min_cost = INFINITY;
      int min_j=0;
      
      for(int j=0; j<LABEL_COUNT; j++)
        if(is_state_on_sphere(j))
          if(min_cost > final_distribution[him][j])
            min_cost = final_distribution[him][j], min_j = j;
      
      labels[him] = min_j;
      //      estimated_dirs[him] = norm_dir_of_label(min_j);
      estimated_dirs[him] = norm_dir_of_label(min_j, true);
      
      //      cout << him << ":" << "(" << min_j << ")," << to_dist(unpack_t0(min_j)) << "," << to_dist(unpack_t1(min_j)) << "," << to_dist(unpack_t2(min_j)) << "," << to_twist(unpack_tw(min_j)) << " ";
      if(!no_labels_g)
        cout << min_j << " " << to_dist(unpack_t0(min_j)) << " " << to_dist(unpack_t1(min_j)) << " " << to_dist(unpack_t2(min_j)) << " " << to_twist(unpack_tw(min_j)) << " " << estimated_dirs[him][0][0] << " " << estimated_dirs[him][0][1] << " " << estimated_dirs[him][0][2] << endl;
    }
  if(!no_labels_g)
    cout << endl;

  return labels;
}

template<class T>
inline int sign(T X)
{
  if(X >= 0) return 1;
  return -1;
}


template<class T>
_DMatrix<T> compute_zerotwist_rot(double X, double Y, double Z, double twist=0)
{
  _DMatrix<T> result(3,3);
  
  double len = sqrt(X*X + Y*Y + Z*Z);
  X /= len; Y /= len; Z /= len;
  double len2 = sqrt(X*X + Z*Z);
  double X_ = X / len2, Z_ = Z / len2;
  
  result[0][0] = -Z_;   result[0][1] = -Y * X_;           result[0][2] = -X;
  result[1][0] = 0;     result[1][1] = Z * Z_ + X * X_;   result[1][2] = -Y;
  result[2][0] = X_;    result[2][1] = -Y * Z_;           result[2][2] = -Z;
  
  return result;
}



// rotation_matrix_abstwist : produce a 3x3 rotation matrix from a viewing direction and twist angle
//  X, Y, Z is the camera viewing direction
//  abstwist represents the twist in terms of an absolute y coordinate (in abs x-y-z space)
//
#define QUANTUM 1e-10
template<class T>
_DMatrix<T> rotation_matrix_abstwist(double X, double Y, double Z, double abstwist=0)
{
  _DMatrix<T> result(3,3);

  // last column of R is just -(X,Y,Z), (normalized just in case)
  double len = sqrt(X*X + Y*Y + Z*Z);
  X /= len; Y /= len; Z /= len;

  // if X and Z are both 0, then the camera is pointing straight up or down.
  // for now, just handle this as if abstwist were 0.
  if(fabs(X) < QUANTUM && fabs(Z) < QUANTUM)
    return compute_zerotwist_rot<T>(X, Y, Z, 0);

  // first column is (A, B, C), where B=asin(twist), |(A,B,C)=1|, and (A,B,C) . (X,Y,Z) = 0.
  // solve for C with quadratic equation:
  // (Z^2/X^2 + 1) * C^2 + 2*Z*B*Y/X^2*C + (B^2*Y^2/X^2 - 1 + B^2)
  double B = abstwist;
  double A, C, A1, C1, A2, C2;
  if(fabs(X) < QUANTUM)
    {
      C2 = C1 = B*Y/Z;
      A1 = sqrt(1-B*B*(1+Y*Y/(Z*Z)));
      A2 = -sqrt(1-B*B*(1+Y*Y/(Z*Z)));
    }
  else
    {
      double X_sqr_inv = 1.0/(X*X);
      double _a = (Z*Z)*X_sqr_inv + 1, _b = 2*Z*B*Y*X_sqr_inv, _c = B*B*Y*Y*X_sqr_inv - 1 + B*B;
      double _sqrt = sqrt(_b*_b - 4 * _a * _c);
      C1 = (-_b + _sqrt)/(2*_a), C2 = (-_b - _sqrt)/(2*_a);
      A1 = (-B*Y-C1*Z) / X, A2 = (-B*Y-C2*Z) / X;
    }

  C = C1, A = A1;
  
  // choose solution that is "90 left" of viewing direction (assumes twist angle 
  //  relatively small, definitely smaller than 90 degrees)
  double signmag = fabs(A2)+fabs(Z) - ( fabs(C2)+fabs(X) );
  if((sign(A2) != sign(Z) && signmag >= 0) || (sign(C2) == sign(X) && signmag < 0 ))
    C = C2, A = A2;

  // second column is -(A,B,C) x (X, Y, Z)
  result[0][0] = A; 	result[0][1] = B*Z - C*Y;     result[0][2] = -X;
  result[1][0] = B; 	result[1][1] = C*X - A*Z;     result[1][2] = -Y;
  result[2][0] = C; 	result[2][1] = A*Y - B*X;     result[2][2] = -Z;

  //  cout << "-----" << Endlx << result << endl << compute_zerotwist_rot<T>(X, Y, Z)<< endl;

  //  _DMatrix<T> m = fabs(result - compute_zerotwist_rot<T>(X, Y, Z));
  //  cout << "ERROR " << pointwise_multiply(m,m).sum() << endl;

  return result;
}

// in this version, twist is an angle about the viewing direction (X,Y,Z)
template<class T>
_DMatrix<T> rotation_matrix(double X, double Y, double Z, double twist=0)
{
  double abs_twist = sin(twist) * cos(Y);
  return rotation_matrix_abstwist<T>(X, Y, Z, abs_twist);
}

template<class T>
double twist_of_rotation_matrix(const _DMatrix<T> &rot)
{
  return acos(std::min((rot[0][0] * rot[2][2] - rot[2][0] * rot[0][2]) / sqrt(1 - rot[1][2]*rot[1][2]), 1.0));
}

template<class T>
double tilt_of_rotation_matrix(const _DMatrix<T> &rot)
{
  return asin(-rot[1][2]);
}

double compute_labeling_cost(map<int, DMatrix> &estimated, Pairs &deltas, vector<Edge> &edge_list, int ii)
{
  double cost = 0;
  // pairwise costs
  for(int i=0; i<(int)edge_list.size(); i++)
    {
      Edge edge = edge_list[i];
      PairwiseDiff this_delta_rot = deltas[edge.first][edge.second];
      int me = edge.first, him = edge.second;

      if(estimated.find(me) == estimated.end() || estimated.find(him) == estimated.end())
        {
          ostringstream ofs;
          ofs << "missing estimate for camera " << me << " or " << him;
          throw(ofs.str());
        }

      DMatrix mystate_rots = rotation_matrix<double>(estimated[me][0][0], estimated[me][0][1], estimated[me][0][2], estimated[me][0][3]);
      
      cost += distance2_correct(this_delta_rot, mystate_rots, estimated[him]);
    }

  // incorporate unary costs
  for(map<int, DMatrix>::iterator iter = estimated.begin(); iter != estimated.end(); ++iter)
    cost += unary_dists[iter->first][pack(from_dist(iter->second[0][0]), from_dist(iter->second[0][1]), from_dist(iter->second[0][2]))];
  
  return cost;
}


DMatrix compute_ls_rotation(const DMatrix &est_matrix, const DMatrix &gt_matrix)
{
  DMatrix corr(3,3);
  
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      {
        DMatrix X = est_matrix.extract_col(i) - est_matrix.extract_col(i).mean();
        DMatrix Y = gt_matrix.extract_col(j) - gt_matrix.extract_col(j).mean();
        corr[i][j] = pointwise_multiply(X,Y).mean() / (sqrt(X.covariance()[0][0]) * sqrt(Y.covariance()[0][0]));
      }
  
  cerr << "CORR = " << corr << endl;

  DMatrix u,s,v;
  corr.svd(u,s,v);
  corr = u * v.transpose();
  
  return corr;
}


double compute_gt_error(double &err, hash_map<int, PairwiseDiff> &gt, const set<int> &reachable, map<int, DMatrix> &estimated_dirs, int photo_count, int ii, double &twist_error)
{
  double err2=0;
  err = 0;
  twist_error = 0;
  
  DMatrix gt_matrix(0,3), est_matrix(0,3); 
  for(int i=0; i<photo_count; i++)
    {
      if(reachable.find(i) == reachable.end() || gt.find(i) == gt.end())
        continue;
 
      //      DMatrix last_col(1,3);
      //      last_col[0][0] = -gt[i].rot[2];
      //      last_col[0][1] = -gt[i].rot[5];
      //      last_col[0][2] = -gt[i].rot[8];
      gt_matrix = vert_concat(gt_matrix, -gt[i].rot.extract_col(2).transpose());
      est_matrix = vert_concat(est_matrix, estimated_dirs[i].transpose().extract(DRect(0,0,0,2)));
    }	

  cout << "GT " << est_matrix << endl << gt_matrix << endl;
  
  // find best (least squared error) rotation matrix between ground truth and estimates
  DMatrix gt_rot = compute_ls_rotation(est_matrix, gt_matrix);

  int ct=0;
  for(int i=0; i<photo_count; i++)
    {
      if(reachable.find(i) == reachable.end() || gt.find(i) == gt.end())
        continue;
      
      DMatrix gt_vec = -gt[i].rot.extract_col(2);
      err += vector_norm(gt_vec - estimated_dirs[i].extract(DRect(0,0,2,0)));
      
      DMatrix est_rot_mat = rotation_matrix<double>(estimated_dirs[i][0][0], estimated_dirs[i][0][1], estimated_dirs[i][0][2], estimated_dirs[i][0][3]);
      cout << "ZZZZ " << gt_rot.transpose() << " " << est_rot_mat << endl;
      DMatrix vec2 = gt_rot.transpose() * est_rot_mat;
      DMatrix vec = gt_rot.transpose() * estimated_dirs[i].extract(DRect(0,0,2,0));
      double V = vector_norm(vec - gt_vec);

      double V_ang = acos((gt_vec[0][0] *  vec[0][0] +  gt_vec[0][1] * vec[0][1] +   gt_vec[0][2] * vec[0][2]) /   (vector_norm(vec) *  vector_norm(gt_vec)));
      double V_ang2= acos((gt_vec[0][0] * -vec2[0][2] + gt_vec[0][1] * -vec2[1][2] + gt_vec[0][2] * -vec2[2][2]) / (vector_norm(vec2.extract_col(2)) * vector_norm(gt_vec)));
      double V_tw_ang = fabs(twist_of_rotation_matrix(gt[i].rot) - estimated_dirs[i][0][3]);

      twist_error += V_tw_ang;
      err2 += V;
      ct++;
      cout << "gtcost " << ii << " " << "cam " << i << " " << gt_vec[0][0] <<  "," << gt_vec[0][1] << "," << gt_vec[0][2] << "," << gt[i].rot[1][0]  << "(" << gt[i].rot[1][0] << " tw " << twist_of_rotation_matrix(gt[i].rot) << ") " << -vec2[0][2] << "," << -vec2[1][2] << "," << -vec2[2][2] << "," << vec2[1][0] << "(" << estimated_dirs[i][0][3] << "," << est_rot_mat[1][0] << ") " << estimated_dirs[i][0][3] << " " << V << " " << V_ang << " " << V_tw_ang << " " << V_ang2 << endl;
    }

  cerr << "count = " << ct << endl;
  err /= double(ct);
  twist_error /= double(ct);
  return err2 / double(ct);
}

void copy_message(Msgtmp_type &msg_tmp, Distribution &H, int H_offset)
{
  // create msg_tmp, a 3-d matrix version of H
  for(int i=0, cp=H_offset; i<GRID_COUNT; i++)
    for(int j=0; j<GRID_COUNT; j++)
      for(int k=0; k<GRID_COUNT; k++, cp++)
        msg_tmp.get(i)[j][k] = H[cp];
}

void copy_message(Message &this_msg, int offset, Msgtmp_type &msg_tmp)
{
  for(int i=0, cp=offset; i<GRID_COUNT; i++)
    for(int j=0; j<GRID_COUNT; j++)
      for(int k=0; k<GRID_COUNT; k++, cp++)
        this_msg[cp] = msg_tmp.get(i)[j][k];
}

inline double twist_dist(double t, double t2)
{
  return (t2-t)*(t2-t) * TWIST_BIN_SPACING*TWIST_BIN_SPACING;
}

static pthread_mutex_t mymutex = PTHREAD_MUTEX_INITIALIZER;

void compute_message(int me, int him, int cc_ref, const Pairs &deltas, const MessageSet &last_ms, MessageSet &current_ms, const vector< _DMatrix<TYPE> > &rotation_matrices)
{
  // FIXME: this constructor probably takes a non-trivial amount of time. 
  vector<_DMultiDMatrix<double> > msg_tmp(TWIST_BIN_COUNT, _DMultiDMatrix<double>(3, GRID_COUNT*TWIST_BIN_COUNT, GRID_COUNT, GRID_COUNT));
  double scratch[GRID_COUNT];
  Distribution H(LABEL_COUNT);
  // make buffer of unary potentials
  make_H(H, me, cc_ref);
  
  // incorporate opinions of neighboring nodes
  read_messages(H, me, him, deltas, last_ms);
  
  // get expected (predicted) transformation from his viewing direction to mine
  Pairs::const_iterator deltas_iter_me = deltas.find(me);
  assert(deltas_iter_me != deltas.end());
  const NeighborPairs &my_neighbors = deltas_iter_me->second;

  NeighborPairs::const_iterator deltas_iter_mehim = my_neighbors.find(him);
  assert(deltas_iter_mehim != my_neighbors.end());
  const PairwiseDiff &mehim_pd = deltas_iter_mehim->second;
  _DMatrix<double> v_rot = -mehim_pd.rot.transpose().extract_col(2); 
  //  _DMatrix<double> v_rot = -deltas[me][him].rot.transpose().extract_col(2); 
  
  TYPE v_rot_row[3];
  v_rot_row[0] = v_rot[0][0]; v_rot_row[1] = v_rot[0][1]; v_rot_row[2] = v_rot[0][2];

  pthread_mutex_lock(&mymutex);
  Message &this_msg = current_ms[him][me];
  pthread_mutex_unlock(&mymutex);  

  if(this_msg.size() == 0)
    this_msg = Message(LABEL_COUNT);
  
  for(int twist=0; twist < TWIST_BIN_COUNT; twist++)
    {
      Msgtmp_type &this_msg_tmp = msg_tmp[twist];
      
      // create msg_tmp, a 3-d matrix version of H (for a single twist)
      copy_message(this_msg_tmp, H, twist * GRID_COUNT * GRID_COUNT * GRID_COUNT);
      
      // do distance transform over each of 3 dimensions
      LOOP_CORRECT(x, y, z, X, Y, Z, 0, 1, 2, v_rot_row, this_msg_tmp, rotation_matrices);
      LOOP_CORRECT(x, z, y, X, Z, Y, 0, 2, 1, v_rot_row, this_msg_tmp, rotation_matrices);
      LOOP_CORRECT(y, z, x, Y, Z, X, 1, 2, 0, v_rot_row, this_msg_tmp, rotation_matrices);
    }
  
  if(TWIST_BIN_COUNT > 1)
    {
      double scratch[TWIST_BIN_COUNT];
      for(int x=0; x<GRID_COUNT; x++)
        for(int y=0; y<GRID_COUNT; y++)
          for(int z=0; z<GRID_COUNT; z++)
            {
              for(int t=0; t<TWIST_BIN_COUNT; t++)
                scratch[t] = msg_tmp[t].get(x)[y][z];
	      
              for(int t=0; t<TWIST_BIN_COUNT; t++)
                {
                  double expected_twist = twist_of_rotation_matrix(rotation_matrix<double>(to_dist(x), to_dist(y), to_dist(z), to_twist(t)) * mehim_pd.rot);
		  
                  double _min = INFINITY;
                  for(int t2=0; t2<TWIST_BIN_COUNT; t2++)
                    {
                      double val = scratch[t2] + twist_dist(t, t2);
                      if(val < _min)
                        _min = val;
                    }
		  
                  msg_tmp[t].get(x)[y][z] = _min;
                }
            }
    }
  
  for(int twist=0; twist < TWIST_BIN_COUNT; twist++)
    copy_message(this_msg, twist * GRID_COUNT * GRID_COUNT * GRID_COUNT, msg_tmp[twist]);
  
  subtract_dc(this_msg);
}


// Andrew's pre-propagation trick: since the only unary constraint in the graph is at node cc_ref,
//  propagate messages out from that node (in depth-first order). 
//
void prepropagate_messages(MessageSet &current_ms, MessageSet &last_ms, const vector<Edge> &edge_list, Pairs &deltas, int cc_ref, int ref_label, const vector< _DMatrix<TYPE> > &rotation_matrices)
{
  stack<int> S;
  map<int, int> parents;
  S.push(ref_label);
  parents.insert(make_pair(ref_label, -1));

  while(!S.empty())
    {
      int node = S.top();
      S.pop();

      // push neighbors of S on stack
      for(NeighborPairs::const_iterator iter = deltas[node].begin(); iter != deltas[node].end(); ++iter)
        if(parents.insert(make_pair(iter->first, node)).second)
          S.push(iter->first);

      // send message from nodes's parent to node
      if(parents[node] != -1)
        compute_message(parents[node], node, cc_ref, deltas, current_ms, current_ms, rotation_matrices);
	
    }

}


void do_iter_edge_subset(int beg_edge, int end_edge, const vector<Edge> &edge_list, int cc_ref, const Pairs &deltas,
                         const MessageSet &last_ms, MessageSet &current_ms, const vector< _DMatrix<TYPE> > &rotation_matrices)
{
  int my_edge_count = end_edge - beg_edge + 1;
  for(int i=beg_edge; i<=end_edge; i++)
    {
      if(my_edge_count < 10000 && !(i%100))
        cerr << ".";
      else if(!(i%1000))
        cerr << ".";
		    
      Edge edge = edge_list[i];
      int me = edge.first, him = edge.second;
		
      compute_message(me, him, cc_ref, deltas, last_ms, current_ms, rotation_matrices);
    }
}

void *bp_thread(void *_p)
{
  BP_Params *p = (BP_Params *)_p;

  do_iter_edge_subset(p->beg_edge, p->end_edge, *p->edge_list, p->cc_ref, *p->deltas, *p->last_ms, *p->current_ms, *p->rotation_matrices);
  return 0;
}

void do_iter_edge_subset(int beg_edge, int end_edge, const vector<Edge> &edge_list, int cc_ref,  Pairs &deltas,
                         MessageSet &last_ms, MessageSet &current_ms, vector< _DMatrix<TYPE> > &rotation_matrices, int thread_count)
{ 
  if(thread_count == 0)
    return;

  pthread_t p;
  BP_Params params(beg_edge, (end_edge-beg_edge)/thread_count + beg_edge, &edge_list, cc_ref, &deltas, &last_ms, &current_ms, &rotation_matrices);
  pthread_create(&p, NULL, &bp_thread, (void *) &params);

  do_iter_edge_subset(params.end_edge+1, end_edge, edge_list, cc_ref,  deltas, last_ms, current_ms, rotation_matrices, thread_count-1);

  pthread_join(p, 0); 
}


void read_pairs_file(const char *pairwise_file, int &photo_count, int &edge_count, Pairs &deltas, vector<Edge> &edge_list)
{
  int field_count=0;
  
  // try to figure out file format (# of fields per pair)
  {
    ifstream ifs(pairwise_file, ios::in);
    ifs >> photo_count >> edge_count;
    if(!ifs.good())
      throw string("cannot open ") + pairwise_file;

    cerr  << "fname = " << pairwise_file << endl;
    double t;
    for(field_count=0; ifs.good(); field_count++) {
      ifs >> t;
      //cerr << t << " ";
      // if (field_count > 0 && field_count % 18 == 0) {
      // 	cerr << "\n";
      // }
    }
    field_count--;

    double i_part;
    if(modf(field_count / double(edge_count), &i_part) != 0 || i_part < 14 )
      {
        cerr << "field_count = " << field_count << endl;
        cerr << "edge_count = " << edge_count << endl;
        cerr << "i_part = " << i_part << endl;
        cerr << "frac part = " << modf(field_count / double(edge_count), &i_part) << endl;
        throw string("can't parse file format of ") + pairwise_file;
      }

    field_count /= edge_count;
  }

  ifstream ifs(pairwise_file_g, ios::in);
  ifs >> photo_count >> edge_count;

  cerr << "pairs file seems to have ";
  if(field_count <= 15) ifs >> PairwiseDiff::IOOptions::NO_CRATIOS; else cerr << "c_ratios, ";
  if(field_count <= 16) ifs >> PairwiseDiff::IOOptions::NO_CONFIDENCE; else cerr << "confidences, ";
  int extra_fields = max(field_count - 17, 0);
  cerr << " and " << extra_fields << " fields ignored " << endl;

  // undirected
  edge_count *=2;
  edge_list = vector<Edge>(edge_count);

  for(int ii=0; ii<edge_count; ii+=2)
    {
      int i, j;
      ifs >> i >> j; // can't combine this line with the next (i and j aren't assigned until after statement returns)
      ifs >> deltas[i][j];

      deltas[j][i].rot = deltas[i][j].rot.transpose();
      deltas[j][i].trans[0] = -deltas[i][j].trans[0];
      deltas[j][i].trans[1] = -deltas[i][j].trans[1];
      deltas[j][i].trans[2] = -deltas[i][j].trans[2];
      deltas[j][i].c1_ratio = deltas[i][j].c2_ratio; 
      deltas[j][i].c2_ratio = deltas[i][j].c1_ratio;

      // burn off extra parameters
      double t;
      for(int k=0; k<extra_fields; k++)
        ifs >> t;
	  
      edge_list[ii] = Edge(i,j);
      edge_list[ii+1] = Edge(j,i);
    }
}


// computes norm(R*t_ij - (t_j - t_i))
// row major matrix, column vectors
// R, t_ij, tj_minus_ti column major (from bundler)
// requires tij and tj_minus_ti be normalized
inline double estimate_geotag_error(double *R, double *t_ij, double *tj_minus_ti) {
  return sqrt(sqr(-tj_minus_ti[0] + R[0] * t_ij[0] + R[1] * t_ij[1] + R[2] * t_ij[2])
              + sqr(-tj_minus_ti[1] + R[3] * t_ij[0] + R[4] * t_ij[1] + R[5] * t_ij[2])
              + sqr(-tj_minus_ti[2] + R[6] * t_ij[0] + R[7] * t_ij[1] + R[8] * t_ij[2]));
}

inline double pan_geotag_error(double *R, double *t_ij, double *tj_minus_ti) {
  double r[3];
  r[0] = R[0] * t_ij[0] + R[1] * t_ij[1] + R[2] * t_ij[2];
  r[1] = R[3] * t_ij[0] + R[4] * t_ij[1] + R[5] * t_ij[2];
  r[2] = R[6] * t_ij[0] + R[7] * t_ij[1] + R[8] * t_ij[2];
  // drop the y coordinate and re-normalize
  double new_r[3];
  new_r[0] = r[0];
  new_r[1] = 0;
  new_r[2] = r[2];
  double norm = sqrt(sqr(new_r[0]) + sqr(new_r[1]) + sqr(new_r[2]));
  new_r[0] /= norm;
  new_r[1] /= norm;
  new_r[2] /= norm;
  // cerr << "old = " << t_ij[0] << " " << t_ij[1] << " "  << t_ij[2] << endl;
  // cerr << "new = " << new_tij[0] << " " << new_tij[1] << " "  << new_tij[2] << endl;
  // ignore y coordinate (should be 0 in both vectors)
  return sqrt(sqr(new_r[0] - tj_minus_ti[0]) + sqr(new_r[2] - tj_minus_ti[2]));
}

// sample t_j - t_i + N where N is a 2-D gaussian with variance geotag_variance_g; stores the result in sampled
// Note: no sampling is done by default, since geotag_variance_g = 0 (i.e. it always returns just t_j - t_i).  
void sample_translation_differences(DMatrix &t_i, DMatrix &t_j, int num_trials, vector<DMatrix> &sampled) {
  assert((uint)num_trials == sampled.size());
  assert(num_trials == 0 || (sampled[0].rows() == 3 && sampled[0].cols() == 1));
  DMatrix mean = t_j - t_i;
  //double difference_stddev = sqrt(2*geotag_variance_g);
  double difference_stddev = sqrt(2)*geotag_stdev_g;
  for (int t = 0; t < num_trials; t++) {
    double mag = rand_normal(0, difference_stddev);
    double theta = drand48()*2*M_PI;
    double x = mean[0][0] + mag*cos(theta);
    double z = mean[2][0] + mag*sin(theta);
    DMatrix &v = sampled[t];
    double norm = sqrt(square(x) + square(z));
    v[0][0] = x/norm;
    v[1][0] = 0;
    v[2][0] = z/norm;
  }
}

// sample geotag directions (t_j - t_i)/norm(t_j - t_i); for each
// sample, compute its distance from tij; add this distance to H
void add_to_unary_with_sampling(vector<_DMatrix<TYPE> > &rotation_matrices, DMatrix &ti,
                                DMatrix &tj, _DMatrix<PD_TYPE> &tij, Distribution &H, double weight, double max_error) {
  static const int NUM_TRIALS = 1;
  static vector<DMatrix> trans_diffs(NUM_TRIALS);
  static vector<double*> rots(rotation_matrices.size());
  static bool first = true;
  assert((int)rotation_matrices.size() == LABEL_COUNT);
  if (first) {
    first = false;
    for (uint i = 0; i < trans_diffs.size(); i++) {
      trans_diffs[i] = DMatrix(3,1);
    }
    for (uint i = 0; i < rotation_matrices.size(); i++) {
      if (is_state_on_sphere((int)i)) {
        double *R = new double[9];
        for (int ii = 0; ii < 3; ii++) {
          for (int j = 0; j < 3; j++) {
            R[ii*3+j] = rotation_matrices[i][ii][j];
          }
        }
        rots[i] = R;
      }
      else {
        rots[i] = 0;
      }
    }
  }
  // sample t_j - t_i assuming t_i and t_j come from gaussians
  sample_translation_differences(ti, tj, num_unary_samples_g, trans_diffs);
  double tij_a[3] = {tij[0][0], tij[1][0], tij[2][0]};
  for (uint i = 0; i < trans_diffs.size(); i++) {
    // take a sample
    double td[3] = {trans_diffs[i][0][0], trans_diffs[i][1][0], trans_diffs[i][2][0]};
    for (int label = 0; label < LABEL_COUNT; label++) {
      if (is_state_on_sphere(label)) {
        assert((uint)label < H.size());
        if (geotag_pan_g) {
          H[label] += weight/trans_diffs.size()*min(max_error, pan_geotag_error(rots[label], tij_a, td));
        }
        else {
          H[label] += weight/trans_diffs.size()*min(max_error, estimate_geotag_error(rots[label], tij_a, td));
        }
      }
    }
  }
}

// struct geotag_thread_params {
//   vector<Edge> *edge_list;
//   int first_edge, num_edges;
// };

// void geotag_unary_thread(void *void_params) {
//   geotag_thread_params *params = (geotag_thread_params*)void_params;
//   for (int i = params->first_edge; i < params->num_edges; i++) {
    
//   }
// }

// void make_geotag_unary_potentials(vector<Edge> &edge_list) {
//   vector<pthread_t> threads(thread_count_g);
//   int edges_per_thread = ((int)edge_list.size())/thread_count_g;
//   for (int i = 0; i < thread_count_g; i++) {
//     geotag_thread_params params;
//     params.first_edge = i*edges_per_thread;
//     params.num_edges = ((i == thread_count_g - 1) ? edge_list.size() - params.first_edge : edges_per_thread);
//     params.edge_list = &edge_list;
//     pthread_create(&threads[i], NULL, &geotag_unary_thread, (void *) &params);
//   }
//   for (int i = 0; i < thread_count_g; i++) {
//     pthread_join(&threads[i], 0);
//   }
// }

int main(int argc, char *argv[])
{
  try {
    parse_opts(argc, argv);

    // initialize look-up table
    is_state_on_sphere(0);

    int photo_count, edge_count;
    Pairs deltas;
    vector<Edge> edge_list;
    read_pairs_file(pairwise_file_g, photo_count, edge_count, deltas, edge_list);

    if(!no_remove_cut_g)
      {
        edge_list = remove_cut_edges(edge_list, deltas);
        edge_count = edge_list.size();
      }

    hash_map<int, DMatrix> known_locations;
    if (geoplanar_g) {
      read_2d_locations_file(geoplanar_g, known_locations);
    }
    int cc_ref=0;
    set<int> reachable = check_reachability(photo_count, edge_count, edge_list, cc_ref, random_ref_g, ccref_highest_degree_g);
    if(ref_node_g != -1)
      cc_ref = ref_node_g;
    cerr << reachable.size() << " reachable, reference node = " << cc_ref << endl;
    cout << "ref " << cc_ref << endl;
    const int ref_label = pack(from_dist(0), from_dist(0), from_dist(-1), from_twist(0));

    // now remove edges that aren't involved in the connected component of interest
    edge_list = remove_unconnected(edge_list, reachable, deltas);
    cerr << "removed " << edge_count - edge_list.size() << " irrelevant edges (" << edge_list.size() << " remaining)" << endl;
    edge_count = edge_list.size();

    if(write_edges_g)
      {
        ofstream ofs((string(pairwise_file_g) + ".cut").c_str());

        for(int i=0; i<edge_count; i++)
          {
            int me = edge_list[i].first, him = edge_list[i].second;
            if(him < me)
              continue;
            ofs << me << " " << him << " ";;
            for(int j=0; j<9; j++)
              ofs << deltas[me][him].rot[j] << " ";
            for(int j=0; j<3; j++)
              ofs << deltas[me][him].trans[j] << " ";
            ofs << deltas[me][him].c1_ratio << " " << deltas[me][him].c2_ratio << endl;
          }
      }

    hash_map<int, PairwiseDiff> gt;
    if(gt_g)
      {
        ifstream ifs_gt(gt_g, ios::in);
        ifs_gt >> PairwiseDiff::IOOptions(PairwiseDiff::IOOptions::NO_CRATIOS);
        while(ifs_gt.good())
          {
            int c;
            ifs_gt >> c;
            if(!ifs_gt.good())
              break;
            ifs_gt >> gt[c];
          }
	
        cerr << "read " << gt.size() << " gt entries " << endl;
      }

    map<int, double> vanish_tilts;
    if(vanish_g)
      {
        cerr << "loading tilts and twists... ";
        ifstream ifs(vanish_g, ios::in);
        cerr << "vanish file = " << vanish_g << endl;
        assert(ifs.good());

        int cid;
        double tilt, twist, conf;

        while(ifs.good())
          {
            ifs >> cid >> tilt >> twist >> conf;
            vanish_tilts[cid] = tilt / 180.0 * M_PI;
          }
        cerr << "done. " << endl;
      }

    if(param_est_g)
      {
        int ii=0;
        DMatrix tilt_errors(1, photo_count);

        for(int cid=0; cid<photo_count; cid++)
          {
            if(gt.find(cid) == gt.end())
              continue;

            if(vanish_tilts.find(cid) != vanish_tilts.end())
              tilt_errors[0][ii] = sin(tilt_of_rotation_matrix(gt[cid].rot)) - sin(vanish_tilts[cid]);

            ii++;
          }

        tilt_errors = tilt_errors.extract(DRect(0, 0, 0, ii-1));
        cout << "tilt stddev = " << (tilt_errors.transpose().covariance()) << endl;

        ii=0;
        DMatrix edge_errors(edge_count, 3);
        for(int i=0; i<edge_count; i++)
          {
            int me = edge_list[i].first, him = edge_list[i].second;
            if(gt.find(me) == gt.end() || gt.find(him) == gt.end())
              continue;

            DMatrix predicted = compute_zerotwist_rot<double>(-gt[me].rot[0][2], -gt[me].rot[1][2], -gt[me].rot[2][2]) * deltas[me][him].rot;
            edge_errors.set_row(i, -predicted.extract_col(2).transpose() - (-gt[him].rot.extract_col(2).transpose()));

            ii++;
          }

        edge_errors = edge_errors.extract(DRect(0, 0, ii-1, 2));
        DMatrix edge_cov = edge_errors.covariance();

        cout << "edge cov = " << edge_cov << endl;

        return 0;
      }


    if(perfect_pairs_g)
      {
        ofstream ofs(perfect_pairs_g);

        for(int i=0; i<edge_count; i++)
          {
            int me = edge_list[i].first, him = edge_list[i].second;
            if(gt.find(me) == gt.end() || gt.find(him) == gt.end())
              continue;

            DMatrix rot = gt[me].rot.inverse() * gt[him].rot;
            ofs << me << " " << him << " ";
            for(int j=0; j<9; j++)
              ofs << rot[0][j] << " ";
            for(int j=0; j<3; j++)
              ofs << 0 <<  " ";
            ofs << deltas[me][him].c1_ratio << " " << deltas[me][him].c2_ratio << endl;
          }

        return 0;
      }


    // set up unary costs
    cerr << "setting up unary potentials..." << endl;
    Distribution generic_unary(LABEL_COUNT), ref_unary(LABEL_COUNT);
    const int ref_x = from_dist(0), ref_y = from_dist(0), ref_z = from_dist(-1);

    for(int i=0; i< LABEL_COUNT; i++)
      {
        if(!is_state_on_sphere(i))
          generic_unary[i] = ref_unary[i] = INFINITY;
        else
          {
            int x = unpack_t0(i), y = unpack_t1(i), z = unpack_t2(i);
            ref_unary[i] = sqr(x-ref_x) + sqr(y-ref_y) + sqr(z-ref_z);
          }
      }

    if(!vanish_g)
      {
        unary_dists = DistributionMap(photo_count);
        if (geoplanar_g) {
          for(int cid=0; cid<photo_count; cid++) {
            unary_dists[cid] = ref_unary;
          }
        }
        else {
          for(int cid=0; cid<photo_count; cid++) {
            if(cid == cc_ref)
              unary_dists[cid] = ref_unary;
            else
              unary_dists[cid] = generic_unary;
          }
        }
      }
    else
      {
        unary_dists = DistributionMap(photo_count);
        for(int cid=0; cid<photo_count; cid++)
          {
            if(vanish_tilts.find(cid) == vanish_tilts.end())
              unary_dists[cid] = generic_unary;
            else
              {
                double vanish_y = sin(vanish_tilts[cid]);
                Distribution &my_unary = unary_dists[cid];
                my_unary = Distribution(LABEL_COUNT);
                for(int i=0; i< LABEL_COUNT; i++)
                  {
                    if(!is_state_on_sphere(i))
                      my_unary[i] = INFINITY;
                    else
                      my_unary[i] = min(sqr(unary_scale_g * (to_dist(unpack_t1(i)) - vanish_y)), unary_trunc_g);
                  }

                if(!geoplanar_g && cid == cc_ref)
                  my_unary = ref_unary;
              }
          }
      }

    // Precompute the map from viewing directions to rotation matricies
    cerr << "precomputing rotation matrices... " << TWIST_BIN_COUNT << "," << GRID_COUNT << endl;
    vector< _DMatrix<TYPE> > rotation_matrices(LABEL_COUNT);
    for(int tw=0; tw < TWIST_BIN_COUNT; tw++)
      for(int x=0; x < GRID_COUNT; x++)
        for(int y=0; y < GRID_COUNT; y++)
          for(int z=0; z < GRID_COUNT; z++)
            rotation_matrices[pack(x, y, z, tw)] = rotation_matrix<TYPE>(to_dist(x), to_dist(y), to_dist(z), to_twist(tw));
    cerr << "done" << endl;
    if (geoplanar_g) {
      cerr << "Building geoplanar..." << endl;
      int edges_with_geotag = 0;
      for (uint ei = 0; ei < edge_list.size(); ei++) {
        int me = edge_list[ei].first;
        int him = edge_list[ei].second;
        if (known_locations.find(me) != known_locations.end()
            && known_locations.find(him) != known_locations.end()) {
          edges_with_geotag++;
        }
      }

      cerr << "fraction of edges with 2 geotags = " << ((double)edges_with_geotag) / edge_list.size() << endl;
      double weight = scale_geotag_unary_g ? (scaled_geotag_weight_g * ((double)edges_with_geotag) / edge_list.size()) : unary_geotag_weight_g;
      cerr << "geotag unary edge weight = " << weight << endl;
      ofstream geotag_priors("rotbp_geotag_priors");
      cerr << "known locations: " << known_locations.size() << endl;
      for (uint ei = 0; ei < edge_list.size(); ei++) {
        int me = edge_list[ei].first;
        int him = edge_list[ei].second;
        if (known_locations.find(me) != known_locations.end()
            && known_locations.find(him) != known_locations.end()) {
          DMatrix t_i = known_locations[me];
          DMatrix t_j = known_locations[him];
          //	  DMatrix &t_ij = deltas[me][him].trans;
          _DMatrix<PD_TYPE>  t_ij(3,1,deltas[me][him].trans);
          //assert(abs(vector_norm_sq(t_ij) - 1) <= 0.0001);
          if (ei % ((int)ceil(edge_list.size()/20.0)) == 0) {
            cerr << ".";
          }
          if (write_geo_priors_g) {
            DMatrix diff = t_j - t_i;
            diff /= vector_norm(diff);
            geotag_priors << me << " " << him << " " << t_ij[0][0] << " " << t_ij[1][0] << " " << t_ij[2][0]
                          << " " << diff[0][0] << " " << diff[1][0] << " " << diff[2][0] << endl;
          }
          add_to_unary_with_sampling(rotation_matrices, t_i, t_j, t_ij, unary_dists[me], weight, unary_geotag_trunc_g);
        }
      }
      geotag_priors.close();
    }

    cerr << "out of geotag priors" << endl;

    if(rescore_g)
      {
        ifstream ifs(rescore_g);
        map<int, DMatrix> estimated_dirs;
        cerr << "reading " << rescore_g << endl;
        while(ifs.good())
          {
            int cid;
	    
            // handle two formats: cid followed by 3 viewing directions, or cid followed by 9-entry rotation matrix
            string str;
            getline(ifs, str);
            istringstream iss(str);
	    
            iss >> cid;
            if(!ifs.good())
              break;
	    
            estimated_dirs[cid] = DMatrix(4,1);
            estimated_dirs[cid] = 0;
            if(count(str.begin(), str.end(), ' ') > 4)
              {
                DMatrix tmp(3,3);
                iss >> tmp[0][0] >> tmp[0][1] >> tmp[0][2] >> tmp[1][0] >> tmp[1][1] >> tmp[1][2] >> tmp[2][0] >> tmp[2][1] >> tmp[2][2];

                estimated_dirs[cid][0][0] = -tmp[0][2];
                estimated_dirs[cid][1][0] = -tmp[1][2];
                estimated_dirs[cid][2][0] = -tmp[2][2];
                estimated_dirs[cid][3][0] = twist_of_rotation_matrix(tmp);
              }
            else
              iss >> estimated_dirs[cid][0][0] >> estimated_dirs[cid][0][1] >> estimated_dirs[cid][0][2];
          }

        double cost = compute_labeling_cost(estimated_dirs, deltas, edge_list, 0);
        double err = 0, err2=0, twist_error=0;
        if(gt_g)
          err2 = compute_gt_error(err, gt, reachable, estimated_dirs, photo_count, 0, twist_error);

        cout << "iter " << 0 << " of " << 0 << " " << cost << " " << err << " " << err2 << " " << twist_error << endl;

        return 0;
      }

    if(discrete_g)
      {
        //      initialize message buffers
        MessageSet last_ms; 
        for(int i=0; i<edge_count; i++)
          {
            if(last_ms[edge_list[i].first][edge_list[i].second].size() == 0)
              last_ms[edge_list[i].first][edge_list[i].second] = Message(LABEL_COUNT);
          }
        MessageSet current_ms = last_ms;
        cerr << "done" << endl;

        // Prepropagate messages using BFS from the reference node
        if(prepropagate_g)
          {
            prepropagate_messages(current_ms, last_ms, edge_list, deltas, cc_ref, ref_label, rotation_matrices);
            last_ms = current_ms;
          }

        for(int ii=0; ii<max_ii_g; ii++)
          {
            cerr << ii << endl;
            int beg_edge = 0, end_edge = edge_count-1;
            if(thread_count_g > 1)
              do_iter_edge_subset(beg_edge, end_edge, edge_list, cc_ref, deltas, last_ms, current_ms, rotation_matrices, thread_count_g);
            else
              do_iter_edge_subset(beg_edge, end_edge, edge_list, cc_ref, deltas, last_ms, current_ms, rotation_matrices);

            last_ms = current_ms;
	    
            map<int, DMatrix> estimated_dirs;
            cerr << endl << "computing MAP estimates and energy..." << endl;
            hash_map<int, vector<TYPE> > final_distribution = compute_final_distribution(current_ms, edge_list, cc_ref); 

            vector<int> labels = find_map_estimates(final_distribution, photo_count, estimated_dirs, ii);
            double cost = compute_labeling_cost(estimated_dirs, deltas, edge_list, ii);
            double err = 0, err2=0, twist_error=0;
            if(gt_g)
              err2 = compute_gt_error(err, gt, reachable, estimated_dirs, photo_count, ii, twist_error);

            cout << "iter " << ii << " of " << max_ii_g << " " << cost << " " << err << " " << err2 << " " << twist_error << endl;
          }
      }
    else
      {
        // continuous code snipped for now
      }
    
  } catch(const string &str)
    {
      cerr << str << endl;
    }
}
