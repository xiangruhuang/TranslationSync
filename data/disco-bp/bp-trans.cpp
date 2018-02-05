#include "DImage.h"
#include <vector>
#include <istream>
#include <fstream>
//#include <io.h>
#include <DistanceTransform.h>
#include <algorithm>
#include <numeric>
#include <map>
#include <DImageIO.h>
#include <DrawText.h>
#include <ext/hash_set>
#include "bp-trans-common.h"
#include <iomanip>
#include <DProfile.h>

#ifdef HADOOP_VERSION
#include "hadoop/Pipes.hh"
#include "hadoop/TemplateFactory.hh"
#include "hadoop/StringUtils.hh"
#include <hadoop_utils.h>
#endif

using namespace std;

// input files
char *pairwise_file_g = 0, *unary_file_g = 0;

// # of BP iterations to run
int max_ii_g = 100;

// relative weighting of unary potentials
double unary_weighting_g = 100;

// truncation points for unary, pairwise costs
double pair_trunc_g = 10, prior_trunc_g = 10;

// ground truth file
char *gt_g = 0;

// if nonzero, initialize unary potentials based on ground truth. value indicates
//  fraction of nodes that should have unary potentials.
double known_from_gt_g = 0;

// hacks to prevent collapsing to trivial solution (all cameras at same point)
double anticollapse_g = 0, anticollapse_penalty_g = 0;

// only use high-confidence geotags for initializing unary potentials
bool geo_highconf_only_g = false;

// specifies whether to view pairwise constraints as rays (false) or lines (true)
bool bidir_g = true;

// set >1 to enable "banding", which improves the accuracy of the objective function
//  estimates.
int bands_g = 1;

// save memory (but increase compute time) by not caching priors
bool precompute_priors_off_g = false;

// set to true to disable removing cut edges
bool no_remove_cut_edges_g = false;

// For debugging...
double stretch_g = 4;
int bufsz_g = 4000;
int node_g = -1;
char *rescore_g = 0;
int thread_count_g = 1;
bool dump_messages_g = false;

// translate entire coordinate frame 
DPoint translate_g(0,0);

// set to true to estimate size of label space from geotags
bool est_label_space_g = false;

// size of label space per dimension
int bins_g = 201;

// selects between different gt file formats
bool gt_trans_only_g = false;

int GRID_COUNT2_g = 71;
double DIST_GRAN2_g = 0.08;

using namespace __gnu_cxx;
using namespace std;


void parse_opts(int argc, char *argv[])
{
  int ii=1;
  for(ii=1; ii<argc; ii++)
    {
      if(!strcmp(argv[ii], "--iters"))
	max_ii_g = atoi(argv[++ii]);
      else if(!strcmp(argv[ii], "--rescore"))
	rescore_g = argv[++ii];
      else if(!strcmp(argv[ii], "--unweight"))
	unary_weighting_g = atof(argv[++ii]);
      else if(!strcmp(argv[ii], "--bands"))
	bands_g = atoi(argv[++ii]);
      else if(!strcmp(argv[ii], "--gt"))
	gt_g = argv[++ii];
      else if (!strcmp(argv[ii], "--gtformat")) {
	ii++;
	if (!strcmp(argv[ii], "trans"))
	  gt_trans_only_g = true;
	else if (!strcmp(argv[ii], "rottrans"))
	  gt_trans_only_g = false;
	else {
	  cerr << "unknown ground truth file format " << argv[ii] << endl;
	  exit(1);
	}
      }
      else if(!strcmp(argv[ii], "--nocutedgefilter"))
	no_remove_cut_edges_g = true;
      else if(!strcmp(argv[ii], "--threads"))
	thread_count_g = atoi(argv[++ii]);
      else if(!strcmp(argv[ii], "--knownfromgt"))	
	known_from_gt_g = atof(argv[++ii]);
      else if(!strcmp(argv[ii], "--node"))
	node_g = atoi(argv[++ii]);
      else if(!strcmp(argv[ii], "--pairtrunc"))
	pair_trunc_g = atof(argv[++ii]);
      else if(!strcmp(argv[ii], "--priortrunc"))
	prior_trunc_g = atof(argv[++ii]);
      else if(!strcmp(argv[ii], "--bidir"))
	bidir_g = false;
      else if(!strcmp(argv[ii], "--highconf"))
	geo_highconf_only_g = true;
      else if(!strcmp(argv[ii], "--png"))
	{
	  stretch_g = atof(argv[++ii]);
	  bufsz_g = atoi(argv[++ii]);
	}
      else if(!strcmp(argv[ii], "--anticollapse"))
	{
	  anticollapse_g = atof(argv[++ii]);
	  anticollapse_penalty_g = atof(argv[++ii]);
	}
      else if(!strcmp(argv[ii], "--translate"))
	{
	  int row = atoi(argv[++ii]);
	  int col = atoi(argv[++ii]);
	  translate_g = DPoint(row, col);
	}
      else if(!strcmp(argv[ii], "--estimatelabelspace"))
	{
	  est_label_space_g = true;
	  bins_g = atoi(argv[++ii]);
	  assert(bins_g % 2 == 1);
	}
      else if(!strcmp(argv[ii], "--labelspace"))
	{
	  GRID_COUNT2_g =  atoi(argv[++ii]); //101; // 21; // 101; // 41; // 131; // 131 breaks it -- WHY????
	  DIST_GRAN2_g = atof(argv[++ii]); // 0.15; // 0.025; // 3.0; // 7.0; // 0.1;
	  assert(GRID_COUNT2_g % 2 == 1);
	}
      else if(!strcmp(argv[ii], "--dumpmessages"))
	dump_messages_g = true;
      else
	break;
    }

  if(ii < argc)
    pairwise_file_g = argv[ii++];

  if(ii < argc)
    unary_file_g = argv[ii++];

  cerr << pairwise_file_g << " "<< unary_file_g << endl;

  if(argc != ii)
    {
      cerr << "bad commandline options." << endl;
      exit(1);
    }
}

void *bp_thread(void *_p)
{
  BPtrans_Params *p = (BPtrans_Params *)_p;

  p->inf->do_iter_edge_subset(p->beg_edge, p->end_edge, *p->edge_list, *p->deltas, *p->last_ms, *p->current_ms, *p->known_locations, *p->D_planes, *p->final, p->ii);
  return 0;
}



void do_iter_edge_subset(int beg_edge, int end_edge, const vector<Edge> &edge_list, Pairs &deltas, 
			 PackedMessageSet &last_ms, PackedMessageSet &current_ms, map<int, DMatrix> &known_locations, PriorsMap &D_planes, 
			 map< int, Message > &final, int ii, int thread_count, const BPTrans_Inference *bp)
{ 
  if(thread_count == 0)
    return;

  pthread_t p;
  BPtrans_Params params(beg_edge, (end_edge-beg_edge)/thread_count + beg_edge, &edge_list, &deltas, &last_ms, &current_ms, &known_locations, &D_planes, &final, ii, bp);
  pthread_create(&p, NULL, &bp_thread, (void *) &params);

  do_iter_edge_subset(params.end_edge+1, end_edge, edge_list, deltas, last_ms, current_ms, known_locations, D_planes, final, ii, thread_count-1, bp);

  pthread_join(p, 0); 
}

void estimate_labelspace_size(const hash_map<int, UnaryPosition> &unary, const BPTrans_Inference &bp)
{
  _DRect<double> bound, bound2;
  
  for(hash_map<int, UnaryPosition>::const_iterator iter = unary.begin(); iter != unary.end(); ++iter)
    {
      if(fabs(iter->second.position).sum() && iter->second.confidence > 0.0)
	bound = bounding_rectangle(bound, _DRect<double>(_DPoint<double>(iter->second.position[0][0], iter->second.position[0][2]), 
							 _DPoint<double>(iter->second.position[0][0], iter->second.position[0][2])));
      bound2 = bounding_rectangle(bound2, _DRect<double>(_DPoint<double>(iter->second.position[0][0], iter->second.position[0][2]), 
							 _DPoint<double>(iter->second.position[0][0], iter->second.position[0][2])));
    }
  
  _DPoint<double> b1 = elementwise_max(-bound.top_left(), bound.bottom_right());
  double max_b1 = max(b1.row(), b1.col()) * 1.2;
  cout << bound << " " << bins_g << " " << bp.grid_count()/2/max_b1 << endl;
  
  _DPoint<double> b2 = elementwise_max(-bound2.top_left(), bound2.bottom_right());
  double max_b2 = max(b2.row(), b2.col()) * 1.2;
  cout << bound2 << " " << bins_g << " " << bp.grid_count()/2/max_b2 << endl;
  
}

void rescore(const char *rescore_file, const BPTrans_State &bps)
{
  hash_map<int, DMatrix> estimates;
  
  ifstream ifs(rescore_file, ios::in);
  while(ifs.good())
    {
      int cid; 
      double r1, r2, r3;
      ifs >> cid >> r1 >> r2 >> r3;
      
      if(!ifs.good())
	break;
      
      DMatrix vec(1,3);
      vec[0][0] = r1; vec[0][1] = r2; vec[0][2] = r3;
      estimates[cid] = vec;
    }
  
  double corr_err;
  double err = compare_with_gt(bps.gt, bps.photo_count, estimates, bps.reachable, bps.known_locations, corr_err);
  cerr << "DONE RESCORE " << err << endl;
  double cost = 0;
  cout << "iter " << 0 << " of " << 0 << " " << cost << " "  << err << " " << corr_err << endl;
}

#ifndef HADOOP_VERSION

int main(int argc, char *argv[])
{
  parse_opts(argc, argv);

  BPTrans_Inference bp(prior_trunc_g, unary_weighting_g, pair_trunc_g, bands_g, anticollapse_g, 
		       anticollapse_penalty_g, bidir_g, precompute_priors_off_g, translate_g, GRID_COUNT2_g, DIST_GRAN2_g);

  try {

    BPTrans_State bps;
    bps.initialize(pairwise_file_g, unary_file_g, gt_g, no_remove_cut_edges_g, translate_g, gt_trans_only_g,
		   known_from_gt_g, geo_highconf_only_g, bp);

    if(est_label_space_g)
      {
	estimate_labelspace_size(bps.unary, bp);
	return 0;
      }
    else if(rescore_g)
      {
	rescore(rescore_g, bps);
	return 0;
      }

    BPTrans_Debug debug(stretch_g, bp, bufsz_g, bps.gt);
    debug.make_gt_image(bps.known_locations, bps.photo_count, bps.reachable);

    // FIXME - this is a hack. if there are a lot of known camera positions, turn off pre-computing the priors.
    //  (saves memory, at the expense of speed)
    if(bps.known_locations.size() > 5000)
      {
	cerr << "WARNING: Turning off prior pre-computation to save memory." << endl;
	precompute_priors_off_g = true;
	bp.disable_precompute_priors();
      }

    // calculate priors
    PriorsMap D_planes;
    if(!precompute_priors_off_g)
      for(map<int, DMatrix>::const_iterator iter = bps.known_locations.begin(); iter != bps.known_locations.end(); ++iter)
	D_planes[iter->first] = bp.compute_priors(iter->second);
      
    // initialize message buffers
    PackedMessageSet last_ms(bps.photo_count);
    PackedMessageSet current_ms = last_ms;
    map< int, Message > final; //(photo_count, vector<double>(LABEL_COUNT2));

    for(int ii=0; ii<max_ii_g; ii++)
      {
	cerr << ii << endl;

	int beg_edge = 0, end_edge = bps.edge_count - 1;
	if(thread_count_g > 1)
	  do_iter_edge_subset(beg_edge, end_edge, bps.edge_list, bps.deltas, last_ms, current_ms, bps.known_locations, D_planes, final, ii, thread_count_g, &bp);
	else
	  bp.do_iter_edge_subset(beg_edge, end_edge, bps.edge_list, bps.deltas, last_ms, current_ms, bps.known_locations, D_planes, final, ii);

	last_ms = current_ms;
	cerr << endl;

	// compute MAP estimates from last iteration
	//   0th iteration labels and costs aren't valid, because state likelihood distributions for iteration ii are
	//    computed on the ii+1-th iteration.
	if(ii > 0)
	  {
	    vector<int> labels = bp.map_estimate(final, bps.photo_count);
	    for(int i=0; i<bps.photo_count; i++)
	      cout << "labels " << ii << " " << i << " " << labels[i] << " " << bp.to_dist2(bp.unpack_v0(labels[i])) << " " << bp.to_dist2(bp.unpack_v1(labels[i])) << " " << bp.to_dist2(bp.unpack_v2(labels[i])) << endl;

	    double corr_err;
	    double cost = bp.compute_iter_energy(labels, bps.deltas, bps.edge_list);
	    double err = bp.compute_gt_error(labels, bps.deltas, bps.gt, corr_err, bps.reachable, bps.known_locations);

	    cout << "iter " << ii << " of " << max_ii_g << " " << cost << " "  << err << " " << corr_err << endl;

	    debug.make_iter_image(bps.known_locations, labels, bps.photo_count, bps.reachable, ii);
	  } 

	if(dump_messages_g)
	  {
	    char temp[1024];
	    sprintf(temp, "iter_%03d_messages.txt", ii);
	    ofstream ofs(temp);

	      for(PackedMessageSet::const_iterator dest_iter = current_ms.begin(); dest_iter != current_ms.end(); ++dest_iter)
		{
		  int new_dest = dest_iter->first;
		  
		  for(hash_map<int, PackedMessage>::const_iterator src_iter = dest_iter->second.begin(); src_iter != dest_iter->second.end(); ++src_iter)
		    {
		      int new_src = src_iter->first;
		      const PackedMessage &pm = src_iter->second;

		      ofs << new_dest << " " << new_src << " " << setprecision(10) << pm << endl;
		    }
		}
	  }
      }
      
  } catch(const string &str)
    {
      cerr << str << endl;
    }
}

#else

void parse_hadoop_cmdline(const string &cmdline)
{
  char **argv = new char *[cmdline.size()];
  
  istringstream iss(cmdline);
  string word;
  int i=1;
  for(i=1; iss.good(); i++)
    {
      iss >> word;
      argv[i] = new char[word.size() + 1];
      strcpy(argv[i], word.c_str());
    }
  int argc = i;
  
  parse_opts(argc, argv);
}

// The reducer takes a set of incoming messages for a node, and outputs a set of outgoing messages.
//
// input: dest node_count src1 msg_len1 msg... src2 msg_len2 msg... ...
// output: (dest, src msg_len msg) ...
//
class BPTrans_Mapper : public HadoopPipes::Mapper  
{
public:
  BPTrans_Mapper(HadoopPipes::TaskContext& context) 
  { 
    cerr << "init mapper" << endl;
    string cmdline = get_init_string("bptrans.cmdline", context);

    parse_hadoop_cmdline(cmdline);

    precompute_priors_off_g = true;
    bp = BPTrans_Inference(prior_trunc_g, unary_weighting_g, pair_trunc_g, bands_g, anticollapse_g, 
			 anticollapse_penalty_g, bidir_g, precompute_priors_off_g, translate_g, GRID_COUNT2_g, DIST_GRAN2_g);

    bps.initialize(pairwise_file_g, unary_file_g, gt_g, no_remove_cut_edges_g, translate_g, gt_trans_only_g,
		   known_from_gt_g, geo_highconf_only_g, bp);
    cerr << "done init mapper" << endl;
  }
  
  void map(HadoopPipes::MapContext& context) 
  {   
    cerr << "in mapper" << endl;
    PriorsMap D_planes;
      
    // initialize message buffers
    PackedMessageSet last_ms;
    std::map< int, Message > final; 

    //    cerr << context.getInputValue() << endl;

    /*
    // FIXME : HACK to remove problematic "inf"'s 
    string str2 = context.getInputValue();
    int pos1=0;
    while(pos1 != string::npos)
      {
	pos1 = str2.find("inf", pos1);
	if(pos1 != string::npos)
	  str2 = str2.replace(pos1, 3, "10e100");
      }
    cerr << str2 << endl;
    */

    istringstream iss(context.getInputValue());
    
    int node_count;
    int dest;
    iss >> dest;
    iss >> node_count;

    int iter=1;
    vector<Edge> edge_list2; //(node_count);

    for(int i=0; i<node_count; i++)
      {
	PackedMessage pm;
	int src, msg_flag;
	iss >> src >> msg_flag;
	
	if(msg_flag)
	  {
	    iss >> pm;
	    last_ms[dest][src] = pm;
	  }

	// make sure we actually care about this edge
	Edge e(dest, src);
	if( binary_search(bps.edge_list.begin(), bps.edge_list.end(), e, compare_first) )
	  edge_list2.push_back(e);

	if(!msg_flag)
	  iter=0;
      }

    PackedMessageSet current_ms;

    // FIXME: could possibly use multi-threaded version here
    int ii=iter; // fake iteration number

    node_count = edge_list2.size();

    cout << "dest " << dest << " node count " << node_count << endl;

    bp.do_iter_edge_subset(0, node_count-1, edge_list2, bps.deltas, last_ms, current_ms, bps.known_locations, D_planes, final, ii);
    cerr << endl;
    
    for(PackedMessageSet::const_iterator dest_iter = current_ms.begin(); dest_iter != current_ms.end(); ++dest_iter)
      {
	int new_dest = dest_iter->first;
	
	assert(dest_iter->second.size() == 1);
	assert(dest_iter->second.begin()->first == dest);
	const PackedMessage &pm = dest_iter->second.begin()->second;

	ostringstream key, val;
	key << new_dest;
	val << dest << " 1 ";
	val << setprecision(10) << pm;

	context.emit(key.str(), val.str());
      }

    if(ii > 0)
      {
	assert(final.size() == 1);
	ostringstream oss1;
	//	oss1 << dest << " " <<  setprecision(10) <<  final[dest] << endl;
	int label = bp.map_estimate_singlecamera(final[dest]);
	oss1 << dest << " " << label;
	context.emit("-1", oss1.str());
      }
  }
  
  BPTrans_State bps;
  BPTrans_Inference bp;
};

// The reducer just collects all of the messages destined for each node into a single
//  key-value pair
//
// input:  dest src msg_len msg...
// output: dest node_count src1 msg_len1 msg... src2 msg_len2 msg... ...
//
class BPTrans_Reducer: public HadoopPipes::Reducer
{
public:
  BPTrans_Reducer(HadoopPipes::TaskContext& context) 
  { 
    string cmdline = get_init_string("bptrans.cmdline", context);

    parse_hadoop_cmdline(cmdline);

    precompute_priors_off_g = true;
    bp = BPTrans_Inference(prior_trunc_g, unary_weighting_g, pair_trunc_g, bands_g, anticollapse_g, 
			 anticollapse_penalty_g, bidir_g, precompute_priors_off_g, translate_g, GRID_COUNT2_g, DIST_GRAN2_g);

    bps.initialize(pairwise_file_g, unary_file_g, gt_g, no_remove_cut_edges_g, translate_g, gt_trans_only_g,
		   known_from_gt_g, geo_highconf_only_g, bp);
  }

  void reduce(HadoopPipes::ReduceContext& context) 
  {
    int ii=get_init_param("bptrans.iternum", context);
    if(atoi(context.getInputKey().c_str()) == -1)
      {
	DProfile prof(10);
	prof.begin(0);
	ofstream ofs("iter_output");
	vector<int> labels(bps.photo_count);
	
	prof.begin(1);
	while(context.nextValue())
	  {
	    int node;
	    int label;
	    //	    Message msg;

	    istringstream iss(context.getInputValue());
	    //	    iss >> node >> msg;
	    iss >> node >> label;
	    labels[node] = label;
	    //	    final[node] = msg;

	  }
	prof.end(1);
	    
	// compute MAP estimates from last iteration
	prof.begin(2);
	//	vector<int> labels = bp.map_estimate(final, bps.photo_count);
	prof.end(2);

	prof.begin(3);
	for(int i=0; i<bps.photo_count; i++)
	  ofs << "labels " << ii << " " << i << " " << labels[i] << " " << bp.to_dist2(bp.unpack_v0(labels[i])) << " " << bp.to_dist2(bp.unpack_v1(labels[i])) << " " << bp.to_dist2(bp.unpack_v2(labels[i])) << endl;
	prof.end(3);
	
	double corr_err;
	prof.begin(4);
	double cost = bp.compute_iter_energy(labels, bps.deltas, bps.edge_list);
	prof.end(4);
	prof.begin(5);
	double err = bp.compute_gt_error(labels, bps.deltas, bps.gt, corr_err, bps.reachable, bps.known_locations, ofs);
	prof.end(5);

	ofs << "iter " << ii << " of " << max_ii_g << " " << cost << " "  << err << " " << corr_err << endl;

	system((string("/usr/local/hadoop/bin/hadoop dfs -copyFromLocal iter_output ") + get_init_string("mapred.output.dir", context)).c_str());
	prof.end(0);
      }
    else
      {
	ostringstream out, out2;
	int cnt=0;
	while(context.nextValue())
	  {
	    out << context.getInputValue() << " ";
	    cnt++;
	  }
	
	out2 << cnt << " " << out.str();
	
	context.emit(context.getInputKey(), out2.str());
      }
  }

  BPTrans_State bps;
  BPTrans_Inference bp;
};



int main(int argc, char *argv[])
{
  try 
    {
      return HadoopPipes::runTask(HadoopPipes::TemplateFactory<BPTrans_Mapper, BPTrans_Reducer>());
    } 
  catch(const string &str)
    {
      cerr << str << endl;
    }
}


#endif
