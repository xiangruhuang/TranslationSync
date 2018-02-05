#include<DistanceTransform.h>
#include <iostream>
#include <DProfile.h>
#include <DImageIO.h>

DProfile prof(10);

using namespace std;

int main(int argc, char *argv[])
{
  int trials = atoi(argv[2]);

#ifdef ASDF
  _DImage<unsigned char> img = LoadDImage(argv[1]);
  _DMatrix<double> mat;
  change_type(img.get_luma_plane(), mat);

  DistanceTransform_2D<double> dd(0.01, true);
  DistanceTransform_2D<double> dd2(0.01, bool(false));

      prof.begin(2);
            for(int j=0; j<trials; j++)
      	dd.do_transform(mat);
	    //      dd.do_transform(mat);
      prof.end(2);
      
      prof.begin(1);
      for(int j=0; j<trials; j++)
	dd2.do_transform(mat);
      //	dd2.do_transform(row.extract_row(j));
      prof.end(1);

      cout << (fabs(dd.do_transform(mat) - dd2.do_transform(mat))).max() << endl;


#else
  int sz = atoi(argv[3]);
  DistanceTransform_2D<double> dd(10, bool(atoi(argv[1])));
  DistanceTransform_2D<double> dd2(10, bool(false));
    _DPlane<double> row = _DMatrix<double>(1, sz, _DMatrix<double>::random) * 1000;

  _DPlane<double> mat(sz,sz);
  int SZ = sz*sz;
  double *cp = mat[0];

      //      for(int j=0; j<SZ; j++)
      //	cp[j] = drand48() * 1000;

      prof.begin(2);
            for(int j=0; j<trials; j++)
      	dd.do_transform(row);
	    //      dd.do_transform(mat);
      prof.end(2);
      
      prof.begin(1);
      for(int j=0; j<trials; j++)
	dd2.do_transform(row);
      //	dd2.do_transform(row.extract_row(j));
      prof.end(1);

#endif
  //  for(int j=0; j<trials; j++)
  //     cout << (fabs(dd.do_transform(row.extract_row(j)) - dd2.do_transform(row.extract_row(j)))).max() << endl;
}
