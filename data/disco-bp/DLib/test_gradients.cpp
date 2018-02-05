#include <DMatrix.h>
#include <DGaussianModel.h>
#include <DPlane.h>
#include <DistanceTransform.h>
#include <iomanip>
#include <DImageIO.h>
#include <DGaussianKernel.h>

using namespace std;

main(int argc, char *argv[])
{
try {
  _DImage<unsigned char> img = LoadDImage(argv[1]);
    _DPlane<double> plane;
    change_type(img.get_luma_plane(), plane);

  // compute x gradients
  _DPlane<double> x_grad = plane.get_x_gradient();
  _DPlane<double> y_grad = plane.get_y_gradient();
  _DPlane<double> mag = sqrt(sqr(x_grad) + sqr(y_grad));

  DMatrix cov(2,2);
  cov=0;
  cov[0][0] = 10;
  cov[1][1] = 5;
  _DPlane<double> conv_result = plane.convolve_gaussian(cov, 20);

  _DPlane<unsigned char> result;
  change_type(conv_result, result);
  SaveDImage("mag.png", result);

} catch(string &str)
{ cout << str << endl;
}
}
