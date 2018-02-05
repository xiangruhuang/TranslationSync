//
// test_gaussian.cpp
//
// This code exercises many parts of DLib, by computing a 2-D gaussian
// probability distribution in several different ways:
//
//  1. using DGaussianModel::get_likelihood_plane(), which computes
//     the distribution directly from the Gaussian PDF
//
//  2. doing DPlane::convolve_gaussian() on an impulse function, which
//     uses image rotations and two 1-D convolutions of Gaussians
//
//  3. using a distance transform, which uses image rotations and 
//     the fast L2 distance transform algorithm of [FH05]
//
//  4. doing DPlane::cross_correlate() on an impulse function and using
//     a DGaussianKernel as the convolution kernel
// 
//  5. same as #4, but using DPlane::cross_correlate_fft()
//


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
  _DMatrix<double> sigma(2,2);

  sigma = 0;
  sigma[0][0] = atof(argv[1]);
  sigma[1][0] = sigma[0][1] = atof(argv[2]);
  sigma[1][1] = atof(argv[3]);

  cout << setw(5) << setprecision(4) << sigma << endl;

  const int sz = 101;

  _DMatrix<double> test(sz,sz);
  test = 0;
  test[sz/2][sz/2] = 1;
  _DMatrix<double> res(5,sz);

  _DMatrix<double> mean(2,1);
  mean[0][0]=mean[0][1]=sz/2;
  DGaussianModel<double> model(sigma, mean);

  cout.flags(ios::fixed | ios::right);
  //  cout << setw(7) << setprecision(2) << model.get_likelihood_plane(sz,sz).extract_row(50) << endl;
  res.set_row(0, model.get_likelihood_plane(sz,sz).extract_row(50));

  const int sigma_count = 7;

  cout.flags(ios::fixed | ios::right);
  //  cout << setw(7) << setprecision(2) << log(_DPlane<double>(test).convolve_gaussian(sigma)).extract_row(50) << endl;
  res.set_row(1, log(_DPlane<double>(test).convolve_gaussian(sigma, sigma_count)).extract_row(50));
  //  cout << setw(7) << setprecision(10) << (_DPlane<double>(test).convolve_gaussian(sigma)) << endl;

  _DMatrix<double> sigmas2 = sigma * 2.0;

  DistanceTransform_2D<double> dt(sigmas2);
  
  _DMatrix<double> test2(sz,sz);
  test2 = 1e100;
  test2[sz/2][sz/2] = 0;
  test2[sz/2][sz/2+1] = 0;
  test2[sz/2+1][sz/2] = 0;
  test2[sz/2+1][sz/2+1] = 0;

  _DMatrix<double> dt_image = dt.do_transform(test2); 
  dt_image= -dt_image - (1/2.0) * 2.0 * log(2.0*M_PI) - (1 / 2.0) * log(sigma.determinant());

  //  cout << setw(7) << setprecision(2) << dt_image.extract_row(50) << endl;
  res.set_row(2, dt_image.extract_row(50));

  res.set_row(3, log(_DPlane<double>(test).cross_correlate(_DGaussianKernel<double>(sigma, sigma_count)).extract_row(50)));

  res.set_row(4, log(_DPlane<double>(test).cross_correlate_fft(_DGaussianKernel<double>(sigma, sz, sz)).extract_row(50)));

  cout << setw(7) << setprecision(2) << res.transpose() << endl;

  _DMatrix<double> m1 = model.get_likelihood_plane(sz,sz);
  m1 = m1-m1.min();

  _DMatrix<unsigned char> ooo;
  change_type(m1/m1.max()*255, ooo);
  SaveDImage("out1.png", ooo);
  cout << (m1/m1.max()*255).extract_row(sz/2) << endl;

  m1 = pointwise_max(log(_DPlane<double>(test).convolve_gaussian(sigma, sigma_count)), _DMatrix<double>(sz,sz,-50));
  m1 = m1-m1.min();

  change_type(m1/m1.max()*255, ooo);
  SaveDImage("out2.png", ooo);
  cout << (m1/m1.max()*255).extract_row(sz/2) << endl;

  m1 = dt_image;
  m1 = m1-m1.min();

  change_type(m1/m1.max()*255, ooo);
  SaveDImage("out3.png", ooo);
  cout << (m1/m1.max()*255).extract_row(sz/2) << endl;
} catch(string &str)
{ cout << str << endl;
}
}
