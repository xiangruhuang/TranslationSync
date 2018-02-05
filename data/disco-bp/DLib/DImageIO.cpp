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
#include "stdio.h"
#include <fstream>
#include <iostream>
#include "DImage.h"
#include <assert.h>

using namespace std;

DImage ReadPPMImage(const char *filename)
{
  FILE *fp = fopen(filename, "rb");
  char temp[1024];
  assert(fp);
  fgets(temp, 1024, fp);
  assert(temp[0] == 'P' && temp[1] == '6');

  int width, height, maxval;
  do {
    fgets(temp, 1024, fp);
  } while(temp[0] == '#');

  sscanf(temp, "%d %d", &width, &height);

  fgets(temp, 1024, fp);
  sscanf(temp, "%d", &maxval);

  DImage image(3, height, width);

  for(int i=0; i<height; i++)
    for(int j=0; j<width; j++)
      for(int k=0; k<3; k++)
      {
        image[k][i][j] = fgetc(fp);
      } 

  fclose(fp);

  return image;
}


DImage ReadPBMImage(const char *filename)
{
  FILE *fp = fopen(filename, "rb");
  char temp[1024];

  fgets(temp, 1024, fp);
  assert(temp[0] == 'P' && temp[1] == '4');

  int width, height;
  do {
    fgets(temp, 1024, fp);
  } while(temp[0] == '#' || temp[0]=='\n');

  sscanf(temp, "%d %d", &width, &height);

  DImage image(1, height, width);

  for(int i=0; i<height; i++)
    for(int j=0; j<ceil(width/8.0); j++)
      {
         int tmp = fgetc(fp);
         for(int k=128, l=0; k>=1 && l<8; k=k>>1, l++)
         {
           if(j*8+l < width)
              image[0][i][j*8+l] = (tmp & k)?255:0;
         }
      } 

  fclose(fp);

  return image;
}


void WritePPMImage(const _DImage<unsigned char> &img, const char *filename)
{
  if(img.planes() == 1)
    return WritePPMImage(img[0], filename);

  assert(img.planes() == 3);

  FILE *fp = fopen(filename, "wb");

  // write magic number
  fprintf(fp, "P6\n");

  // write dimensions
  fprintf(fp, "%d %d\n", img.cols(), img.rows());
  
  // write max pixel value
  fprintf(fp, "255\n");

  for(int i=0; i<img.rows(); i++)
    for(int j=0; j<img.cols(); j++)
      for(int k=0; k<3; k++)
	fputc(int(img[k][i][j]), fp);

  fclose(fp);

  return;
}

void WritePPMImage(const _DPlane<unsigned char> &img, const char *filename)
{
  FILE *fp = fopen(filename, "wb");

  // write magic number
  fprintf(fp, "P6\n");

  // write dimensions
  fprintf(fp, "%d %d\n", img.cols(), img.rows());
  
  // write max pixel value
  fprintf(fp, "255\n");

  for(int i=0; i<img.rows(); i++)
    for(int j=0; j<img.cols(); j++)
      for(int k=0; k<3; k++)
	fputc(int(img[i][j]), fp);

  fclose(fp);

  return;
}

#ifdef IMGIO_SUPPORT
#include <corona.h>


static _DImage<unsigned char> load_helper(corona::Image *image, bool want_alpha);



_DImage<unsigned char> LoadDImage(const char *fname, bool want_alpha)
{
  corona::Image* image = corona::OpenImage(fname, corona::PF_R8G8B8A8);
  if (!image) 
    throw(std::string("couldn't open image " + std::string(fname)));

  return load_helper(image, want_alpha);
}


_DImage<unsigned char> LoadDImage(const void *buffer, int buf_size, bool want_alpha)
{
  corona::File* file = corona::CreateMemoryFile(buffer, buf_size);
  if (!file) 
    throw(std::string("couldn't create memory file "));
  corona::Image* image = corona::OpenImage(file, corona::PF_R8G8B8A8);
  if (!image) 
    throw(std::string("couldn't open image "));

  _DImage<unsigned char> img = load_helper(image, want_alpha);
  delete file;
  return img;
}


static _DImage<unsigned char> load_helper(corona::Image *image, bool want_alpha)
{
  int width  = image->getWidth();
  int height = image->getHeight();
  void* pixels = image->getPixels();

  int planes = want_alpha ? 4 : 3;

  _DImage<unsigned char> result(planes, height, width);

  // we're guaranteed that the first eight bits of every pixel is red,
  // the next eight bits is green, and so on...
  
  unsigned char *out_cp1 = result[0][0], *out_cp2 = result[1][0], *out_cp3 = result[2][0], *out_cp4;;
  unsigned char *in_cp = (unsigned char *) pixels;
  if(want_alpha) out_cp4 = result[3][0];

  int sz = height * width;
  for(int i=0; i < sz; i++)
    {
      *(out_cp1++) = *(in_cp++);
      *(out_cp2++) = *(in_cp++);
      *(out_cp3++) = *(in_cp++);
      if(want_alpha)
	*(out_cp4++) = *in_cp;

      in_cp++; // skip over alpha channel
    }

  delete image;

  return result;
}


void SaveDImage(const char *filename, const _DImage<unsigned char> &img)
{
  char *buf = new char[img.rows() * img.cols() * 3];
  unsigned char *in_cp1, *in_cp2, *in_cp3;

  if(img.planes() == 1)
    in_cp1 = in_cp2 = in_cp3 = img[0][0];
  else if(img.planes() == 3)
    in_cp1 = img[0][0], in_cp2 = img[1][0], in_cp3 = img[2][0];
  else
    throw std::string("SaveDImage only supports images with 1 or 3 planes");

  unsigned char *out_cp = (unsigned char *) buf;

  for(int i=0; i < img.rows() * img.cols(); i++)
    {
      *(out_cp++) = *(in_cp1++);
      *(out_cp++) = *(in_cp2++);
      *(out_cp++) = *(in_cp3++);
    }

  corona::Image *image = corona::CreateImage(img.cols(), img.rows(), corona::PF_R8G8B8, buf);

  corona::SaveImage(filename, corona::FF_PNG, image);

  delete image;
  delete[] buf;
}

void SaveDImage(const char *filename, const _DPlane<unsigned char> &img)
{
  SaveDImage(filename, _DImage<unsigned char>(img));
}

#endif
