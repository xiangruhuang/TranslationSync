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
#include "DRect.h"

template<class T>
_DRect<T> bounding_rectangle(const _DRect<T> &r1, const _DRect<T> &r2)
{
  return _DRect<T>(std::min(r1.top(), r2.top()),
	       std::min(r1.left(), r2.left()),
	       std::max(r1.bottom(), r2.bottom()),
	       std::max(r1.right(), r2.right()));
}

template<class T>
_DRect<T> intersection(const _DRect<T> &r1, const _DRect<T> &r2)
{
  return _DRect<T>(std::max(r1.top(), r2.top()),
	       std::max(r1.left(), r2.left()),
	       std::min(r1.bottom(), r2.bottom()),
	       std::min(r1.right(), r2.right()));
}

template<class T>
T area_of_intersection(const _DRect<T> &r1, const _DRect<T> &r2)
{
  T area = intersection(r1, r2).area();

  if(area < 0) 
    area = 0;

  return area;
}

template<class T>
T area_of_union(const _DRect<T> &r1, const _DRect<T> &r2)
{
  return r1.area() + r2.area() - area_of_intersection(r1, r2);
}

template<class T>
std::ostream & operator<<(std::ostream &os, const _DRect<T> &rect)
{
  os << "[" << rect.top_left() << ", " << rect.bottom_right() << "]";
  return os;
}

template<class T>
std::istream & operator>>(std::istream &is, _DRect<T> &rect)
{
  char c;

  do { is >> c;  } while(c == ' ');

  assert(c=='[');

  is >> rect._top_left;

  do { is >> c; } while(c == ' ');

  assert(c==',');

  is >> rect._bottom_right;

  do { is >> c; } while(c == ' ');

  assert(c==']');
  
  return is;
}



#define DECLARE(x) \
  template  _DRect<x> bounding_rectangle(const _DRect<x> &r1, const _DRect<x> &r2); \
  template  _DRect<x> intersection(const _DRect<x> &r1, const _DRect<x> &r2); \
  template  x area_of_union(const _DRect<x> &r1, const _DRect<x> &r2); \
  template   x area_of_intersection(const _DRect<x> &r1, const _DRect<x> &r2); \
  template  std::ostream & operator<<(std::ostream &os, const _DRect<x> &rect); \
  template std::istream & operator>>(std::istream &is, _DRect<x> &rect);

DECLARE(double)
DECLARE(int)
DECLARE(float)

