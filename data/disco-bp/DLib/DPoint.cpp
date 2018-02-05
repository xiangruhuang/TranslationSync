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
#include "DPoint.h"
#include <iostream>

using namespace std;

template<class T>
_DPoint<T> operator+(const _DPoint<T> &p1, const _DPoint<T> &p2)
{
  return _DPoint<T>(p1.row() + p2.row(), p1.col() + p2.col());
}

template<class T>
_DPoint<T> operator-(const _DPoint<T> &p1, const _DPoint<T> &p2)
{
  return _DPoint<T>(p1.row() - p2.row(), p1.col() - p2.col());
}

template<class T>
_DPoint<T> operator-(const _DPoint<T> &p1, double p2)
{
  return _DPoint<T>(T(p1.row() - p2), T(p1.col() - p2));
}

template<class T>
_DPoint<T> operator*(const _DPoint<T> &p1, const _DPoint<T> &p2)
{
  return _DPoint<T>(p1.row() * p2.row(), p1.col() * p2.col());
}

template<class T>
_DPoint<T> operator*(const _DPoint<T> &p1, double factor)
{
  return _DPoint<T>(T(p1.row() * factor), T(p1.col() * factor));
}

template<class T>
_DPoint<T> operator/(const _DPoint<T> &p1, double p2)
{
  return _DPoint<T>(T(p1.row() / p2), T(p1.col() / p2));
}

template<class T>
_DPoint<T> operator/(const _DPoint<T> &p1, const _DPoint<T> &p2)
{
  return _DPoint<T>(T(p1.row() / p2.row()), T(p1.col() / p2.col()));
}

template<class T>
bool operator<(const _DPoint<T> &p1, const _DPoint<T> &p2)
{
  return p1.row() < p2.row() || (p1.row() == p2.row() && p1.col() < p2.col());
}

template<class T>
bool operator>(const _DPoint<T> &p1, const _DPoint<T> &p2)
{
  return p1.row() > p2.row() || (p1.row() == p2.row() && p1.col() > p2.col());
}

template<class T>
bool operator<=(const _DPoint<T> &p1, const _DPoint<T> &p2)
{
  return p1 < p2 || p1 == p2;
}

template<class T>
bool operator>=(const _DPoint<T> &p1, const _DPoint<T> &p2)
{
  return p1 > p2 || p1 == p2;
}

template<class T>
bool operator==(const _DPoint<T> &p1, const _DPoint<T> &p2)
{
  return p1.row() == p2.row() && p1.col() == p2.col();
}

template<class T>
bool operator!=(const _DPoint<T> &p1, const _DPoint<T> &p2)
{
  return p1.row() != p2.row() || p1.col() != p2.col();
}

template<class T>
ostream &operator<<(ostream &os, const _DPoint<T> &p)
{
  os << "(" << p.row() << ", " << p.col() << ")";
  return os;
}

template<class T>
_DPoint<T> sqrt(const _DPoint<T> &p)
{
  return _DPoint<T>((T)sqrt((double) p.row()), (T)sqrt((double) p.col()));
}

template<class T>
_DPoint<T> exp(const _DPoint<T> &p)
{
  return _DPoint<T>((T)exp((double) p.row()), (T)exp((double) p.col()));
}

template<class T>
_DPoint<T> log(const _DPoint<T> &p)
{
  return _DPoint<T>((T)log((double) p.row()), (T)log((double) p.col()));
}

template<class T>
double atan(const _DPoint<T> &p)
{
  return 0; //_DPoint<T>((int)atan(p.row()), (int)atan(p.col()));
}

template<class T>
std::istream &operator>>(std::istream &is, _DPoint<T> &p){
  
  char c;

  do { is >> c; } while(c == ' ');

  assert(c=='(');

  is >> p._row;

  do { is >> c; } while(c == ' ');

  assert(c==',');

  is >> p._col;

  do { is >> c; } while(c == ' ');

  assert(c==')');

  return is;
}


#define DECLARE(x) \
  template _DPoint<x> operator+(const _DPoint<x> &p1, const _DPoint<x> &p2); \
  template _DPoint<x> operator-(const _DPoint<x> &p1, const _DPoint<x> &p2); \
  template _DPoint<x> operator-(const _DPoint<x> &p1, double p2); \
  template _DPoint<x> operator*(const _DPoint<x> &p1, const _DPoint<x> &p2); \
  template _DPoint<x> operator*(const _DPoint<x> &p1, double factor); \
  template bool operator<(const _DPoint<x> &p1, const _DPoint<x> &p2); \
  template bool operator>(const _DPoint<x> &p1, const _DPoint<x> &p2); \
  template bool operator<=(const _DPoint<x> &p1, const _DPoint<x> &p2); \
  template bool operator>=(const _DPoint<x> &p1, const _DPoint<x> &p2); \
  template bool operator==(const _DPoint<x> &p1, const _DPoint<x> &p2); \
  template bool operator!=(const _DPoint<x> &p1, const _DPoint<x> &p2); \
  template _DPoint<x> operator/(const _DPoint<x> &p1, double p2); \
  template _DPoint<x> operator/(const _DPoint<x> &p1, const _DPoint<x> &p2); \
  template _DPoint<x> sqrt(const _DPoint<x> &p); \
  template _DPoint<x> log(const _DPoint<x> &p); \
  template _DPoint<x> exp(const _DPoint<x> &p); \
  template std::ostream &operator<<(std::ostream &os, const _DPoint<x> &p); \
  template std::istream &operator>>(std::istream &is, _DPoint<x> &p);

DECLARE(double)
DECLARE(int)
DECLARE(float)



