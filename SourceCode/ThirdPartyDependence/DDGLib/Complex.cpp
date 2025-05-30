#include "Complex.h"
#include <iostream>
#include <cmath>
#include <corecrt_math_defines.h>
using namespace std;

namespace DDG
{
   Complex::Complex( double a, double b )
   // constructs number a+bi
   : re( a ), im( b )
   {}

   double& Complex :: operator[]( int p )
   {
      return ( &re )[ p ];
   }

   const double& Complex :: operator[]( int p ) const
   {
      return ( &re )[ p ];
   }

   void Complex::operator+=( const Complex& z )
   // add z
   {
      re += z.re;
      im += z.im;
   }

   void Complex::operator-=( const Complex& z )
   // subtract z
   {
      re -= z.re;
      im -= z.im;
   }

   void Complex::operator*=( const Complex& z )
   // Complex multiply by z
   {
      double a = re;
      double b = im;
      double c = z.re;
      double d = z.im;

      re = a*c-b*d;
      im = a*d+b*c;
   }

   void Complex::operator*=( double r )
   // scalar multiply by r
   {
      re *= r;
      im *= r;
   }

   void Complex::operator/=( double r )
   // scalar divide by r
   {
      re /= r;
      im /= r;
   }

   void Complex::operator/=( const Complex& z )
   // scalar divide by r
   {
      *this *= z.inv();
   }

   Complex Complex::operator-( void ) const
   {
      return Complex( -re, -im );
   }

   Complex Complex::bar( void ) const
   // returns Complex conjugate
   {
      return Complex( re, -im );
   }

   Complex Complex::inv( void ) const
   // returns inverse
   {
      return this->bar() / this->norm2();
   }

   double Complex::arg( void ) const
   // returns argument
   {
      return atan2( im, re );
   }

   double Complex::norm( void ) const
   // returns norm
   {
      return sqrt( re*re + im*im );
   }

   double Complex::norm2( void ) const
   // returns norm squared
   {
      return re*re + im*im;
   }

   Complex Complex::unit( void ) const
   // returns complex number with unit norm and same modulus
   {
      return *this / this->norm();
   }

   void Complex :: normalize( void )
   {
      *this /= norm();
   }

   Complex Complex::exp( void ) const
   // complex exponentiation
   {
      return ::exp( re ) * Complex( cos( im ), sin( im ));
   }

   Complex Complex::root( double k, double n ) const
   {
      double r = norm();
      double theta = arg() / k;
      double phi = theta + 2.*n*M_PI/k;
      return Complex( r*cos(phi), r*sin(phi) );
   }

   Complex operator+( const Complex& z1, const Complex& z2 )
   // binary addition
   {
      Complex z = z1;
      z += z2;
      return z;
   }

   Complex operator-( const Complex& z1, const Complex& z2 )
   // binary subtraction
   {
      Complex z = z1;
      z -= z2;
      return z;
   }

   Complex operator*( const Complex& z1, const Complex& z2 )
   // binary Complex multiplication
   {
      Complex z = z1;
      z *= z2;
      return z;
   }

   Complex operator*( const Complex& z, double r )
   // right scalar multiplication
   {
      Complex zr = z;
      zr *= r;
      return zr;
   }

   Complex operator*( double r, const Complex& z )
   // left scalar multiplication
   {
      return z*r;
   }

   Complex operator/( const Complex& z, double r )
   // scalar division
   {
      Complex zr = z;
      zr /= r;
      return zr;
   }

   Complex operator/( const Complex& z1, const Complex& z2 )
   // complex division
   {
      Complex z = z1;
      z /= z2;
      return z;
   }

   double dot( const Complex& z1, const Complex& z2 )
   {
      return z1.re*z2.re + z1.im*z2.im;
   }

   double cross( const Complex& z1, const Complex& z2 )
   {
      return z1.re*z2.im - z1.im*z2.re;
   }

   std::ostream& operator<<( std::ostream& os, const Complex& z )
   // prints components
   {
      if( z.im > 0 )
      {
         os << z.re << " + " << z.im << "i";
      }
      else if( z.im < 0 )
      {
         os << z.re << " - " << -z.im << "i";
      }
      else
      {
         os << z.re;
      }

      return os;
   }
}

