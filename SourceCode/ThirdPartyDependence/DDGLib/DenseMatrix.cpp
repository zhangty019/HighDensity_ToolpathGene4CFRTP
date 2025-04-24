#include "DenseMatrix.h"

namespace DDG
{
   template <class T>
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
       DenseMatrix<T> :: to_Eigen( void )
   {
      return *cData;
   }


   template <>
   const DenseMatrix<Real>& DenseMatrix<Real> :: operator=
       (Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B )
   // copies a cholmod_dense* into a DenseMatrix;
   // takes responsibility for deallocating B
   {
       throw std::invalid_argument("This kind of template should not used!");
      assert( B.size() );
      cData = &B;

      m = B.rows();
      n = B.cols();
      data.resize( m*n );

      for (int i = 0; i < m; i++) {
          for (int j = 0; j < n; j++) {
              data[i * n + j] = B(i, j);
          }
      }
      return *this;
   }

   template <>
   const DenseMatrix<Complex>& DenseMatrix<Complex> :: operator=
       (Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B)
   // copies a cholmod_dense* into a DenseMatrix;
   // takes responsibility for deallocating B
   {
       throw std::invalid_argument("This kind of template should not used!");
       assert(B.size());
       cData = &B;

       m = B.rows();
       n = B.cols();
       data.resize(m * n);

       for (int i = 0; i < m; i++) {
           for (int j = 0; j < n; j++) {
               data[i * n + j] = B(i, j);
           }
       }
       return *this;
   }

   template <>
   const DenseMatrix<Quaternion>& DenseMatrix<Quaternion> :: operator=
       (Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B)
   // copies a cholmod_dense* into a DenseMatrix;
   // takes responsibility for deallocating B
   {
       throw std::invalid_argument("This kind of template should not used!");
       assert(B.size());
       cData = &B;

       m = B.rows();
       n = B.cols();
       data.resize(m * n);

       for (int i = 0; i < m; i++) {
           for (int j = 0; j < n; j++) {
               data[i * n + j] = B(i, j);
           }
       }
       return *this;
   }

   template <>
   std::ostream& operator<< (std::ostream& os, const DenseMatrix<Real>& o)
   {
      const int p = 3;
      os.precision( p );
      os << scientific;

      for( int i = 0; i < o.nRows(); i++ )
      {
         os << "[ ";
         for( int j = 0; j < o.nColumns(); j++ )
         {
            double x = o(i,j);

            if( x == 0. )
            {
               os << " 0";
               for( int k = 0; k < p+6; k++ )
               {
                  os << " ";
               }
            }
            else if( x > 0. )
            {
               os << " " << x << " ";
            }
            else
            {
               os << x << " ";
            }
         }
         os << "]" << endl;
      }
   
      return os;
   }

   template <>
   std::ostream& operator<< (std::ostream& os, const DenseMatrix<Complex>& o)
   {
      const int p = 2;
      os.precision( p );
      os << scientific;

      for( int i = 0; i < o.nRows(); i++ )
      {
         os << "[ ";
         for( int j = 0; j < o.nColumns(); j++ )
         {
            Complex z = o(i,j);

            if( z.re == 0. )
            {
               os << " 0";
               for( int k = 0; k < p+5; k++ )
               {
                  os << " ";
               }
            }
            else if( z.re > 0. )
            {
               os << " " << z.re;
            }
            else
            {
               os << z.re;
            }

            if( z.im == 0 )
            {
               os << " ";
            }
            else if( z.im >= 0 )
            {
               os << "+";
            }
            else
            {
               os << "-";
            }

            if( z.im == 0. )
            {
               for( int k = 0; k < p+8; k++ )
               {
                  os << " ";
               }
            }
            else
            {
               os << abs( z.im ) << "i ";
            }
         }
         os << " ]" << endl;
      }
   
      return os;
   }

   template <>
   std::ostream& operator<< (std::ostream& os, const DenseMatrix<Quaternion>& o)
   {
      const int p = 2;
      os.precision( p );
      os << scientific;

      for( int i = 0; i < o.nRows(); i++ )
      {
         os << "[";
         for( int j = 0; j < o.nColumns(); j++ )
         {
            Quaternion q = o(i,j);

            os << " " << q;
         }
         os << " ]" << endl;
      }
   
      return os;
   }

   template <>
   void DenseMatrix<Real> :: randomize( void )
   // replaces entries with uniformly distributed real random numbers in the interval [-1,1]
   {
      for( int i = 0; i < m*n; i++ )
      {
         data[i] = 2.*unitRand() - 1.;
      }
   }

   template <>
   void DenseMatrix<Complex> :: randomize( void )
   // replaces entries with uniformly distributed real random numbers in the interval [-1,1]
   {
      for( int i = 0; i < m*n; i++ )
      {
         data[i].re = 2.*unitRand() - 1.;
         data[i].im = 2.*unitRand() - 1.;
      }
   }

   template <>
   void DenseMatrix<Quaternion> :: randomize( void )
   // replaces entries with uniformly distributed real random numbers in the interval [-1,1]
   {
      for( int i = 0; i < m*n; i++ )
      {
         for( int k = 0; k < 4; k++ )
         {
            data[i][k] = 2.*unitRand() - 1.;
         }
      }
   }





   // ÌØ»¯°æ±¾ 
   template <>
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
       fromRealDenseMatrix(DenseMatrix<Real>& A) {
       Eigen::MatrixXd M(A.nRows(), A.nColumns());
       for (int i = 0; i < A.nRows(); i++) {
           for (int j = 0; j < A.nColumns(); j++) {
               M(i, j) = (double)A(i, j);
           }
       }
       return M;
   }

   template <>
   Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
       fromComplexDenseMatrix(DenseMatrix<Complex>& A) {
       Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> M(A.nRows(), A.nColumns());
       for (int i = 0; i < A.nRows(); i++) {
           for (int j = 0; j < A.nColumns(); j++) {
               M(i, j) = std::complex<double>(A(i, j).re, A(i, j).im);
           }
       }
       return M;
   }

   template <>
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
       fromQuaternionDenseMatrix(DenseMatrix<Quaternion>& A) {
       Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M(A.nRows()*4, A.nColumns()*4);
       for (int i = 0; i < A.nRows(); i++) {
           for (int j = 0; j < A.nColumns(); j++) {
               Quaternion& q = A(i, j);
               M(i * 4 + 0, j * 4 + 0) = q[0];
               M(i * 4 + 0, j * 4 + 1) = -q[1];
               M(i * 4 + 0, j * 4 + 2) = -q[2];
               M(i * 4 + 0, j * 4 + 3) = -q[3];
               M(i * 4 + 1, j * 4 + 0) = q[1];
               M(i * 4 + 1, j * 4 + 1) = q[0];
               M(i * 4 + 1, j * 4 + 2) = -q[3];
               M(i * 4 + 1, j * 4 + 3) = q[2];
               M(i * 4 + 2, j * 4 + 0) = q[2];
               M(i * 4 + 2, j * 4 + 1) = q[3];
               M(i * 4 + 2, j * 4 + 2) = q[0];
               M(i * 4 + 2, j * 4 + 3) = -q[1];
               M(i * 4 + 3, j * 4 + 0) = q[3];
               M(i * 4 + 3, j * 4 + 1) = -q[2];
               M(i * 4 + 3, j * 4 + 2) = q[1];
               M(i * 4 + 3, j * 4 + 3) = q[0];
           }
       }
       return M;
   }

}


