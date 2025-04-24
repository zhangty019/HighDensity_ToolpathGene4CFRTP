#include "SparseMatrix.h"
#include "DenseMatrix.h"

namespace DDG
{
   

    //special version will used
    template <>
    Eigen::SparseMatrix<double>
        fromRealSparseMatrix(const SparseMatrix<Real>& A) {
        std::vector<Eigen::Triplet<double>> triplets;
        for (SparseMatrix<Real>::const_iterator e = A.begin();
            e != A.end();
            e++)
        {
            int i = e->first.second;
            int j = e->first.first;
            const Real& q(e->second);
            triplets.push_back(Eigen::Triplet<double>(i, j, (double)q));
        }
        Eigen::SparseMatrix<double> _M(A.nRows(), A.nColumns());
        _M.setFromTriplets(triplets.begin(), triplets.end());
        return _M;
    }
    
    template <>
    Eigen::SparseMatrix<std::complex<double>>
        fromComplexSparseMatrix(const SparseMatrix<Complex>& A) {
        std::vector<Eigen::Triplet<std::complex<double>>> triplets;
        for (SparseMatrix<Complex>::const_iterator e = A.begin();
            e != A.end();
            e++)
        {
            int i = e->first.second;
            int j = e->first.first;
            const Complex& q(e->second);
            std::complex<double> _q(q.re, q.im);
            triplets.push_back(Eigen::Triplet<std::complex<double>>(i, j, _q));
        }

        Eigen::SparseMatrix<std::complex<double>>
            _M(A.nRows(), A.nColumns());
        _M.setFromTriplets(triplets.begin(), triplets.end());
        return _M;
    }
    
    template <>
    Eigen::SparseMatrix<double>
        fromQuaternionSparseMatrix(const SparseMatrix<Quaternion>& A) {
        std::vector<Eigen::Triplet<double>> triplets;
        for (SparseMatrix<Quaternion>::const_iterator e = A.begin();
            e != A.end();
            e++)
        {
            int i = e->first.second;
            int j = e->first.first;
            const Quaternion& q(e->second);
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 0, j * 4 + 0, q[0]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 0, j * 4 + 1, -q[1]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 0, j * 4 + 2, -q[2]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 0, j * 4 + 3, -q[3]));

            triplets.push_back(Eigen::Triplet<double>(i * 4 + 1, j * 4 + 0, q[1]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 1, j * 4 + 1, q[0]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 1, j * 4 + 2, -q[3]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 1, j * 4 + 3, q[2]));

            triplets.push_back(Eigen::Triplet<double>(i * 4 + 2, j * 4 + 0, q[2]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 2, j * 4 + 1, q[3]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 2, j * 4 + 2, q[0]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 2, j * 4 + 3, -q[1]));

            triplets.push_back(Eigen::Triplet<double>(i * 4 + 3, j * 4 + 0, q[3]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 3, j * 4 + 1, -q[2]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 3, j * 4 + 2, q[1]));
            triplets.push_back(Eigen::Triplet<double>(i * 4 + 3, j * 4 + 3, q[0]));
        }

        Eigen::SparseMatrix<double> _M(A.nRows() * 4, A.nColumns() * 4);
        _M.setFromTriplets(triplets.begin(), triplets.end());
        return _M;
    }


   template <>
   const SparseMatrix<Real>& SparseMatrix<Real> :: 
       operator=( Eigen::SparseMatrix<double> B )
   {
      assert( B.size());
      cData = &B;
      m = B.rows();
      n = B.cols();
      resize( m, n );

      for (int k = 0; k < B.outerSize(); ++k) {
          for (Eigen::SparseMatrix<double>::InnerIterator it(B, k); it; ++it)
          {
              (*this) (it.row(), it.col()) = it.value();
          }
      }
      return *this;
   }

   template <>
   const SparseMatrix<Complex>& SparseMatrix<Complex> :: 
       operator=(Eigen::SparseMatrix<double> B )
   {
      assert(B.size());
      cData = &B;

      m = B.rows();
      n = B.cols();
      resize(m, n);

      for (int k = 0; k < B.outerSize(); ++k) {
          for (Eigen::SparseMatrix<double>::InnerIterator it(B, k); it; ++it)
          {
              (*this) (it.row(), it.col()) = it.value();
          }
      }
      return *this;
   }


   void SparseMatrix<Quaternion> ::allocateSparse(void) {
       //do not need do anything for eigen
   }


   template <>
   void SparseMatrix<Real> :: setEntry( const_iterator e, int i, double* pr )
   {
      pr[i] = e->second;
   }

   template <>
   void SparseMatrix<Complex> :: setEntry( const_iterator e, int i, double* pr )
   {
      pr[i*2+0] = e->second.re;
      pr[i*2+1] = e->second.im;
   }

   template <>
   void solve( SparseMatrix<Real>& A,
                DenseMatrix<Real>& x,
                DenseMatrix<Real>& b )
   // solves the sparse linear system Ax = b using sparse QR factorization
   {
#ifdef SP_DEBUG
      int t0 = clock();
#endif
    //wrapper
      auto MA = fromRealSparseMatrix(A);
      auto Mb = fromRealDenseMatrix(b);
      Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> _SparseQR(MA);
      auto _x = _SparseQR.solve(Mb);
      x = DenseMatrix<Real>(x.nRows(), x.nColumns());
      for (int i = 0; i < x.nColumns(); i++) {
          for (int j = 0; j < x.nRows(); j++) {
              x(i, j) = Real(_x(i, j));
          }
      }
#ifdef SP_DEBUG
      int t1 = clock();
      cout << "[qr] time: " << seconds( t0, t1 ) << "s" << "\n";
      cout << "[qr] max residual: " << residual( A, x, b ) << "\n";
      cout << "[qr] size: " << A.nRows() << " x " << A.nColumns() << "\n";
      cout << "[qr] rank: " << (*context).SPQR_istat[4] << "\n";
#endif
   }

   template <>
   void solve( SparseMatrix<Complex>& A,
                DenseMatrix<Complex>& x,
                DenseMatrix<Complex>& b )
   // solves the sparse linear system Ax = b using sparse QR factorization
   {
#ifdef SP_DEBUG
      int t0 = clock();
#endif
      auto MA = fromComplexSparseMatrix(A);
      auto Mb = fromComplexDenseMatrix(b);
      Eigen::SparseQR<Eigen::SparseMatrix<std::complex<double>>, Eigen::COLAMDOrdering<int>> _SparseQR(MA);
      auto _x = _SparseQR.solve(Mb);
      x = DenseMatrix<Complex>(x.nRows(), x.nColumns());
      for (int i = 0; i < x.nColumns(); i++) {
          for (int j = 0; j < x.nRows(); j++) {
              x(i, j) = Complex(_x(i, j).real(), _x(i, j).imag());
          }
      }
#ifdef SP_DEBUG
      int t1 = clock();
      cout << "[qr] time: " << seconds( t0, t1 ) << "s" << "\n";
      cout << "[qr] max residual: " << residual( A, x, b ) << "\n";
      cout << "[qr] size: " << A.nRows() << " x " << A.nColumns() << " (complex)" << "\n";
      cout << "[qr] rank: " << (*context).SPQR_istat[4] << "\n";
#endif
   }

   template <>
   void solve( SparseMatrix<Quaternion>& A,
                DenseMatrix<Quaternion>& x,
                DenseMatrix<Quaternion>& b )
   // solves the sparse linear system Ax = b using sparse QR factorization
   {
#ifdef SP_DEBUG
      int t0 = clock();
#endif
      auto MA = fromQuaternionSparseMatrix(A);
      auto Mb = fromQuaternionDenseMatrix(b);
      Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> _SparseQR(MA);
      auto _x = _SparseQR.solve(Mb);
      x = DenseMatrix<Quaternion>(x.nRows(), x.nColumns()/4);
      for (int i = 0; i < x.nColumns(); i++) {
          for (int j = 0; j < x.nRows(); j++) {
              x(i, j) = Quaternion(_x(i, j*4), _x(i, j * 4+1), 
                  _x(i, j * 4+2), _x(i, j * 4+3));
          }
      }
#ifdef SP_DEBUG
      int t1 = clock();
      cout << "[qr] time: " << seconds( t0, t1 ) << "s" << "\n";
      cout << "[qr] max residual: " << residual( A, x, b ) << "\n";
      cout << "[qr] size: " << A.nRows() << " x " << A.nColumns() << " (quaternion)" << "\n";
      cout << "[qr] rank: " << (*context).SPQR_istat[4]/4 << "\n";
#endif
   }

   template <>
   void solveSymmetric( SparseMatrix<Complex>& A,
                        DenseMatrix<Complex>& x,
                        DenseMatrix<Complex>& b )
   // solves the sparse linear system Ax = b using sparse LU factorization
   {
#ifdef SP_DEBUG
      int t0 = clock();
#endif
      auto MA = fromComplexSparseMatrix(A);
      auto Mb = fromComplexDenseMatrix(b);
      Eigen::SparseQR<Eigen::SparseMatrix<std::complex<double>>, 
          Eigen::COLAMDOrdering<int>> _SparseQR(MA);
      auto _x = _SparseQR.solve(Mb);
      x = DenseMatrix<Complex>(x.nColumns(), x.nRows());
      for (int i = 0; i < x.nColumns(); i++) {
          for (int j = 0; j < x.nRows(); j++) {
              x(i, j) = Complex(_x(i, j).real(), _x(i, j).imag());
          }
      }
#ifdef SP_DEBUG
      int t1 = clock();
      cout << "[lu] time: " << seconds( t0, t1 ) << "s" << "\n";
      cout << "[lu] max residual: " << residual( A, x, b ) << "\n";
#endif
   }

   
   template <>
   void smallestEigPositiveDefinite(SparseMatrix<Real>& A,
       SparseMatrix<Real>& B,
       DenseMatrix<Real>& x)
       // solves A x = lambda x for the smallest nonzero eigenvalue lambda
       // A must be positive (semi-)definite; x is used as an initial guess
   {
#ifdef SP_DEBUG
       int t0 = clock();
#endif
       auto MA = fromRealSparseMatrix(A);
       auto MB = fromRealSparseMatrix(B);

       Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol;
       chol.compute(MA);

       x.randomize();
       auto _x = fromRealDenseMatrix(x);

       for (int iter = 0; iter < maxEigIter; iter++)
       {
           _x = MB * _x;

           //backsolvePositiveDefinite(L, x, x);
           _x = chol.solve(_x);
           
           _x /= sqrt(((MB * _x).transpose() * (_x)).norm());
           //x /= sqrt(inner(B * x, x).norm());
       }

       for (int i = 0; i < x.nRows(); i++) {
           for (int j = 0; j < x.nColumns(); j++) {
               x(i, j) = Real(_x(i, j));
           }
       }
#ifdef SP_DEBUG
       int t1 = clock();
       cout << "[eig] time: " << seconds(t0, t1) << "s" << "\n";
       cout << "[eig] max residual: " << residual(A, B, x) << "\n";
       cout << "[eig] energy: " << inner(A * x, x) << "\n";
#endif
   }
   
   template <>
   void smallestEigPositiveDefinite(SparseMatrix<Complex>& A,
       SparseMatrix<Complex>& B,
       DenseMatrix<Complex>& x)
       // solves A x = lambda x for the smallest nonzero eigenvalue lambda
       // A must be positive (semi-)definite; x is used as an initial guess
   {
#ifdef SP_DEBUG
       int t0 = clock();
#endif
       auto MA = fromComplexSparseMatrix(A);
       auto MB = fromComplexSparseMatrix(B);

       Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> chol;
       chol.compute(MA);

       x.randomize();
       auto _x = fromComplexDenseMatrix(x);

       for (int iter = 0; iter < maxEigIter; iter++)
       {
           _x = MB * _x;

           //backsolvePositiveDefinite(L, x, x);
           _x = chol.solve(_x);

           _x /= sqrt(((MB * _x).conjugate().transpose() * (_x)).norm());
           //x /= sqrt(inner(B * x, x).norm());
       }

       for (int i = 0; i < x.nRows(); i++) {
           for (int j = 0; j < x.nColumns(); j++) {
               x(i, j) = Complex(_x(i, j).real(), _x(i,j).imag());
           }
       }
#ifdef SP_DEBUG
       int t1 = clock();
       cout << "[eig] time: " << seconds(t0, t1) << "s" << "\n";
       cout << "[eig] max residual: " << residual(A, B, x) << "\n";
       cout << "[eig] energy: " << inner(A * x, x) << "\n";
#endif
   }

   template <class T>
   void smallestEigPositiveDefinite(SparseMatrix<T>& A,
       SparseMatrix<T>& B,
       DenseMatrix<T>& E,
       DenseMatrix<T>& x)
       // solves A x = lambda B x for the smallest nonzero eigenvalue lambda orthogonal
       // to the vectors spanned by the columns of E, which must be orthonormal. A must
       // be positive (semi-)definite, B must be symmetric; x is used as an initial guess
   {
#ifdef SP_DEBUG
       int t0 = clock();
#endif
       throw std::invalid_argument("This kind of template should not used!");
       SparseFactor<T> L;
       L.build(A);

       x.randomize();

       for (int iter = 0; iter < maxEigIter; iter++)
       {
           x = B * x;
           backsolvePositiveDefinite(L, x, x);

           // project out components of E
           for (int k = 0; k < E.nColumns(); k++)
           {
               double xDotEk = 0.;
               for (int i = 0; i < x.nRows(); i++)
               {
                   xDotEk += x(i) * E(i, k);
               }
               for (int i = 0; i < x.nRows(); i++)
               {
                   x(i) -= xDotEk * E(i, k);
               }
           }

           x /= sqrt(inner(B * x, x).norm());
       }
#ifdef SP_DEBUG
       int t1 = clock();
       cout << "[eig] time: " << seconds(t0, t1) << "s" << "\n";
       cout << "[eig] max residual: " << residual(A, B, x) << "\n";
       cout << "[eig] energy: " << inner(A * x, x) << "\n";
#endif
   }

   template <class T>
   void nSmallestEigsPositiveDefinite(SparseMatrix<T>& A,
       SparseMatrix<T>& B,
       std::vector<DenseMatrix<T> >& V,
       std::vector<double>& D,
       int nEigs)
       // computes the n smallest eigenvectors and eigenvalues of A with
       // respect to the inner product B; V stores the eigenvectors as
       // columns and D is a list of eigenvalues (in the same order); n
       // is the requested number of eigenpairs
   {
       throw std::invalid_argument("This kind of template should not used!");
       int nRows = A.nRows();
       V.resize(nEigs);
       D.resize(nEigs);

       SparseFactor<T> L;
       L.build(A);

       for (int k = 0; k < nEigs; k++)
       {
           V[k] = DenseMatrix<T>(nRows);
           V[k].randomize();

           for (int iter = 0; iter < maxEigIter; iter++)
           {
               V[k] = B * V[k];
               backsolvePositiveDefinite(L, V[k], V[k]);

               // project out previous vectors
               for (int j = 0; j < k; j++)
               {
                   Real Cjk = inner(V[j], V[k]);
                   V[k] -= Cjk * V[j];
               }

               V[k] /= sqrt(inner(B * V[k], V[k]).norm());
           }

           D[k] = inner(A * V[k], V[k]);
       }
   }

}
