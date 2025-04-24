//============================================================
// Connection.h
// 
// A Connection is used to compute a trivial connection (with
// specified singularities) on a Mesh.
//

#ifndef TCODS_CONNECTION_H
#define TCODS_CONNECTION_H

#include "HalfEdge.h"
#include "Mesh.h"
//#include "CMWrapper.h"
//#include <SuiteSparseQR.hpp>
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "Eigen/Core"
#include <vector>
#include <set>

namespace DDG
{
    class Connection
    {
    public:
        Connection(Mesh& mesh);
        ~Connection(void);
        bool update(void); // recompute connection; returns true iff connection changed

    protected:
        void build(void);                                                             // build and factorize constraint matrix
        void buildCycleMatrix(SparseMatrix<Real>& A, std::vector<Cycle>& cycles) const; // builds matrix A where each row encodes a cycle in "cycles"
        void append1RingBases(std::vector<Cycle>& cycles);                      // constructs a cycle around the dual cell of every vertex
        void setupRHS(void);                                                          // modify right hand side to take singularities into account
        void resetRHS(void);                                                          // reset right hand side to nonsingular values
        void applyCotanWeights(SparseMatrix<Real>& A);                                        // replace A with sqrt(star1)*A, where star1 is the Hodge star on primal 1-forms
        void applyCotanWeights(DenseMatrix<Real>& x);                                         // replace x with sqrt(star1)*x

        Mesh& mesh;

        //cm::Common common;                       // common interface to CHOLMOD/SPQR
        DenseMatrix<Real> K;                             // vector containing Riemannian holonomy around basis cycles
        DenseMatrix<Real> b;                             // vector encoding right hand side of linear system
        //SuiteSparseQR_factorization<double>* QR; // QR factorization of constraint matrix
        Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > QR;

        std::vector<int> vertex2row;             // maps vertex indices to matrix row indices
        unsigned int nBasisCycles;               // number of basis cycles (excludes directional constraints)
        std::vector<bool> generatorOnBoundary;   // indicates whether each generator is part of the surface boundary

    };
};

#endif
