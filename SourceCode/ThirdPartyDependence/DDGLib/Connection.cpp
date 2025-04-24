#include "Connection.h"
#include <cmath>
#include <iostream>
#include "qdebug.h"
//#include <Eigen/PardisoSupport>

using namespace std;
//using namespace cm;

namespace DDG
{
	Connection::Connection(Mesh& _mesh)
		: mesh(_mesh)
	{
		build();
	}

	Connection :: ~Connection(void)
	{
		//if( QR != NULL ) SuiteSparseQR_free<double>( &QR, common );
	}

	void Connection::build(void)
		// build the system of constraint equations
	{
		// construct a list of basis cycles, represented by a vector
		// of oriented edges (really: halfedges)
		std::vector<Cycle> basisCycles;

		// append contractible basis cycles
		append1RingBases(basisCycles);

		unsigned int nContractibleCycles = basisCycles.size();

		// append noncontractible basis cycles
		mesh.appendDualGenerators(basisCycles);

		nBasisCycles = basisCycles.size();

		// append derivative constraints used to specify directional constraints in faces
		vector<double> directionalHolonomies;
		mesh.appendDirectionalConstraints(basisCycles, directionalHolonomies);

		// build constraint matrix
		SparseMatrix<Real> A(mesh.edges.size(), basisCycles.size());
		buildCycleMatrix(A, basisCycles);
		applyCotanWeights(A);

		// prefactor [Q,R,E] = qr( A' )

		//if( QR != NULL ) { SuiteSparseQR_free<double>( &QR, common ); QR = NULL; }
		Eigen::SparseMatrix<double, Eigen::ColMajor> _A = fromRealSparseMatrix<Real>(A.transpose());
		_A.makeCompressed();
		
		qDebug() << "Row\tCol\tVal" << endl;
		for (int k = 0; k < _A.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(_A, k); it; ++it)
			{
				//qDebug() << it.row() << "\t"; // row index
				//qDebug() << it.col() << "\t"; // col index (here it is equal to k)
				//qDebug() << it.value() << endl;
				qDebug() << it.row() << " " << it.col() << " " << it.value();
			}
		}


		qDebug() << "start QR";
		QR.compute(_A);
		qDebug() << "end QR";
		//QR = SuiteSparseQR_factorize <double> (7, -2., *A, common);

		// compute Riemannian holonomy of basis cycles
		K = DenseMatrix<Real>(basisCycles.size(), 1);
		b = DenseMatrix<Real>(basisCycles.size(), 1);
		generatorOnBoundary.resize(mesh.nGenerators());
		for (unsigned int i = 0; i < nContractibleCycles; i++)
		{
			K(i) = -mesh.vertices[i].defect();
		}
		for (unsigned int i = nContractibleCycles;
			i < nContractibleCycles + mesh.nGenerators();
			i++)
		{
			if (Mesh::isDualBoundaryLoop(basisCycles[i]))
			{
				K(i) = -mesh.boundaryLoopCurvature(basisCycles[i]);
				generatorOnBoundary[i - nContractibleCycles] = true;
			}
			else
			{
				K(i) = -mesh.defect(basisCycles[i]);
				generatorOnBoundary[i - nContractibleCycles] = false;
			}
		}

		// specify change in angle for each directional constraint
		for (unsigned int i = 0; i < directionalHolonomies.size(); i++)
		{
			K(i + nBasisCycles) = -directionalHolonomies[i];
		}

		// setup the right hand side using the Riemannian holonomy, ignoring
		// singularities for now; also make sure rhsChanged() is initially true
		for (int i = 0; i < b.length(); i++)
		{
			b(i) = K(i);
		}
	}

	void Connection::setupRHS(void)
		// add 2*pi*k to the right hand side, where k is the vector of singularity/generator indices
	{
		double indexSum = 0;

		// iterate over vertices
		for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
		{
			int i = vertex2row[v->index]; // get the row index of the current vertex

			b(i) = K(i) + 2 * M_PI * v->k;

			indexSum += v->k;
		}

		// iterate over generators
		int nContractibleCycles = nBasisCycles - mesh.generatorIndices.size();
		for (int i = 0; i < (int)mesh.generatorIndices.size(); i++)
		{
			int j = nContractibleCycles + i;
			b(j) = K(j) + 2. * M_PI * mesh.generatorIndices[i];

			if (generatorOnBoundary[i])
			{
				indexSum += mesh.generatorIndices[i];
			}
		}

		// display a warning if the sum of singular indices does not
		// add up to the Euler characteristic
		if (abs(indexSum - mesh.eulerCharacteristic()) > 1e-7)
		{
			/*
			cerr << endl;
			cerr << "  *************************************************************" << endl;
			cerr << "    Warning: indices do not add up to Euler characteristic!" << endl;
			cerr << "             (solution may have unwanted singularities)" << endl;
			cerr << endl;
			cerr << "             Euler characteristic: " << mesh.eulerCharacteristic() << endl;
			cerr << "                   sum of indices: " << indexSum << endl;
			cerr << "  *************************************************************" << endl;
			cerr << endl;
			*/
			qDebug() << endl;
			qDebug() << "  *************************************************************" << endl;
			qDebug() << "    Warning: indices do not add up to Euler characteristic!" << endl;
			qDebug() << "             (solution may have unwanted singularities)" << endl;
			qDebug() << endl;
			qDebug() << "             Euler characteristic: " << mesh.eulerCharacteristic() << endl;
			qDebug() << "                   sum of indices: " << indexSum << endl;
			qDebug() << "  *************************************************************" << endl;
			qDebug() << endl;
		}
	}

	void Connection::resetRHS(void)
	{
		// make a copy of the right hand side used in the most recent solve
		for (int i = 0; i < b.length(); i++)
		{
			b(i) = K(i);
		}
	}

	bool Connection::update(void)
	{
		// specify singular values in right hand side b
		setupRHS();

		// Eigen::MatrixXd _b = fromRealDenseMatrix(b);
		// Eigen::MatrixXd _y = QR.solve(_b);
		// Eigen::MatrixXd _x = QR.matrixQ() * _y;
		// solve y = R'\(E'*b)
		// y = SuiteSparseQR_solve (SPQR_RTX_EQUALS_ETB, A, *b, common) ;

		// compute w = Q*y
		//x = SuiteSparseQR_qmult (SPQR_QX, QR, *y, common) ;

		Eigen::MatrixXd _b = fromRealDenseMatrix(b);
		Eigen::MatrixXd _x = QR.solve(_b);
		DenseMatrix<Real> x(_x.rows(), _x.cols());
		for (int i = 0; i < _x.rows(); i++) {
			for (int j = 0; j < _x.cols(); j++) {
				x(i, j) = Real(_x(i, j));
			}
		}
		applyCotanWeights(x);

		for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
		{
			e->theta = x(e->index);
		}

		// restore original right hand side
		resetRHS();

		return true;
	}

	void Connection::applyCotanWeights(SparseMatrix<Real>& A)
	{
		for (SparseMatrix<Real>::iterator e = A.begin(); e != A.end(); e++)
		{
			int row = e->first.second;
			//double& val(e->second);
			double s = mesh.edges[row].star;
			e->second *= sqrt(s);
			//val *= sqrt( s );
		}
	}

	void Connection::applyCotanWeights(DenseMatrix<Real>& x)
	{
		for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
		{
			x(e->index) *= sqrt(e->star);
		}
	}

	void Connection::buildCycleMatrix(SparseMatrix<Real>& A, vector<Cycle>& cycles) const
	{
		for (unsigned int l = 0; l < cycles.size(); l++)
		{
			for (Cycle::iterator h = cycles[l].begin();
				h != cycles[l].end();
				h++)
			{
				int k = (*h)->edge->index;
				int i = (*h)->vertex->index;
				int j = (*h)->flip->vertex->index;

				if (i > j) A(k, l) = -1.;
				else        A(k, l) = 1.;
			}
		}
	}

	void Connection::append1RingBases(vector<Cycle>& cycles)
	{
		// contractible bases
		vertex2row.resize(mesh.vertices.size());
		for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
		{
			if (v->onBoundary())
			{
				vertex2row[v->index] = -1;
				continue;
			}

			Cycle c;
			HalfEdgeIter he = v->he;
			do
			{
				c.push_back(he);
				he = he->flip->next;
			} while (he != v->he);

			vertex2row[v->index] = cycles.size();
			cycles.push_back(c);
		}
	}
}

