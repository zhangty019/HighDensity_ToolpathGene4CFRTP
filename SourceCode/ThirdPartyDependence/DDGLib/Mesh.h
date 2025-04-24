// -----------------------------------------------------------------------------
// libDDG -- Mesh.h
// -----------------------------------------------------------------------------
//
// Mesh represents a polygonal surface mesh using the halfedge data structure.
// It is essentially a large collection of disjoint vertices, edges, and faces
// that are ``glued together'' by halfedges which encode connectivity (see
// the documentation for an illustration).  By construction, the halfedge data
// structure cannot represent nonorientable surfaces or meshes with nonmanifold
// edges.
//
// Mesh elements are referenced using iterators -- common usage of these
// iterators is to either traverse an entire vector of mesh elements:
//
//    // visit all vertices
//    for( VertexIter i = vertices.begin(); i != vertices.end(); i++ )
//    {
//       //...
//    }
//
// or to perform a local traversal over the neighborhood of some mesh element:
//
//    // visit both halfedges of edge e
//    HalfEdgeIter he = e->he;
//    do
//    {
//       // ...
//
//       he = he->flip;
//    }
//    while( he != e->he );
//
// (See Types.h for an explicit definition of iterator types.)
//
// Meshes with boundary are handled by creating an additional face for each
// boundary loop (the method Face::isBoundary() determines whether a given
// face is a boundary loop).  Isolated vertices (i.e., vertiecs not contained
// in any edge or face) reference a dummy halfedge and can be checked via
// the method Vertex::isIsolated().
//

#ifndef DDG_MESH_H
#define DDG_MESH_H


#include <vector>
#include <string>

#include "HalfEdge.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "Complex.h"

//#include "Connection.h"

namespace DDG
{
	class Connection;

	class Mesh
	{
	public:

		Mesh(void);
		// constructs an empty mesh

		Mesh(const Mesh& mesh);
		// constructs a copy of mesh

		const Mesh& operator=(const Mesh& mesh);
		// copies mesh

		int read(const std::string& filename, bool postProcessing=true);
		// reads a mesh from a Wavefront OBJ file; return value is nonzero
		// only if there was an error

		int write(const std::string& filename);
		// writes a mesh to a Wavefront OBJ file; return value is nonzero
		// only if there was an error

		bool reload(void);
		// reloads a mesh from disk using the most recent input filename

		void normalize(void);
		// centers around the origin and rescales to have unit radius

		int eulerCharacteristic(void) const;
		// returns V-E+F

		std::vector<HalfEdge> halfedges;
		std::vector<Vertex>   vertices;
		std::vector<Edge>     edges;
		std::vector<Face>     faces;
		// storage for mesh elements

		unsigned int fieldDegree;
		double fieldOffsetAngle;
		void computeSmoothestSection(void);
		void computeCurvatureAlignedSection(void);
		void computeTrivialSection(void);
		void extractSingularities(void);
		void alignTrivialSection(void);
		void computeParameterization(int coordinate);
		void glueParameterization(void);
		void parameterize(void);
		void buildEnergy(SparseMatrix<Real>& A, int coordinate);
		void buildDirichletEnergy(SparseMatrix<Real>& A);

		double energy(const SparseMatrix<Real>& A,
			const DenseMatrix<Real>& x,
			double eps);

		double lambda;
		// target line frequency

		int nCoordinateFunctions;
		// specifies whether to compute just the stripe
		// pattern (n=1) or both the stripe pattern and
		// an orthogonal coordinate (n=2), so that we
		// can construct a 2D parameterization



	protected:
		std::string inputFilename;

		void indexElements(void);
		void initializeSmoothStructure(void);
		void buildMassMatrix(void);
		void buildFieldEnergy(void);
		void assignTextureCoordinates(int coordinate);
		void buildDualLaplacian(void);

		SparseMatrix<Complex> massMatrix;
		SparseMatrix<Real> realMassMatrix;
		SparseMatrix<Complex> energyMatrix;
		SparseMatrix<Real> dualLaplacian;

		int nComputedCoordinateFunctions;
		// keeps track of how many coordinate functions were
		// computed the last time we ran Mesh::parameterize()


		// for paper Trivial Connections on Discrete Surfaces
	public:
		void writeEOBJ(std::string path);

		void topologyChange(void);                             // must be called by any method that modifies mesh connectivity
		void geometryChange(void);                             // must be called by any method that modifies mesh geometry
		//int nFaces(void) const;                              // number of faces, excluding virtual faces
		int nGenerators(void) const;                           // returns the number of independent noncontractible loops
		bool hasBoundary(void) const;                          // returns true iff the mesh has boundary
		void computeFrameAngles(double initialAngle);          // computes a global direction field using the trivial connection
		void appendDualGenerators(std::vector<Cycle>& cycles); // appends generators on the dual graph to "cycles"
		static bool isDualBoundaryLoop(const Cycle& cycle);    // returns true only if cycle is along boundary
		double boundaryLoopCurvature(const Cycle& cycle);      // computes the Riemannian holonomy around a boundary loop
		void computeTrivialConnection(void);                   // compute trivial connection with specified singularities
		void appendDirectionalConstraints(std::vector<Cycle>& cycles, std::vector<double>& holonomies);

		double parallelTransport(double phi, HalfEdgeCIter he); // transport a direction across an edge using Levi-Civita
		double defect(const Cycle& c); // computes angle defect resulting from parallel transport around a dual cycle c

		void integralCurve(FaceIter initialFace,
			const Vector& initialPoint,
			double initialAngle,
			std::vector<Vector>& curve,
			std::vector<Vector>& normals,
			int maxPts = 3000);
		// computes an integral curve of the current direction field starting
		// with the specified tangent direction.  Stops after maxPts points have
		// been computed or when something "bad" happens numerically.  The vector
		// "normals" stores the normal to the surface at each point on the curve.
		std::vector<double> generatorIndices; // target holonomy around generators divided by 2pi
		double fieldAngle;                    // initial angle for global direction field

	private:
		Connection* connection;                   // represents a trivial connection on the mesh

	protected:

		//void indexElements(void);                // assigns a unique ID to each mesh element
		void buildTreeCotreeDecomposition(void); // updates the tree-cotree decomposition
		void buildPrimalSpanningTree(void);      // builds a spanning tree of primal edges
		void buildDualSpanningCoTree(void);      // builds a spanning tree of dual edges that do not cross the primal tree

		// transportOrder caches data needed to update direction field from connection
		struct TransportData
		{
			double delta;
			double sign;
			double* omega;
			double* alphaI;
			double* alphaJ;
		};
		std::vector<TransportData> transportOrder;
		FaceIter transportRoot;
		
	};
}

#endif

