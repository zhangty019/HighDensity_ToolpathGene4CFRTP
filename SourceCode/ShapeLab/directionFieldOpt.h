#pragma once
#include "../QMeshLib/PolygenMesh.h"
#include <set>
#include <Eigen/Eigen>
//#include <Eigen/PardisoSupport>
#include <fstream>  // For file handling
#include <sstream>  // For string stream to create unique file names

class directionFieldOpt
{
public:
	directionFieldOpt() {}
	~directionFieldOpt() {};

	void initialize(PolygenMesh* isoLayerSet, std::string model_name);
	void run();

private:

	void _installCoordinateSystem(QMeshPatch* materialMesh);
	void _installSlopes(QMeshPatch* materialMesh);
	void _classifyBoundary(QMeshPatch* materialMesh);
	void _smoothVectorField(QMeshPatch* materialMesh);
	void _buildLaplacianMatrix(QMeshPatch* materialMesh, 
		Eigen::SparseMatrix<double>& LaplacianMatrix);
	void _normalizeWeights(QMeshPatch* materialMesh);
	void _createSolutionVector(QMeshPatch* materialMesh,
		int column, Eigen::VectorXd& u);
	void _resolveSingularity(QMeshPatch* materialMesh);
	void _getNodeDirectionField(QMeshPatch* materialMesh);

	double calculateAngleChange(QMeshEdge* startEdge, QMeshEdge* nextEdge);
	void classifyEdges(std::set<QMeshEdge*>& edgeSet, double angleChange);

	void _setOptParameters(std::string model_name);

	void _calAngleDiff_dirctionField_vs_stressField(QMeshPatch* materialMesh);

private:

	PolygenMesh* m_isoLayerSet;

	//initial value
	//face->outBFace <-- face->isBoundaryConstraint
	//face->inBFace  <-- face->isBoundaryConstraint
	Eigen::Vector2d b_Weighting = { 1.0 * pow(10,-3), 1.0 * pow(10,-3) };
	//Weighting of stress = stressWeight_ratio * pow(10,[0,-stressWeight_range])
	double stressWeight_range = 2.0;
	double stressWeight_ratio = 0.1;
};