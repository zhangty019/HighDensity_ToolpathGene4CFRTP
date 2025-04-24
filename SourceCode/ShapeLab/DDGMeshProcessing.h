#pragma once
#include <Mesh.h>

#include "../QMeshLib/PolygenMesh.h"

class DDGMeshProcessing {
public:
	DDGMeshProcessing() {};
	~DDGMeshProcessing() {};

	void setMesh(std::string one_layerPath) {
		mesh.read(one_layerPath);
		mesh.nCoordinateFunctions = m_nCoordinateFunctions;
		mesh.lambda = m_frequence;
	}

	void setInitField(std::vector<std::vector<double>>& initField);
	void meshProcessing();
	QMeshPatch* toQMeshPatch();

	std::vector<std::vector<double>> perVerticeGradient;

private:
	int m_nCoordinateFunctions = 2;
	float m_frequence = 6; // frameBA/rest (2.7/5)
	DDG::Mesh mesh;
};