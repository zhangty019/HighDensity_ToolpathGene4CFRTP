#pragma once

#include "PolygenMesh.h"

class meshOperation
{

public:
	meshOperation() {};
	~meshOperation() {};

	void setBoundaryNodes(PolygenMesh* polymesh);
	void initialiseIndices(PolygenMesh* polymesh);
	void initial(PolygenMesh* sourceLayerSet);
	void offsetMesh(PolygenMesh* sourceLayerSet, PolygenMesh* offsetLayerSet);

private:

	double offset_value = 0.1;

	

};

