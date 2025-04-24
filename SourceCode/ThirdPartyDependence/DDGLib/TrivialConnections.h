#pragma once

// Like the problem.h in "Trivial Connections on Discrete Surfaces"

#include <vector>
#include <iosfwd>
#include <string>
#include "Mesh.h"
#include "../QMeshLib/PolygenMesh.h"

class TrivialConnections
{
public:
    TrivialConnections(void);               // default constructor
    TrivialConnections(std::istream& in);   // construct a problem from a valid istream
    void read(std::istream& in); // load a problem from a valid istream
    void solve(std::string path);      // solve the problem


    PolygenMesh* toPolygenMesh(DDG::Mesh& mesh);
    PolygenMesh* polygenmesh;


protected:
    std::string inputPath;
    std::string outputPath;

    double fieldAngle;

    typedef std::pair<int, double> Singularity;
    std::vector<Singularity> singularities;

    typedef std::pair<int, double> Generator;
    std::vector<Generator> generators;
};
