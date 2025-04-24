#include "TrivialConnections.h"

#include "Mesh.h"
#include "Types.h"
#include "qdebug.h"

#include <iostream>
#include <sstream>

using namespace std;


    TrivialConnections::TrivialConnections(void)
    {}

    TrivialConnections::TrivialConnections(istream& in)
    {
        read(in);
    }

    void TrivialConnections::solve(std::string path)
    {
        DDG::Mesh mesh;

        mesh.read(path, false);
        // init parent
        mesh.topologyChange();
        mesh.geometryChange();
        qDebug() << "start solve";

        int V = mesh.vertices.size();
        for (vector<Singularity>::const_iterator v = singularities.begin();
            v != singularities.end();
            v++)
        {
            int i = v->first;
            double k = v->second;

            if (i < 0 || i >= V)
            {
                qDebug()<< "Warning: singularity requested at vertex " << QString::number(i + 1) << " (mesh has only " << QString::number(V) << " vertices!)" << endl;
                //cerr << "Warning: singularity requested at vertex " << i + 1 << " (mesh has only " << V << " vertices!)" << endl;
            }
            else
            {
                mesh.vertices[i].k = k;
            }
        }
        qDebug() << "nGenerators" << endl;

        int G = mesh.nGenerators();
        for (vector<Generator>::const_iterator g = generators.begin();
            g != generators.end();
            g++)
        {
            int i = g->first;
            double k = g->second;

            if (i < 0 || i >= G)
            {
                qDebug() << "Warning: additional holonomy requested around generator " << i + 1 << " (mesh has only " << G << " generators!)" << endl;
                //cerr << "Warning: additional holonomy requested around generator " << i + 1 << " (mesh has only " << G << " generators!)" << endl;
            }
            else
            {
                mesh.generatorIndices[i] = k;
            }
        }

        mesh.fieldAngle = fieldAngle;

        qDebug() << "computeTrivialConnection" << endl;
        mesh.computeTrivialConnection();

        mesh.writeEOBJ("out.eobj");
        
        toPolygenMesh(mesh);
    }

    void TrivialConnections::read(istream& in)
    {
        singularities.clear();
        generators.clear();

        string s;
        while (getline(in, s))
        {
            stringstream line(s);
            string token;

            line >> token;
            transform(token.begin(), token.end(), token.begin(), ::tolower);

            if (token == "in")
            {
                line >> inputPath;
            }
            else if (token == "out")
            {
                line >> outputPath;
            }
            else if (token == "vertex")
            {
                Singularity s;
                line >> s.first >> s.second;
                singularities.push_back(s);
            }
            else if (token == "generator")
            {
                Generator g;
                line >> g.first >> g.second;
                generators.push_back(g);
            }
            else if (token == "angle")
            {
                line >> fieldAngle;
            }
        }
    }


    PolygenMesh* TrivialConnections::toPolygenMesh(DDG::Mesh& mesh) {
        polygenmesh = new PolygenMesh(CURVED_LAYER);
        QMeshPatch* newMesh = new QMeshPatch;
        QMeshNode* node;

        int nv = mesh.vertices.size();
        int nf = mesh.faces.size();

        float* v = new float[mesh.vertices.size() * 3];
        uint* f = new uint[mesh.faces.size() * 3];

        std::vector<float> fn;
        std::vector<float> ft;

        for (int i = 0; i < nv; i++) {
            v[i * 3] = (float)(mesh.vertices[i].position.x);
            v[i * 3 + 1] = (float)(mesh.vertices[i].position.y);
            v[i * 3 + 2] = (float)(mesh.vertices[i].position.z);
        }

        for (int i = 0; i < mesh.faces.size(); i++) {
            auto face = mesh.faces[i];
            f[i * 3] = (uint)(face.he->vertex->index);
            f[i * 3 + 1] = (uint)(face.he->next->vertex->index);
            f[i * 3 + 2] = (uint)(face.he->next->next->vertex->index);

            double alpha = face.alpha;
            DDG::Vector w(cos(alpha), sin(alpha), 0.);
            DDG::Vector u = face.toGlobal(w);
            DDG::Vector du = u - face.he->vertex->position;

            du.normalize();
            fn.push_back(du.x);
            fn.push_back(du.y);
            fn.push_back(du.z);

            DDG::HalfEdgeIter he = face.he;

            for (int j = 0; j < 3; j++)
            {
                DDG::Complex z = he->texcoord / (2. * M_PI) + DDG::Complex(.5, .5);
                ft.push_back(z.re);
                ft.push_back(z.im);
                he = he->next;
            }
        }

        newMesh->constructionFromVerFaceTable(nv, v, nf, f);

        int faceIndex = 0;
        GLKPOSITION Pos;
        for (Pos = newMesh->GetFaceList().GetHeadPosition(); Pos != NULL; faceIndex++) {
            QMeshFace* face = (QMeshFace*)(newMesh->GetFaceList().GetNext(Pos));
            face->SetNormal(fn[3 * faceIndex], fn[3 * faceIndex + 1], fn[3 * faceIndex + 2]);

            for (int i = 0; i < 3; i++) {
                face->SetVerticeNormal(i, fn[3 * faceIndex], fn[3 * faceIndex + 1], fn[3 * faceIndex + 2]);
                face->SetTextureCoord(i, ft[6 * faceIndex + 2 * i], ft[6 * faceIndex + 2 * i + 1]);
            }
        }

        polygenmesh->meshList.AddTail(newMesh);
        polygenmesh->computeRange();
        polygenmesh->setModelName("modelName");

        return polygenmesh;
    }