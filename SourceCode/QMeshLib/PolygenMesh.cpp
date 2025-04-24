// PMBody.cpp: implementation of the PMBody class.
//
//////////////////////////////////////////////////////////////////////
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include <memory.h>

#include "PolygenMesh.h"

#include <QDebug>

#define PI		3.141592654
#define DEGREE_TO_ROTATE(x)		0.0174532922222*x
#define ROTATE_TO_DEGREE(x)		57.295780490443*x

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PolygenMesh::PolygenMesh(mesh_type type)
{
    ClearAll();
    m_drawListID=-1;
    m_bVertexNormalShading=false;
    isTransparent = false;
    m_drawListNumber = 6;
    meshType = type;
}

PolygenMesh::~PolygenMesh()
{
    ClearAll();
    if (m_drawListID!=-1) glDeleteLists(m_drawListID, m_drawListNumber);
}

//////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////

void PolygenMesh::CompBoundingBox(double boundingBox[])
{
    GLKPOSITION PosMesh;
    GLKPOSITION Pos;
    double xx,yy,zz;

    boundingBox[0]=boundingBox[2]=boundingBox[4]=1.0e+32;
    boundingBox[1]=boundingBox[3]=boundingBox[5]=-1.0e+32;

    for(PosMesh=meshList.GetHeadPosition();PosMesh!=NULL;) {
        QMeshPatch *mesh=(QMeshPatch *)(meshList.GetNext(PosMesh));
        for(Pos=mesh->GetNodeList().GetHeadPosition();Pos!=NULL;) {
            QMeshNode *node=(QMeshNode *)(mesh->GetNodeList().GetNext(Pos));
            node->GetCoord3D(xx,yy,zz);

            if (xx<boundingBox[0]) boundingBox[0]=xx;
            if (xx>boundingBox[1]) boundingBox[1]=xx;
            if (yy<boundingBox[2]) boundingBox[2]=yy;
            if (yy>boundingBox[3]) boundingBox[3]=yy;
            if (zz<boundingBox[4]) boundingBox[4]=zz;
            if (zz>boundingBox[5]) boundingBox[5]=zz;
        }
    }
}

void PolygenMesh::DeleteGLList()
{
    if (m_drawListID!=-1) {
        glDeleteLists(m_drawListID, m_drawListNumber);
        m_drawListID=-1;
    }
}

void PolygenMesh::BuildGLList(bool bVertexNormalShading)
{
    if (m_drawListID!=-1) glDeleteLists(m_drawListID, m_drawListNumber);
    m_drawListID = glGenLists(m_drawListNumber);

    _buildDrawShadeList(bVertexNormalShading);
    _buildDrawMeshList();
    _buildDrawNodeList();
    _buildDrawProfileList();
    _buildDrawFaceNormalList();
    _buildDrawNodeNormalList();
    computeRange();
}

void PolygenMesh::_buildDrawShadeList(bool bVertexNormalShading)
{
    glNewList(m_drawListID, GL_COMPILE);

    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);

    drawOriginalCoordinate();
    if (isTransparent) {
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_FALSE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_CULL_FACE);
    }

    if (this->meshType == CURVED_LAYER) {
        for (GLKPOSITION Pos = meshList.GetHeadPosition(); Pos != NULL; ) {
            QMeshPatch* mesh = (QMeshPatch*)(meshList.GetNext(Pos));

            if (mesh->drawThisPatch == false) continue;
            float rr, gg, bb;

            glBegin(GL_TRIANGLES);
            for (GLKPOSITION PosFace = (mesh->GetFaceList()).GetHeadPosition(); PosFace != NULL;) {
                QMeshFace* face = (QMeshFace*)((mesh->GetFaceList()).GetNext(PosFace));

                _changeValueToColor(mesh->GetIndexNo() + 1, rr, gg, bb);
                //if (face->stessField_flag) { rr = 0.0; gg = 0.0; bb = 0.0; }
                if (face->singularWarning) { rr = 0.0; gg = 0.0; bb = 1.0; }

                glColor3f(rr, gg, bb);

                this->drawSingleFace(face);
            }
            glEnd();
        }
    }

    if (isTransparent) {
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
    }

    glEndList();
}

void PolygenMesh::_changeValueToColor(int nType, float & nRed, float & nGreen, float & nBlue)
{
    float color[][3]={
        {220,20,60},
        {107,200,35},
        {30,144,255},
        {255,105,180},
        {244,164,96},
        {176,196,222},
        {255,100,70},
        {128,255,128},
        {128,128,255},
        {255,255,128},
        {0,128,0},
        {255,128,255},
        {255,214,202},
        {128,128,192},
        {255,165,0}, //orange
        {255,128,192},
//		{39, 64, 139},//RoyalBlue
        {128,128,64},
        {0,255,255},
        {238,130,238},//violet
        {220,220,220},//gainsboro
        {188, 143, 143}, // rosy brown
        {46, 139, 87},//sea green
        {210, 105, 30 },//chocolate
        {237, 150, 100},
        {100, 149, 237},//cornflower blue
        {243, 20, 100},
        // 26th
        {0,0,0}
    };

//	printf("%d ",nType);
    nRed=color[nType%25][0]/255.0f;
    nGreen=color[nType%25][1]/255.0f;
    nBlue=color[nType%25][2]/255.0f;
}

void PolygenMesh::_buildDrawMeshList()
{
    if (meshList.GetCount() == 0) return;

    glNewList(m_drawListID + 1, GL_COMPILE);
    glDisable(GL_LIGHTING);
    if (edgeColor)  glLineWidth(1.0);

    if (this->meshType == CURVED_LAYER) {
        glLineWidth(0.5);
        for (GLKPOSITION Pos = meshList.GetHeadPosition(); Pos != NULL; ) {
            QMeshPatch* mesh = (QMeshPatch*)(meshList.GetNext(Pos));

            if (mesh->drawThisPatch == false) continue;
            float rr, gg, bb;

            glBegin(GL_LINES);

            for (GLKPOSITION PosEdge = (mesh->GetEdgeList()).GetHeadPosition(); PosEdge != NULL;) {
                QMeshEdge* edge = (QMeshEdge*)((mesh->GetEdgeList()).GetNext(PosEdge));

                glColor3f(0.75, 0.75, 0.75);

                //this->drawSingleEdge(edge);
            }

            glEnd();
        }
    }

    if (this->meshType == TOOL_PATH) {
        glLineWidth(3.0);
        for (GLKPOSITION Pos = meshList.GetHeadPosition(); Pos != NULL; ) {
            QMeshPatch* mesh = (QMeshPatch*)(meshList.GetNext(Pos));

            if (mesh->drawThisPatch == false) continue;
            float rr, gg, bb;

            glBegin(GL_LINES);

            for (GLKPOSITION PosEdge = (mesh->GetEdgeList()).GetHeadPosition(); PosEdge != NULL;) {
                QMeshEdge* edge = (QMeshEdge*)((mesh->GetEdgeList()).GetNext(PosEdge));

                if (edge->is_shortChain) continue;
                //_changeValueToColor(mesh->GetIndexNo() + 3, rr, gg, bb);

                _changeValueToColor(edge->seg_Idx, rr, gg, bb);

                glColor3f(rr, gg, bb);
                if (edge->drawThisEdge) this->drawSingleEdge(edge);
            }
            glEnd();
        }
    }

    glEndList();
}

void PolygenMesh::_buildDrawNodeList()
{
    if (meshList.GetCount() == 0) return;

    glNewList(m_drawListID + 2, GL_COMPILE);
    glDisable(GL_LIGHTING);

    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(5.0);

    if (this->meshType == CURVED_LAYER) {
        glBegin(GL_POINTS);
        for (GLKPOSITION Pos = meshList.GetHeadPosition(); Pos != NULL; ) {
            QMeshPatch* mesh = (QMeshPatch*)(meshList.GetNext(Pos));

            if (mesh->drawThisPatch == false) continue;
            float rr, gg, bb;

            for (GLKPOSITION PosNode = (mesh->GetNodeList()).GetHeadPosition(); PosNode != NULL;) {
                QMeshNode* node = (QMeshNode*)((mesh->GetNodeList()).GetNext(PosNode));

                float rr = gg = bb = 0.8;

                glColor3f(rr, gg, bb);
                //drawSingleNode(node);
            }
        }
        glEnd();
    }

    if (this->meshType == TOOL_PATH) {
        glPointSize(4.0);
        glBegin(GL_POINTS);
        for (GLKPOSITION Pos = meshList.GetHeadPosition(); Pos != NULL; ) {
            QMeshPatch* mesh = (QMeshPatch*)(meshList.GetNext(Pos));

            if (mesh->drawThisPatch == false) continue;
            float rr, gg, bb;

            for (GLKPOSITION PosNode = (mesh->GetNodeList()).GetHeadPosition(); PosNode != NULL;) {
                QMeshNode* node = (QMeshNode*)((mesh->GetNodeList()).GetNext(PosNode));

                if (node->is_shortChain) continue;
                
                float rr = gg = bb = 0.8;

                if (node->seg_Idx >= 0)_changeValueToColor(node->seg_Idx, rr, gg, bb);
                if (node->is_endPnt || node->Jump_nextSecStart || node->Jump_preSecEnd || node->cut_info){
                    rr = 0.0f; gg = 0.0f; bb = 0.0f;
                }
                if (node->is_endPnt) { rr = gg = bb = 0.0; }

                glColor3f(rr, gg, bb);
                drawSingleNode(node);
            }
        }
        glEnd();
    }

    glEndList();
}

void PolygenMesh::_buildDrawProfileList()
{
    glNewList(m_drawListID + 3, GL_COMPILE);

    glEnable(GL_TEXTURE_2D);
    glGenTextures(1, texid);
    glBindTexture(GL_TEXTURE_2D, texid[0]);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, m_texture_size[0], m_texture_size[1], 0, GL_BGRA, GL_UNSIGNED_BYTE, m_texture_data);

    if (this->meshType == CURVED_LAYER) {

        for (GLKPOSITION PosMesh = meshList.GetHeadPosition(); PosMesh != NULL;) {
            QMeshPatch* mesh = (QMeshPatch*)(meshList.GetNext(PosMesh));

            if (mesh->drawThisPatch == false) continue;

            double edgeLength = 0;
            for (GLKPOSITION Pos = mesh->GetEdgeList().GetHeadPosition(); Pos != NULL;) {
                QMeshEdge* Edge = (QMeshEdge*)(mesh->GetEdgeList().GetNext(Pos));
                edgeLength += Edge->CalLength();
            }
            edgeLength /= mesh->GetEdgeNumber(); // average edge length

            glLineWidth(2.5);

            if (mesh->draw_stripPattern) {
                
                glBegin(GL_TRIANGLES);
                for (GLKPOSITION PosFace = (mesh->GetFaceList()).GetHeadPosition(); PosFace != NULL;) {
                    QMeshFace* face = (QMeshFace*)((mesh->GetFaceList()).GetNext(PosFace));
                    this->drawSingleFaceWithTexture(face);
                }
                glEnd();
            
            }
            else {

                if (mesh->draw_Field_onFace) {
                    for (GLKPOSITION Pos = mesh->GetFaceList().GetHeadPosition(); Pos;) {
                        QMeshFace* face = (QMeshFace*)mesh->GetFaceList().GetNext(Pos);

                        if (mesh->draw_directionField_onFace && face->GetIndexNo() % 5 != 0) continue;

                        float rr, gg, bb;   rr = gg = bb = 0.35;
                        double length = edgeLength * 2.5;
                        Eigen::Vector3d n = Eigen::Vector3d::Zero();

                        if (mesh->draw_directionField_onFace) {
                            n = face->vectorDir;
                            n = n / 2.0;
                        }
                        else {
                            n = face->stessField_vector;
                            _changeValueToColor(1.0, 0.0, face->stessField_vector.norm(), rr, gg, bb);
                        }

                        double x, y, z;
                        face->CalCenterPos(x, y, z);

                        /*--------Draw vector field--------*/
                        if (n.norm() < 0.01) continue;

                        //----1. Draw Array (BODY) ----//
                        glColor3f(rr, gg, bb);

                        glBegin(GL_LINES);
                        //glVertex3d(x, y, z); 
                        glVertex3d(x - n[0] * length, y - n[1] * length, z - n[2] * length);
                        glVertex3d(x + n[0] * length, y + n[1] * length, z + n[2] * length);
                        glEnd();

                        //---- 2. Draw Array (TIP) ----//   
                        double endPoint[3] = { x + n[0] * length, y + n[1] * length, z + n[2] * length };
                        double n_temp[3]; for (int i = 0; i < 3; i++) n_temp[i] = n(i);
                        //drawSingleArrayTip(endPoint, n_temp, length);
                    }
                }
                else {
                    for (GLKPOSITION Pos = mesh->GetNodeList().GetHeadPosition(); Pos;) {
                        auto* node = (QMeshNode*)mesh->GetNodeList().GetNext(Pos);

                        if (node->GetIndexNo() % 20 != 0) continue;
                        float rr, gg, bb;   rr = 0.9;  gg = bb = 0.1;
                        double length = edgeLength * 10.0;
                        Eigen::Vector3d n = Eigen::Vector3d::Zero();

                        n = node->directionField;

                        double x, y, z;
                        node->GetCoord3D(x, y, z);

                        /*--------Draw vector field--------*/
                        if (n.norm() < 0.01) continue;

                        //----1. Draw Array (BODY) ----//
                        glColor3f(rr, gg, bb);

                        glBegin(GL_LINES);
                        glVertex3d(x, y, z); glVertex3d(x + n[0] * length, y + n[1] * length, z + n[2] * length);
                        glEnd();

                        //---- 2. Draw Array (TIP) ----//   
                        double endPoint[3] = { x + n[0] * length, y + n[1] * length, z + n[2] * length };
                        double n_temp[3]; for (int i = 0; i < 3; i++) n_temp[i] = n(i);
                        //drawSingleArrayTip(endPoint, n_temp, length);
                    }
                }
            }
        }
    }

    if (this->meshType == TOOL_PATH) {
        for (GLKPOSITION PosMesh = meshList.GetHeadPosition(); PosMesh != NULL;) {
            QMeshPatch* mesh = (QMeshPatch*)(meshList.GetNext(PosMesh));

            if (mesh->drawThisPatch == false) continue;

            double edgeLength = 0;
            for (GLKPOSITION Pos = mesh->GetEdgeList().GetHeadPosition(); Pos != NULL;) {
                QMeshEdge* Edge = (QMeshEdge*)(mesh->GetEdgeList().GetNext(Pos));
                edgeLength += Edge->CalLength();
            }
            edgeLength /= mesh->GetEdgeNumber(); // average edge length

            glLineWidth(0.5);
            for (GLKPOSITION Pos = mesh->GetNodeList().GetHeadPosition(); Pos;) {
                QMeshNode* node = (QMeshNode*)mesh->GetNodeList().GetNext(Pos);
                float rr, gg, bb;
                rr = 0.0;  gg = bb = 1.0;

                double x, y, z;
                node->GetCoord3D(x, y, z);

                double n[3], length;
                length = edgeLength;

                /*--------Draw vector field--------*/
                if (node->m_printTan.norm() < 0.01) continue;
                for (int i = 0; i < 3; i++) n[i] = node->m_printTan(i);

                //---- Draw Array (BODY) ----//
                glColor3f(rr, gg, bb);

                glBegin(GL_LINES);
                glVertex3d(x, y, z); glVertex3d(x + n[0] * length, y + n[1] * length, z + n[2] * length);
                glEnd();

                //---- Draw Array (TIP) ----//   
                double endPoint[3] = { x + n[0] * length, y + n[1] * length, z + n[2] * length };
                drawSingleArrayTip(endPoint, n, length);
            }
        }
    }

    glDisable(GL_TEXTURE_2D);
    glEndList();
}

void PolygenMesh::_buildDrawFaceNormalList()
{
    if (meshList.GetCount() == 0) return;

    glNewList(m_drawListID + 4, GL_COMPILE);
    glDisable(GL_LIGHTING);

    glColor3f(0.5, 0.0, 0.5);

    glLineWidth(1.0);

    if (this->meshType == CURVED_LAYER) {

        glBegin(GL_LINES);
        for (GLKPOSITION meshPos = meshList.GetHeadPosition(); meshPos != NULL;) {
            QMeshPatch* mesh = (QMeshPatch*)(meshList.GetNext(meshPos));

            if (mesh->drawThisPatch == false) continue;

            QMeshEdge* edge = (QMeshEdge*)mesh->GetEdgeList().GetHead();
            double length = edge->CalLength();
            for (GLKPOSITION Pos = mesh->GetFaceList().GetHeadPosition(); Pos != NULL;) {
                QMeshFace* face = (QMeshFace*)(mesh->GetFaceList().GetNext(Pos));
                double x, y, z, nx, ny, nz;
                face->CalCenterPos(x, y, z);
                face->CalPlaneEquation();
                face->GetNormal(nx, ny, nz);
                glVertex3d(x, y, z);
                glVertex3d(x + nx * length, y + ny * length, z + nz * length);
            }
        }
        glEnd();
    }

    glEndList();
}

void PolygenMesh::_buildDrawNodeNormalList()
{
    if (meshList.GetCount() == 0) return;

    glNewList(m_drawListID + 5, GL_COMPILE);
    glDisable(GL_LIGHTING);

    glColor3f(0.0, 0.5, 0.0);

    glLineWidth(1.0);

    if (this->meshType == TOOL_PATH) {
        glBegin(GL_LINES);
        for (GLKPOSITION meshPos = meshList.GetHeadPosition(); meshPos != NULL;) {
            QMeshPatch* mesh = (QMeshPatch*)(meshList.GetNext(meshPos));


            if (mesh->drawThisPatch == false) continue;

            QMeshEdge* edge = (QMeshEdge*)mesh->GetEdgeList().GetHead();
            //double length = edge->CalLength();
            double length = 1.0;
            for (GLKPOSITION Pos = mesh->GetNodeList().GetHeadPosition(); Pos != NULL;) {
                QMeshNode* node = (QMeshNode*)(mesh->GetNodeList().GetNext(Pos));

                //if (node->resampleChecked == false) continue;

                double x, y, z;
                node->GetCoord3D(x, y, z);
                double n[3];
                node->GetNormal(n[0], n[1], n[2]);

                glVertex3d(x, y, z);
                glVertex3d(x + n[0] * length, y + n[1] * length, z + n[2] * length);
            }
        }
        glEnd();
    }
    else {

        glBegin(GL_LINES);
        for (GLKPOSITION meshPos = meshList.GetHeadPosition(); meshPos != NULL;) {
            QMeshPatch* mesh = (QMeshPatch*)(meshList.GetNext(meshPos));

            if (this->meshType == CURVED_LAYER && mesh->drawThisPatch == false) continue;

            QMeshEdge* edge = (QMeshEdge*)mesh->GetEdgeList().GetHead();
            double length = edge->CalLength();
            //double length = 1.0;
            for (GLKPOSITION Pos = mesh->GetNodeList().GetHeadPosition(); Pos != NULL;) {
                QMeshNode* node = (QMeshNode*)(mesh->GetNodeList().GetNext(Pos));
                double x, y, z;
                node->GetCoord3D(x, y, z);
                double n[3];
                node->CalNormal(n);

                glVertex3d(x, y, z);
                glVertex3d(x + n[0] * length, y + n[1] * length, z + n[2] * length);
            }
        }
        glEnd();
    }

    glEndList();
}

void PolygenMesh::drawOriginalCoordinate() {
	double axisLeng = 30.0;
	glLineWidth(2.0);
	glBegin(GL_LINES);

	// X-axis - Red Color
	glColor3f(1.0, 0.0, 0.0); glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(axisLeng, 0.0, 0.0);

	// Y-axis - green Color
	glColor3f(0.0, 1.0, 0.0);
	glVertex3d(0.0, 0.0, 0.0); glVertex3d(0.0, axisLeng, 0.0);

	// Z-axis - black Color
	glColor3f(0.0, 0.0, 0.0);
	glVertex3d(0.0, 0.0, 0.0); glVertex3d(0.0, 0.0, axisLeng);

	glEnd();
}

void PolygenMesh::drawShade()
{
    if (meshList.IsEmpty()) {glDeleteLists(m_drawListID, m_drawListNumber); m_drawListID=-1; return;}
    glCallList(m_drawListID);
}

void PolygenMesh::drawMesh()
{
    if (meshList.IsEmpty()) {glDeleteLists(m_drawListID, m_drawListNumber); m_drawListID=-1; return;}
    glCallList(m_drawListID+1);
}

void PolygenMesh::drawNode()
{
    if (meshList.IsEmpty()) {glDeleteLists(m_drawListID, m_drawListNumber); m_drawListID=-1; return;}
    glCallList(m_drawListID+2);
}

void PolygenMesh::drawProfile()
{
    if (meshList.IsEmpty()) {glDeleteLists(m_drawListID, m_drawListNumber); m_drawListID=-1; return;}
    glCallList(m_drawListID+3);
}

void PolygenMesh::drawFaceNormal()
{
    if (meshList.IsEmpty()) {glDeleteLists(m_drawListID, m_drawListNumber); m_drawListID=-1; return;}
    glCallList(m_drawListID+4);
}

void PolygenMesh::drawNodeNormal()
{
    if (meshList.IsEmpty()) {glDeleteLists(m_drawListID, m_drawListNumber); m_drawListID=-1; return;}
    glCallList(m_drawListID+5);
}

void PolygenMesh::ClearAll()
{
    GLKPOSITION Pos;

    for(Pos=meshList.GetHeadPosition();Pos!=NULL;) {
        QMeshPatch *mesh=(QMeshPatch *)(meshList.GetNext(Pos));
        delete mesh;
    }
    meshList.RemoveAll();
}

void PolygenMesh::computeRange()
{
    double range=0.0,ll,xx,yy,zz;
    GLKPOSITION Pos;
    GLKPOSITION PosNode;

    for(Pos=meshList.GetHeadPosition();Pos!=NULL;) {
        QMeshPatch *mesh=(QMeshPatch *)(meshList.GetNext(Pos));
        for(PosNode=(mesh->GetNodeList()).GetHeadPosition();PosNode!=NULL;) {
            QMeshNode *node=(QMeshNode *)((mesh->GetNodeList()).GetNext(PosNode));

            node->GetCoord3D(xx,yy,zz);
            ll=xx*xx+yy*yy+zz*zz;

            if (ll>range) range=ll;
        }
    }

    m_range=(float)(sqrt(range));
}

void PolygenMesh::_changeValueToColor(double maxValue, double minValue, double Value,
                                 float & nRed, float & nGreen, float & nBlue)
{
//	Value=fabs(Value);

    if (Value<=minValue)
    {
        nRed=0.0;
        nGreen=0.0;
        nBlue=0.0;
        return;
    }

    if ((maxValue-minValue)<0.000000000001)
    {
        nRed=0.0;
        nGreen=0.0;
        nBlue=1.0;
        return;
    }

    double temp=(Value-minValue)/(maxValue-minValue);

//    nRed=(float)(1.0-temp);	nGreen=(float)(1.0-temp); nBlue=(float)(1.0-temp);	return;

    if (temp>0.75)
    {
        nRed=1;
        nGreen=(float)(1.0-(temp-0.75)/0.25);
        if (nGreen<0) nGreen=0.0f;
        nBlue=0;
        return;
    }
    if (temp>0.5)
    {
        nRed=(float)((temp-0.5)/0.25);
        nGreen=1;
        nBlue=0;
        return;
    }
    if (temp>0.25)
    {
        nRed=0;
        nGreen=1;
        nBlue=(float)(1.0-(temp-0.25)/0.25);
        return;
    }
    else
    {
        nRed=0;
        nGreen=(float)(temp/0.25);
        nBlue=1;
    }
}

void PolygenMesh::ImportOBJFile(char *filename, std::string modelName)
{
    QMeshPatch *newMesh = new QMeshPatch;
    if (newMesh->inputOBJFile(filename)){
        meshList.AddTail(newMesh);
        computeRange();
        setModelName(modelName);
    }
    else
        delete newMesh;
}

void PolygenMesh::ImportTETFile(char *filename, std::string modelName)
{
	QMeshPatch *newMesh = new QMeshPatch;
	if (newMesh->inputTETFile(filename, false)) {
		meshList.AddTail(newMesh);
		computeRange();
		setModelName(modelName);
	}
	else
		delete newMesh;
}

void PolygenMesh::drawSingleEdge(QMeshEdge* edge) {

    double xx, yy, zz;
    edge->GetStartPoint()->GetCoord3D(xx, yy, zz);
    glVertex3d(xx, yy, zz);
    edge->GetEndPoint()->GetCoord3D(xx, yy, zz);
    glVertex3d(xx, yy, zz);

}

void PolygenMesh::drawSingleNode(QMeshNode* node) {

    double nx, ny, nz, xx, yy, zz;
    node->GetNormal(nx, ny, nz);
    node->GetCoord3D(xx, yy, zz);
    glNormal3d(nx, ny, nz);
    glVertex3d(xx, yy, zz);

}

void PolygenMesh::drawSingleFace(QMeshFace* face) {

    double xx, yy, zz, dd;
    for (int i = 0; i < 3; i++) {
        QMeshNode* node = face->GetNodeRecordPtr(i);

        if (m_bVertexNormalShading) {
            double normal[3];
            node->CalNormal(normal); 
            glNormal3dv(normal);
        }
        else {
            face->CalPlaneEquation();
            face->GetPlaneEquation(xx, yy, zz, dd);
            glNormal3d(xx, yy, zz);
        }

        node->GetCoord3D(xx, yy, zz);
        glVertex3d(xx, yy, zz);
    }
}

void PolygenMesh::drawSingleArrayTip(double pp[3], double dir[3], double arrowLength) {

    double bone_hight = 0.3 * arrowLength;
    double bone_radius = 0.06 * arrowLength;

    Eigen::Vector3d endPP = { pp[0],pp[1],pp[2] };

    Eigen::Vector3d A = { dir[0],dir[1],dir[2] };
    Eigen::Vector3d B = { 0, 1.0, 0 };

    Eigen::Matrix3d rotationMatrix;
    rotationMatrix = Eigen::Quaterniond().setFromTwoVectors(B, A);

    Eigen::Vector3d pp1 = { bone_radius * sin(0) , 0.0 , bone_radius * cos(0) };
    Eigen::Vector3d pp2 = { bone_radius * sin(120 * 3.14 / 180) , 0.0 , bone_radius * cos(120 * 3.14 / 180) };
    Eigen::Vector3d pp3 = { bone_radius * sin(240 * 3.14 / 180) , 0.0 , bone_radius * cos(240 * 3.14 / 180) };
    Eigen::Vector3d ppCenter = { 0.0, bone_hight, 0.0 };

    pp1 = rotationMatrix * pp1 + endPP;
    pp2 = rotationMatrix * pp2 + endPP;
    pp3 = rotationMatrix * pp3 + endPP;
    ppCenter = rotationMatrix * ppCenter + endPP;

    glBegin(GL_TRIANGLES);

    glColor3f(0.9, 0.2, 0.4);

    glVertex3d(pp1(0), pp1(1), pp1(2));
    glVertex3d(pp2(0), pp2(1), pp2(2));
    glVertex3d(ppCenter(0), ppCenter(1), ppCenter(2));

    glVertex3d(pp2(0), pp2(1), pp2(2));
    glVertex3d(pp3(0), pp3(1), pp3(2));
    glVertex3d(ppCenter(0), ppCenter(1), ppCenter(2));

    glVertex3d(pp3(0), pp3(1), pp3(2));
    glVertex3d(pp1(0), pp1(1), pp1(2));
    glVertex3d(ppCenter(0), ppCenter(1), ppCenter(2));

    glEnd();
}

void PolygenMesh::SetTexture(int nx, int ny, uchar* data, int data_size)
{
    m_texture_size[0] = nx;
    m_texture_size[1] = ny;
    if (m_texture_data != NULL) {
        delete m_texture_data;
    }
    m_texture_data = new uchar[data_size];
    memcpy(m_texture_data, data, data_size * sizeof(uchar));
}

void PolygenMesh::drawSingleFaceWithTexture(QMeshFace* face) {
    double xx, yy, zz, dd;
    for (int i = 0; i < 3; i++) {
        QMeshNode* node = face->GetNodeRecordPtr(i);

        if (m_bVertexNormalShading) {
            double normal[3];
            node->CalNormal(normal); glNormal3dv(normal);
        }
        else {
            face->CalPlaneEquation();
            face->GetPlaneEquation(xx, yy, zz, dd);
            glNormal3d(xx, yy, zz);
        }

        node->GetCoord3D(xx, yy, zz);
        glTexCoord2d(face->GetTextureCoordu(i), face->GetTextureCoordv(i));
        glVertex3d(xx, yy, zz);
    }
}