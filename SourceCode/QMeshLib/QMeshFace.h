
// QMeshFace.h: interface for the QMeshFace class
//
//////////////////////////////////////////////////////////////////////

#ifndef _QMESHFACE
#define _QMESHFACE
#include <stddef.h>

#include "../GLKLib/GLKObList.h"

#define MAX_EDGE_NUM	10

class QMeshPatch;
class QMeshEdge;
class QMeshNode;
class QMeshTetra;

class QMeshFace : public GLKObject
{
public:
	QMeshFace();
	virtual ~QMeshFace();

public:
	bool GetAttribFlag( const int whichBit );
	void SetAttribFlag( const int whichBit, const bool toBe = true );

	int GetIndexNo();		//from 1 to n
	void SetIndexNo( const int _index = 1 );

	bool IsNormalDirection( const int whichEdge );
	void SetDirectionFlag( const int whichEdge, const bool toBe = true );

	QMeshEdge * GetEdgeRecordPtr( const int whichEdge );
    void SetEdgeRecordPtr( const int whichEdge, QMeshEdge * _edge = nullptr );
	int GetEdgeNum();
	void SetEdgeNum(int num);
		
	void GetNodePos( const int whichNode, double &xx, double &yy, double &zz );
	QMeshNode * GetNodeRecordPtr( const int whichNode );

	double GetNodeAngle( const int whichNode );
	void SetNodeAngle( const int whichNode, double angleInRadian );

	void SetColor(float r, float g, float b);
	void GetColor(float &r, float &g, float &b);

	void SetHermiteData(double pos[], double normal[]);
	void GetHermiteData(double pos[], double normal[]);

	void GetPlaneEquation( double &A, double &B, double &C, double &D );
	void CalPlaneEquation( ); 
						// to calculate the plane equation parameter
	                    // Plane equation:  Ax + By + Cz + D = 0, and
                        // Vector(A,B,C) is positive unit normal vector of this plane

	void CalCenterPos(double &xx, double &yy, double &zz);
	void CalCenterPos();
    void GetCenterPos(double &xx, double &yy, double &zz) {xx=center[0];yy=center[1];zz=center[2];};
    void SetCenterPos(double xx, double yy, double zz) {center[0]=xx; center[1]=yy; center[2]=zz;};
	void CalCenterPos_last(double &xx, double &yy, double &zz);
	void CalBoundingBox(double &xmin, double &ymin, double &zmin,
						double &xmax, double &ymax, double &zmax);
	double CalArea();
	double GetArea() {return m_area;};
	double CalArea2D();

	void SetMeshPatchPtr(QMeshPatch* _mesh);
	QMeshPatch* GetMeshPatchPtr();

    GLKObList& GetAttachedList() {return attachedList;};
	int nSplittedNode;

	void SetNormal(double nx, double ny, double nz) {abcd[0]=nx;abcd[1]=ny;abcd[2]=nz;};
    void GetNormal(double &nx, double &ny, double &nz){nx=abcd[0]; ny=abcd[1]; nz=abcd[2];};
	void SetPlane(double nx, double ny, double nz, double dd) {abcd[0]=nx;abcd[1]=ny;abcd[2]=nz;abcd[3]=dd;};
	double m_desiredNormal[3];

    double m_GaussArea;

	int m_nIdentifiedPatchIndex;
    bool selected;
	bool inner, i_inner;
	bool isFixedDraw = false;
	bool isHandleDraw = false;
	double weight;

	void *attachedPointer;
    bool visited;
    double m_weight;
    bool m_centroid;
    double value1;
    double value2;
    double boundaryDis;
    bool m_bndface;
    bool IsMoreWeight;
    bool visible;

    int identifiedIndex;
	GLKObList attachedList;	// a list of attached object

private:
	int indexno;
	bool flags[8];
		// 0 - for boundary faces
		// 1 - for sharp-feature region faces
		// 2 - for this face need to be subdivided from center
		// 3 - for temp use
	
		// 0, 1, 2, 3, ... edges construct an anti-clockwise closed loop
                                      //**********************************
                                      //             point1              *
                                      //              /\                 *
                                      //    1 edge  /    \ _             *
                                      //         |_       |\ 0 edge      *
                                      //        /            \           *
                                      //      /                \         *
                                      //  point2 ------>-------- point0  *
	                                  //              2 edge             *
	                                  //**********************************
	QMeshEdge * edges[MAX_EDGE_NUM];	//	edges
	bool edgeDir[MAX_EDGE_NUM];			//	edge directions
	int edgeNum;
	double nodeAngle[MAX_EDGE_NUM];

	QMeshPatch *meshSurface;	// MESHSURFACE contain this triangle
	double  abcd[4];	// plane equation
	float	rgb[3];		// the color of the face
	double m_area;
    double center[3];

	double m_HermitePos[3],m_HermiteNormal[3];
	double m_texture_coordinate[3][2];
	double m_vertice_normal[3][3];

	//for volume mesh
public:
	QMeshTetra * GetLeftTetra();
	void SetLeftTetra(QMeshTetra * _pLeftTetra = NULL);
	QMeshTetra * GetRightTetra();
	void SetRightTetra(QMeshTetra * _pRightTetra = NULL);

	//ddg class
	void SetVerticeNormal(int n, double nx, double ny, double nz) {
		m_vertice_normal[n][0] = nx;
		m_vertice_normal[n][1] = ny;
		m_vertice_normal[n][2] = nz;
	}

	void GetVerticeNormal(int n, double& nx, double& ny, double& nz) {
		nx = m_vertice_normal[n][0];
		ny = m_vertice_normal[n][1];
		nz = m_vertice_normal[n][2];
	}

	double GetTextureCoordu(int n) {
		return m_texture_coordinate[n][0];
	}
	double GetTextureCoordv(int n) {
		return m_texture_coordinate[n][1];
	}
	void SetTextureCoord(int n, double u, double v) {
		m_texture_coordinate[n][0] = u;
		m_texture_coordinate[n][1] = v;
	}

private:
	//for volume mesh
	QMeshTetra *pLeftTetra;
	QMeshTetra *pRightTetra;

public:
	Eigen::Vector3d stessField_vector = Eigen::Vector3d::Zero();
	bool stessField_flag = false;
	bool boundary = false;

	bool isConstraint = false;//(all = b + s)
	bool isBoundaryConstraint = false;

	Eigen::Vector3d vectorDir = Eigen::Vector3d::Zero(); //(all = b + s)
	Eigen::Vector3d stressField = Eigen::Vector3d::Zero();

	double stressMag = 0;
	bool singularWarning = false;

	int FieldIndexNumber = 0;

	Eigen::Vector3d xDir, yDir, zDir;
	double slope1 = 0;
	double slope2 = 0;

	bool outBFace = false;
	bool inBFace = false;

	//ddg class
	double stripeIndices[2] = { 0,0 };
	double fieldIndices = 0;

	//strip toolpath
	bool isBoundary = false;
};

#endif
