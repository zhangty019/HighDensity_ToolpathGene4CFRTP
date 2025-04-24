#include "meshOperation.h"

void meshOperation::setBoundaryNodes(PolygenMesh* polymesh) {

    //This needs to be improved on priority!!

    for (GLKPOSITION PatchPos = polymesh->GetMeshList().GetHeadPosition(); PatchPos;) {
        QMeshPatch* Patch = (QMeshPatch*)polymesh->GetMeshList().GetNext(PatchPos);

        for (GLKPOSITION pos = Patch->GetNodeList().GetHeadPosition(); pos;) {
            QMeshNode* node = (QMeshNode*)Patch->GetNodeList().GetNext(pos);

            int num = node->GetEdgeNumber();
            for (int i = 0; i < num; i++) {

                QMeshEdge* edge = (QMeshEdge*)node->GetEdgeRecordPtr(i + 1);

                if (edge->GetLeftFace() == NULL || edge->GetRightFace() == NULL) {
                    node->boundary = true;
                    edge->boundary = true;

                    if (edge->GetLeftFace()) {

                        edge->GetLeftFace()->boundary = true;
                        edge->GetLeftFace()->isBoundaryConstraint = true;

                        QMeshNode* sNode = edge->GetStartPoint();
                        QMeshNode* eNode = edge->GetEndPoint();
                        double x1, x2, y1, y2, z1, z2;
                        sNode->GetCoord3D(x1, y1, z1);
                        eNode->GetCoord3D(x2, y2, z2);

                        QMeshFace* thisFace = edge->GetLeftFace();

                        Eigen::Vector3d dirVec(x2 - x1, y2 - y1, z2 - z1);
                        dirVec.stableNormalize();
                        //thisFace->CalPlaneEquation();
                        //double nx, ny, nz;
                        //thisFace->GetNormal(nx, ny, nz);
                        //Eigen::Vector3d f_normal(nx, ny, nz);
                        //dirVec = dirVec.cross(f_normal);
                        dirVec.stableNormalize();
                        thisFace->vectorDir =  dirVec;
                        thisFace->isConstraint = true;

                    }
                    if (edge->GetRightFace()) {
                        edge->GetRightFace()->boundary = true;
                        edge->GetRightFace()->isBoundaryConstraint = true;

                        QMeshNode* sNode = edge->GetStartPoint();
                        QMeshNode* eNode = edge->GetEndPoint();
                        double x1, x2, y1, y2, z1, z2;
                        sNode->GetCoord3D(x1, y1, z1);
                        eNode->GetCoord3D(x2, y2, z2);

                        QMeshFace* thisFace = edge->GetRightFace();

                        Eigen::Vector3d dirVec(x2 - x1, y2 - y1, z2 - z1);
                        dirVec.stableNormalize();
                        thisFace->CalPlaneEquation();
                        double nx, ny, nz;
                        thisFace->GetNormal(nx, ny, nz);
                        Eigen::Vector3d f_normal(nx, ny, nz);
                        dirVec = dirVec.cross(f_normal);
                        dirVec.stableNormalize();
                        thisFace->vectorDir[0] = dirVec.x();
                        thisFace->vectorDir[1] = dirVec.y();
                        thisFace->vectorDir[2] = dirVec.z();
                        thisFace->isConstraint = true;
                    }

                }

            }
        }
    }
}

void meshOperation::initialiseIndices(PolygenMesh* polymesh){

    int patchId = 0;
    for (GLKPOSITION PatchPos = polymesh->GetMeshList().GetHeadPosition(); PatchPos;) {
        QMeshPatch* Patch = (QMeshPatch*)polymesh->GetMeshList().GetNext(PatchPos);

        Patch->SetIndexNo(patchId);

        int count = 0;
        for (GLKPOSITION nodePos = Patch->GetNodeList().GetHeadPosition(); nodePos;) {
            QMeshNode* thisNode = (QMeshNode*)Patch->GetNodeList().GetNext(nodePos);
            thisNode->FieldIndexNumber = count++;
        }
        //std::cout << "Node Count " << count << std::endl;
        count = 0;
        for (GLKPOSITION edgePos = Patch->GetEdgeList().GetHeadPosition(); edgePos;) {
            QMeshEdge* thisEdge = (QMeshEdge*)Patch->GetEdgeList().GetNext(edgePos);
            thisEdge->FieldIndexNumber = count++;
        }
        //std::cout << "Edge Count " << count << std::endl;
        count = 0;
        for (GLKPOSITION facePos = Patch->GetFaceList().GetHeadPosition(); facePos;) {
            QMeshFace* thisFace = (QMeshFace*)Patch->GetFaceList().GetNext(facePos);
            thisFace->FieldIndexNumber = count++;
        }
        //std::cout << "Face Count " << count << std::endl;

        patchId++;
    }
}

void meshOperation::initial(PolygenMesh* sourceLayerSet) {

    for (GLKPOSITION posMesh = sourceLayerSet->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* patch = (QMeshPatch*)sourceLayerSet->GetMeshList().GetNext(posMesh);
        int index = 0;
        for (GLKPOSITION Pos = patch->GetEdgeList().GetHeadPosition(); Pos;) {
            QMeshEdge* Edge = (QMeshEdge*)patch->GetEdgeList().GetNext(Pos);
            Edge->SetIndexNo(index);
            index++;
        }
        index = 0;
        for (GLKPOSITION Pos = patch->GetFaceList().GetHeadPosition(); Pos;) {
            QMeshFace* Face = (QMeshFace*)patch->GetFaceList().GetNext(Pos);
            Face->SetIndexNo(index);
            Face->CalPlaneEquation(); // pre-compute the normal of face
            index++;
        }
        index = 0;
        for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

            Node->SetIndexNo(index);
            index++;
        }
    }

    for (GLKPOSITION posMesh = sourceLayerSet->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* each_patch = (QMeshPatch*)sourceLayerSet->GetMeshList().GetNext(posMesh);

        for (GLKPOSITION Pos = each_patch->GetEdgeList().GetHeadPosition(); Pos;) {
            QMeshEdge* Edge = (QMeshEdge*)each_patch->GetEdgeList().GetNext(Pos);

            Edge->inner = true; // clean flag
        }
        for (GLKPOSITION Pos = each_patch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)each_patch->GetNodeList().GetNext(Pos);

            Node->inner = true; // clean flag
        }
        for (GLKPOSITION Pos = each_patch->GetEdgeList().GetHeadPosition(); Pos;) {
            QMeshEdge* Edge = (QMeshEdge*)each_patch->GetEdgeList().GetNext(Pos);

            if (Edge->GetLeftFace() == NULL || Edge->GetRightFace() == NULL) {

                Edge->inner = false;
                Edge->GetStartPoint()->inner = false;
                Edge->GetEndPoint()->inner = false;
            }
        }
    }
}

void meshOperation::offsetMesh(PolygenMesh* sourceLayerSet, PolygenMesh* offsetLayerSet) {

    for (GLKPOSITION posMesh = sourceLayerSet->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* each_patch = (QMeshPatch*)sourceLayerSet->GetMeshList().GetNext(posMesh);

        /*check if it is open mesh*/
        //// node loop
        int temp = 0;
        for (GLKPOSITION Pos = each_patch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)each_patch->GetNodeList().GetNext(Pos);

            if (Node->inner)  temp++;
        }
        if (temp == each_patch->GetNodeNumber()) {
            std::cout << "this is a closed mesh, not considered!" << std::endl;
            return;
        }

        /* get the VERTEX NUM of offset Layer */
        int structuredMesh_NodeNum = each_patch->GetNodeNumber(); // nodes that same as the source mesh itself
        //// node loop
        for (GLKPOSITION Pos = each_patch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)each_patch->GetNodeList().GetNext(Pos);

            if (Node->inner == false) {
                structuredMesh_NodeNum++;
            }
        }

        /* get the VERTEX Table of offset Layer */
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(structuredMesh_NodeNum, 3);
        int V_row_index = 0;
        //// node loop
        for (GLKPOSITION Pos = each_patch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)each_patch->GetNodeList().GetNext(Pos);

            // VERTEX set which will be kept
            double xx, yy, zz;
            Node->GetCoord3D(xx, yy, zz);
            V.row(V_row_index) << xx, yy, zz;
            V_row_index++;
        }
        //// node loop
        for (GLKPOSITION Pos = each_patch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)each_patch->GetNodeList().GetNext(Pos);

            if (Node->inner == false) { //boudnary nodes
                Eigen::Vector3d coord3D_Node, coord3D_offsetNode = { 0.0,0.0,0.0 };
                Node->GetCoord3D(coord3D_Node[0], coord3D_Node[1], coord3D_Node[2]);

                //double connerArea = 0.0;
                //for (GLKPOSITION Pos_neighbor = Node->GetFaceList().GetHeadPosition(); Pos_neighbor;) {
                //    QMeshFace* oneRing_Face = (QMeshFace*)Node->GetFaceList().GetNext(Pos_neighbor);

                //    oneRing_Face->CalArea();
                //    connerArea += oneRing_Face->GetArea();

                //}

                for (GLKPOSITION Pos_neighbor = Node->GetFaceList().GetHeadPosition(); Pos_neighbor;) {
                    QMeshFace* oneRing_Face = (QMeshFace*)Node->GetFaceList().GetNext(Pos_neighbor);

                    Eigen::Vector3d coord3D_FaceCenter_1ring;
                    oneRing_Face->CalCenterPos(coord3D_FaceCenter_1ring[0], coord3D_FaceCenter_1ring[1], coord3D_FaceCenter_1ring[2]);

                    coord3D_offsetNode += offset_value * (coord3D_Node - coord3D_FaceCenter_1ring).normalized();
                    //coord3D_offsetNode += (oneRing_Face->GetArea()/ connerArea) * offset_value * (coord3D_Node - coord3D_FaceCenter_1ring).normalized();

                    //std::cout << "coord3D_offsetNode " << coord3D_offsetNode.transpose() << std::endl;

                }

                V.row(V_row_index) = coord3D_Node + coord3D_offsetNode / (Node->GetFaceList().GetCount());

                //std::cout << "V.row(V_row_index) " << V.row(V_row_index) << std::endl;

                Node->nodeInd_ofOffsetNode_onThisNode = V_row_index;
                V_row_index++;
            }
        }


        /* get the FACE NUM of tight support Layer */
        int structuredMesh_FaceNum = each_patch->GetFaceNumber(); //the kept face number
        // edge loop
        for (GLKPOSITION Pos = each_patch->GetEdgeList().GetHeadPosition(); Pos;) {
            QMeshEdge* Edge = (QMeshEdge*)each_patch->GetEdgeList().GetNext(Pos);

            if (Edge->inner == false)
                structuredMesh_FaceNum += 2; // 2 means two triangle will generate
            //             |-----|
            //             |  /  |
            //             |-----|
        }

        /* get the FACE Table of tight support Layer */
        Eigen::MatrixXi F = Eigen::MatrixXi::Zero(structuredMesh_FaceNum, 3);
        // face loop
        int F_row_index = 0;
        for (GLKPOSITION Pos = each_patch->GetFaceList().GetHeadPosition(); Pos;) {
            QMeshFace* Face = (QMeshFace*)each_patch->GetFaceList().GetNext(Pos);

            for (int k = 0; k < 3; k++) {
                F(F_row_index, k) = Face->GetNodeRecordPtr(k)->GetIndexNo();
            }
            F_row_index++;
        }
        // edge loop
        for (GLKPOSITION Pos = each_patch->GetEdgeList().GetHeadPosition(); Pos;) {
            QMeshEdge* Edge = (QMeshEdge*)each_patch->GetEdgeList().GetNext(Pos);

            if (Edge->inner == false) {

                QMeshNode* node1 = Edge->GetStartPoint();
                QMeshNode* node2 = Edge->GetEndPoint();

                // the 1st added face
                F(F_row_index, 0) = node1->GetIndexNo();
                F(F_row_index, 1) = node2->nodeInd_ofOffsetNode_onThisNode;
                F(F_row_index, 2) = node2->GetIndexNo();
                F_row_index++;

                // the 2nd added face
                F(F_row_index, 0) = node1->GetIndexNo();
                F(F_row_index, 1) = node1->nodeInd_ofOffsetNode_onThisNode;
                F(F_row_index, 2) = node2->nodeInd_ofOffsetNode_onThisNode;
                F_row_index++;

            }
        }

        // build new mesh from vertex table and face table
        float* nodeTable;
        nodeTable = (float*)malloc(sizeof(float) * V.rows() * 3);
        for (int j = 0; j < V.rows(); j++) {
            for (int i = 0; i < 3; i++) nodeTable[j * 3 + i] = (float)V(j, i);
        }
        unsigned int* faceTable;
        faceTable = (unsigned int*)malloc(sizeof(unsigned int) * F.rows() * 3);
        for (int j = 0; j < F.rows(); j++) {
            for (int i = 0; i < 3; i++) faceTable[j * 3 + i] = F(j, i);
        }

        QMeshPatch* offsetMesh = new QMeshPatch;
        offsetMesh->SetIndexNo(offsetLayerSet->GetMeshList().GetCount()); //index begin from 0
        offsetLayerSet->GetMeshList().AddTail(offsetMesh);
        offsetMesh->constructionFromVerFaceTable(V.rows(), nodeTable, F.rows(), faceTable);
        offsetMesh->patchName = each_patch->patchName;

    }
}