#include "toolpathGeneration_stripe.h"

void toolpathGeneration_stripe::_initial_index(QMeshPatch* patch) {

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

    for (GLKPOSITION Pos = patch->GetEdgeList().GetHeadPosition(); Pos;) {
        QMeshEdge* Edge = (QMeshEdge*)patch->GetEdgeList().GetNext(Pos);

        Edge->inner = true; // clean flag
    }
    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        Node->inner = true; // clean flag
    }
    for (GLKPOSITION Pos = patch->GetEdgeList().GetHeadPosition(); Pos;) {
        QMeshEdge* Edge = (QMeshEdge*)patch->GetEdgeList().GetNext(Pos);

        if (Edge->GetLeftFace() == NULL || Edge->GetRightFace() == NULL) {

            Edge->inner = false;
            Edge->GetStartPoint()->inner = false;
            Edge->GetEndPoint()->inner = false;
        }
    }
}

void toolpathGeneration_stripe::generate_toolPath_marchingSquare() {

    for (GLKPOSITION PatchPos = m_Slices->GetMeshList().GetHeadPosition(); PatchPos;) {
        QMeshPatch* layer = (QMeshPatch*)m_Slices->GetMeshList().GetNext(PatchPos);

        this->_initial_index(layer);
        this->_mark_boundary_face(layer);

        //build a node pos table of layer
        std::vector<Eigen::Vector3d> node_pos;
        for (GLKPOSITION Pos = layer->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)layer->GetNodeList().GetNext(Pos);

            Eigen::Vector3d _pos;
            Node->GetCoord3D(_pos[0], _pos[1], _pos[2]);
            node_pos.push_back(_pos);
        }

        std::vector<Eigen::Vector3d> vertices;
        std::vector<Eigen::Vector2i> edges;
        // list of per-edge(curved layer) indices pointing to the vertices(iso-node) list
        // In other words, the iso-nodes attached on the each edge of curved layer 
        std::vector<std::vector<size_t>> polylineIndices(layer->GetEdgeNumber());

        for (GLKPOSITION Pos = layer->GetFaceList().GetHeadPosition(); Pos;) {
            QMeshFace* Face = (QMeshFace*)layer->GetFaceList().GetNext(Pos);

            if (Face->stripeIndices[0] != 0
                || Face->stripeIndices[1] != 0
                || Face->fieldIndices != 0
                || Face->isBoundary
                ) continue; // singularities are ignored in this function

            std::vector<std::vector<size_t>> edgeIndices(3);

            Eigen::Vector3i fv(
                Face->GetNodeRecordPtr(0)->GetIndexNo(),
                Face->GetNodeRecordPtr(1)->GetIndexNo(),
                Face->GetNodeRecordPtr(2)->GetIndexNo());
            Eigen::Vector3d fs(
                Face->GetTextureCoordu(0),
                Face->GetTextureCoordu(1),
                Face->GetTextureCoordu(2));
            Eigen::Vector3i fe(
                _getOrderedEdgeIndexOfFace(Face, 0),
                _getOrderedEdgeIndexOfFace(Face, 1),
                _getOrderedEdgeIndexOfFace(Face, 2));

            /*std::cout
                << "e0s: " << Face->GetEdgeRecordPtr(1)->GetStartPoint()->GetIndexNo()
                << " e0e: " << Face->GetEdgeRecordPtr(1)->GetEndPoint()->GetIndexNo()
                << " e1s: " << Face->GetEdgeRecordPtr(2)->GetStartPoint()->GetIndexNo()
                << " e1e: " << Face->GetEdgeRecordPtr(2)->GetEndPoint()->GetIndexNo()
                << " e2s: " << Face->GetEdgeRecordPtr(3)->GetStartPoint()->GetIndexNo()
                << " e2e: " << Face->GetEdgeRecordPtr(3)->GetEndPoint()->GetIndexNo()
                << std::endl;
            std::cout << "fv: " << fv.transpose() << std::endl;
            std::cout << "fs: " << fs.transpose() << std::endl;
            std::cout << "fe: " << fe.transpose() << std::endl;*/

            for (int _id = 0; _id < 3; _id++) {
                // list all the crossings along edge
                size_t v0 = fv[_id]; size_t v1 = fv[(_id + 1) % 3];
                double s0 = fs[_id]; double s1 = fs[(_id + 1) % 3];
                size_t e0 = fe[_id];

                std::vector<double> isoPoints = _calculateIsoPoints(s0, s1);

                if (polylineIndices[e0].empty()) {
                    for (const auto& bary : isoPoints) {
                        polylineIndices[e0].push_back(vertices.size());
                        Eigen::Vector3d newVertex = bary * node_pos[v1] + (1 - bary) * node_pos[v0];
                        vertices.push_back(newVertex);
                    }

                    edgeIndices[_id] = polylineIndices[e0];
                }
                else {
                    std::vector<size_t> pi_e0_reverse = polylineIndices[e0];
                    std::reverse(pi_e0_reverse.begin(), pi_e0_reverse.end());
                    edgeIndices[_id].insert(edgeIndices[_id].end(), pi_e0_reverse.begin(), pi_e0_reverse.end());
                }

            }

            std::vector<Eigen::Vector2i> matchings = this->_matchCrossings(edgeIndices);
            edges.insert(edges.end(), matchings.begin(), matchings.end());

        }

        if (vertices.size() == 0 || edges.size() == 0) continue; // skip the empty toolpath

        //build a new path to install the nodes and edges
        QMeshPatch* singlePath = new QMeshPatch;
        singlePath->SetIndexNo(layer->GetIndexNo());
        m_Waypoints->GetMeshList().AddTail(singlePath);
        singlePath->patchName = layer->patchName;// set the path name to be the same as layer name
        _buildPath(vertices, edges, singlePath);

        //record the iso-node on which edge of curved layer 
        //input: node list of path / edge list of layer / polylineIndices
        _build_connection_between_isoNode_with_relatedEdge(singlePath, layer, polylineIndices);
        //set the isoNodeNormal for isoNode
        for (GLKPOSITION Pos = singlePath->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* isoNode = (QMeshNode*)singlePath->GetNodeList().GetNext(Pos);

            QMeshEdge* edge_attachedOn = isoNode->relatedLayerEdge;

            Eigen::Vector3d n1, n2, n3;
            edge_attachedOn->GetLeftFace()->GetNormal(n1[0], n1[1], n1[2]);
            edge_attachedOn->GetRightFace()->GetNormal(n2[0], n2[1], n2[2]);
            n3 = (n1 + n2).normalized();
            isoNode->SetNormal(n3[0], n3[1], n3[2]);
        }

        //handle the singular face for better connection
        _build_edges_in_singular_region(singlePath, layer);

        //record which facet of isoLayer  the edges of toolPath is on (edge(path) -> face(layer))
        _record_edgeOfPath_on_which_faceOFLayer(singlePath);
    }
}

size_t toolpathGeneration_stripe::_getOrderedEdgeIndexOfFace(QMeshFace* Face, size_t _id_edge) {

    size_t aNode_idx = Face->GetNodeRecordPtr(_id_edge)->GetIndexNo();
    size_t bNode_idx = Face->GetNodeRecordPtr((_id_edge + 1) % 3)->GetIndexNo();

    size_t edge_idx = -1;

    for (int i = 0; i < 3; i++) {

        size_t aNode_idx_edge = Face->GetEdgeRecordPtr(i + 1)->GetStartPoint()->GetIndexNo();
        size_t bNode_idx_edge = Face->GetEdgeRecordPtr(i + 1)->GetEndPoint()->GetIndexNo();

        if (((aNode_idx == aNode_idx_edge) && (bNode_idx == bNode_idx_edge)) ||
            ((aNode_idx == bNode_idx_edge) && (bNode_idx == aNode_idx_edge))) {
            edge_idx = Face->GetEdgeRecordPtr(i + 1)->GetIndexNo();
        }

    }

    if (edge_idx == -1) std::cout << "Error: @_getOrderedEdgeIndexOfFace" << std::endl;

    return edge_idx;
}

std::vector<double> toolpathGeneration_stripe::_calculateIsoPoints(double s0, double s1) {
    std::vector<double> isoPoints;

    // Generate the range of points
    for (int i = std::ceil(std::min(s0, s1)); i <= std::floor(std::max(s0, s1)); ++i) {
        isoPoints.push_back(i);
    }

    // Normalize the points
    for (double& ips : isoPoints) {
        ips = (ips - std::min(s0, s1)) / std::abs(s1 - s0);
    }

    // Adjust points if s0 > s1
    if (s0 > s1) {
        for (double& ips : isoPoints) {
            ips = 1 - ips;
        }
        std::reverse(isoPoints.begin(), isoPoints.end());
    }

    return isoPoints;
}

std::vector<Eigen::Vector2i> toolpathGeneration_stripe::_matchCrossings(const std::vector<std::vector<size_t>>& crossings) {
    assert(crossings.size() == 3);

    size_t idxIJ = 2;
    if (crossings[0].size() >= crossings[1].size() && crossings[0].size() >= crossings[2].size()) {
        idxIJ = 0;
    }
    else if (crossings[1].size() >= crossings[2].size() && crossings[1].size() >= crossings[0].size()) {
        idxIJ = 1;
    }

    size_t idxJK = (idxIJ + 1) % 3;
    size_t idxKI = (idxIJ + 2) % 3;

    const std::vector<size_t>& IJ = crossings[idxIJ];
    const std::vector<size_t>& JK = crossings[idxJK];
    const std::vector<size_t>& KI = crossings[idxKI];

    size_t nIJ = IJ.size();
    size_t nJK = JK.size();
    size_t nKI = KI.size();

    assert(nIJ >= nJK && nIJ >= nKI);
    assert(nIJ <= nJK + nKI);
    assert((nIJ + nJK + nKI) % 2 == 0);

    std::vector<Eigen::Vector2i> matchings;

    if (nIJ == nJK + nKI) {
        // Case 1: all edges intersecting ijk cross a common edge ij
        // match IJ with KI
        for (size_t m = 0; m < nKI; ++m) {
            matchings.push_back(Eigen::Vector2i(IJ[m], KI[nKI - m - 1]));
        }
        // match IJ with JK
        for (size_t m = 0; m < nJK; ++m) {
            matchings.push_back(Eigen::Vector2i(IJ[nKI + m], JK[nJK - m - 1]));
        }
    }
    else {
        // Case 2: there is no common edge
        size_t nRemainingCrossings = (nIJ + nJK + nKI) / 2;
        size_t m = 0;
        while (nRemainingCrossings > nJK) {
            matchings.push_back(Eigen::Vector2i(IJ[m], KI[nKI - m - 1]));
            ++m;
            --nRemainingCrossings;
        }

        size_t l = 0;
        while (nRemainingCrossings > nKI - m) {
            matchings.push_back(Eigen::Vector2i(IJ[nIJ - 1 - l], JK[l]));
            --nRemainingCrossings;
            ++l;
        }

        size_t p = 0;
        while (nRemainingCrossings > 0) {
            matchings.push_back(Eigen::Vector2i(JK[nJK - 1 - p], KI[p]));
            ++p;
            --nRemainingCrossings;
        }
    }

    return matchings;
}

void toolpathGeneration_stripe::_buildPath(
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<Eigen::Vector2i>& edges,
    QMeshPatch* singlePath) {
    std::vector<QMeshNode*> nodeList;

    // Create nodes for each vertex
    for (size_t i = 0; i < vertices.size(); ++i) {
        const Eigen::Vector3d& vertex = vertices[i];

        // Create a new QMeshNode
        QMeshNode* isoNode = new QMeshNode;
        isoNode->SetMeshPatchPtr(singlePath);
        isoNode->SetCoord3D(vertex[0], vertex[1], vertex[2]);
        isoNode->SetIndexNo(static_cast<int>(singlePath->GetNodeList().GetCount()));

        nodeList.push_back(isoNode);
        singlePath->GetNodeList().AddTail(isoNode);
    }

    // Create edges between the nodes
    for (const auto& edge : edges) {
        int index1 = edge[0];
        int index2 = edge[1];

        QMeshNode* startNode = nodeList[index1];
        QMeshNode* endNode = nodeList[index2];

        // Create a new QMeshEdge
        QMeshEdge* isoEdge = new QMeshEdge;
        isoEdge->SetStartPoint(startNode);
        isoEdge->SetEndPoint(endNode);
        isoEdge->SetMeshPatchPtr(singlePath);

        isoEdge->SetIndexNo(singlePath->GetEdgeList().GetCount());

        (startNode->GetEdgeList()).AddTail(isoEdge);
        (endNode->GetEdgeList()).AddTail(isoEdge);
        singlePath->GetEdgeList().AddTail(isoEdge);

    }
}

void toolpathGeneration_stripe::_mark_boundary_face(QMeshPatch* surfaceMesh) {

    for (GLKPOSITION Pos = surfaceMesh->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)surfaceMesh->GetNodeList().GetNext(Pos);

        //initialize
        Node->isBoundaryNode = false;
    }

    //get boundary node and install
	for (GLKPOSITION Pos = surfaceMesh->GetEdgeList().GetHeadPosition(); Pos;) {
		QMeshEdge *Edge = (QMeshEdge*)surfaceMesh->GetEdgeList().GetNext(Pos);
		if (Edge->IsBoundaryEdge()) {
			Edge->GetStartPoint()->isBoundaryNode = true;
            Edge->GetEndPoint()->isBoundaryNode = true;
		}
	}
    for (GLKPOSITION Pos = surfaceMesh->GetFaceList().GetHeadPosition(); Pos;) {
        QMeshFace* Face = (QMeshFace*)surfaceMesh->GetFaceList().GetNext(Pos);

        bool isBoundaryFace = false;
        for (int i = 0; i < 3; i++) {
            if (Face->GetNodeRecordPtr(i)->isBoundaryNode){
                isBoundaryFace = true;
                break;
            }
        }

        Face->isBoundary = isBoundaryFace;
    }
}

void toolpathGeneration_stripe::_build_connection_between_isoNode_with_relatedEdge(
    QMeshPatch* singlePath, QMeshPatch* layer, const std::vector<std::vector<size_t>>& polylineIndices) {

    // Create iso-node list for singlePath
    std::vector<QMeshNode*> isoNodeList;
    for (GLKPOSITION Pos = singlePath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)singlePath->GetNodeList().GetNext(Pos);

        isoNodeList.push_back(Node);
    }
    // Create edge list for layer
    std::vector<QMeshEdge*> edgeList;
    for (GLKPOSITION Pos = layer->GetEdgeList().GetHeadPosition(); Pos;) {
        QMeshEdge* Edge = (QMeshEdge*)layer->GetEdgeList().GetNext(Pos);

        edgeList.push_back(Edge);
    }

    // Ensure polylineIndices is valid and corresponds to edgeList size
    assert(polylineIndices.size() == edgeList.size());

    // Iterate over polylineIndices to set up relationships
    for (size_t i = 0; i < polylineIndices.size(); ++i) {
        QMeshEdge* edge = edgeList[i];
        const std::vector<size_t>& nodeIndices = polylineIndices[i];

        // Vector to store isoNodes for this edge
        std::vector<QMeshNode*> installedIsoNode_layerEdge;

        for (size_t index : nodeIndices) {
            QMeshNode* isoNode = isoNodeList[index];

            // Set the related edge for the isoNode
            isoNode->relatedLayerEdge = edge;

            // Add the isoNode to the edge's installedIsoNode_layerEdge vector
            installedIsoNode_layerEdge.push_back(isoNode);
        }

        // Set the installedIsoNode_layerEdge vector for the edge
        edge->installedIsoNode_layerEdge = installedIsoNode_layerEdge;
    }

}

void toolpathGeneration_stripe::_build_edges_in_singular_region(QMeshPatch* singlePath, QMeshPatch* layer) {

    //mark the endNodes of the singlePath
    for (GLKPOSITION Pos = singlePath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* isoNode = (QMeshNode*)singlePath->GetNodeList().GetNext(Pos);

        isoNode->isBoundaryNode = false;
        if (isoNode->GetEdgeNumber() < 2) isoNode->isBoundaryNode = true;
    }

    for (GLKPOSITION Pos = layer->GetFaceList().GetHeadPosition(); Pos;) {
        QMeshFace* Face = (QMeshFace*)layer->GetFaceList().GetNext(Pos);

        if (Face->isBoundary) continue; // boundary face are ignored in this function

        // singularities are considered in this function
        if (Face->stripeIndices[0] != 0 || Face->stripeIndices[1] != 0 || Face->fieldIndices != 0) {

            //collect all of the end nodes whose flag  isoNode->isBoundaryNode = true;
            std::vector<std::vector<QMeshNode*>> boundary_isoNodes_on_Edges_of_one_Face(3);

            for (int i = 0; i < 3; i++) {

                std::vector<QMeshNode*> _isoNodes_on_oneEdges
                    = Face->GetEdgeRecordPtr(i + 1)->installedIsoNode_layerEdge;

                for (int j = 0; j < _isoNodes_on_oneEdges.size(); j++) {

                    if (_isoNodes_on_oneEdges[j]->isBoundaryNode)
                        boundary_isoNodes_on_Edges_of_one_Face[i].push_back(_isoNodes_on_oneEdges[j]);
                }
            }

            // Check for all combinations of edges with iso-nodes
            std::vector<int> edges_with_isoNodes;
            for (int i = 0; i < 3; i++) {
                if (!boundary_isoNodes_on_Edges_of_one_Face[i].empty()) {
                    edges_with_isoNodes.push_back(i);
                }
            }

            // If there are exactly two edges with iso-nodes, create an edge between the iso-nodes
            if (edges_with_isoNodes.size() == 2) {
                int edge1 = edges_with_isoNodes[0];
                int edge2 = edges_with_isoNodes[1];

                QMeshNode* startNode = boundary_isoNodes_on_Edges_of_one_Face[edge1][0];
                QMeshNode* endNode = boundary_isoNodes_on_Edges_of_one_Face[edge2][0];

                // Create a new QMeshEdge
                QMeshEdge* isoEdge = new QMeshEdge;
                isoEdge->SetStartPoint(startNode);
                isoEdge->SetEndPoint(endNode);
                isoEdge->SetMeshPatchPtr(singlePath);

                isoEdge->SetIndexNo(singlePath->GetEdgeList().GetCount());

                (startNode->GetEdgeList()).AddTail(isoEdge);
                (endNode->GetEdgeList()).AddTail(isoEdge);
                singlePath->GetEdgeList().AddTail(isoEdge);
            }

            // New case: If each edge of a triangle face owns iso-Node, connect the nearest iso-nodes
            if (edges_with_isoNodes.size() == 3) {
                // Find the nearest iso-nodes between the edges
                double minDistance = DBL_MAX;
                QMeshNode* nearestNode1 = nullptr;
                QMeshNode* nearestNode2 = nullptr;

                for (int i = 0; i < 3; i++) {
                    for (int j = i + 1; j < 3; j++) {
                        for (auto node1 : boundary_isoNodes_on_Edges_of_one_Face[i]) {
                            for (auto node2 : boundary_isoNodes_on_Edges_of_one_Face[j]) {

                                Eigen::Vector3d p1, p2;
                                node1->GetCoord3D(p1[0], p1[1], p1[2]);
                                node2->GetCoord3D(p2[0], p2[1], p2[2]);

                                double distance = (p1 - p2).norm();
                                if (distance < minDistance) {
                                    minDistance = distance;
                                    nearestNode1 = node1;
                                    nearestNode2 = node2;
                                }
                            }
                        }
                    }
                }

                if (nearestNode1 && nearestNode2) {
                    // Create a new QMeshEdge
                    QMeshEdge* isoEdge = new QMeshEdge;
                    isoEdge->SetStartPoint(nearestNode1);
                    isoEdge->SetEndPoint(nearestNode2);
                    isoEdge->SetMeshPatchPtr(singlePath);

                    isoEdge->SetIndexNo(singlePath->GetEdgeList().GetCount());

                    (nearestNode1->GetEdgeList()).AddTail(isoEdge);
                    (nearestNode2->GetEdgeList()).AddTail(isoEdge);
                    singlePath->GetEdgeList().AddTail(isoEdge);
                }
            }
        }
    }

    //remark the endNodes of the singlePath
    for (GLKPOSITION Pos = singlePath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* isoNode = (QMeshNode*)singlePath->GetNodeList().GetNext(Pos);

        isoNode->isBoundaryNode = false;
        if (isoNode->GetEdgeNumber() < 2) isoNode->isBoundaryNode = true;
    }
}

void toolpathGeneration_stripe::_record_edgeOfPath_on_which_faceOFLayer(QMeshPatch* singlePath) {

    for (GLKPOSITION Pos = singlePath->GetEdgeList().GetHeadPosition(); Pos;) {
        QMeshEdge* Edge = (QMeshEdge*)singlePath->GetEdgeList().GetNext(Pos);
        
        QMeshNode* sNode = Edge->GetStartPoint();
        QMeshNode* eNode = Edge->GetEndPoint();

        QMeshEdge* layerEdge_of_sNode = sNode->relatedLayerEdge;
        QMeshEdge* layerEdge_of_eNode = eNode->relatedLayerEdge;

        QMeshFace* lFace_layerEdge_of_sNode = layerEdge_of_sNode->GetLeftFace();
        QMeshFace* rFace_layerEdge_of_sNode = layerEdge_of_sNode->GetRightFace();

        QMeshFace* lFace_layerEdge_of_eNode = layerEdge_of_eNode->GetLeftFace();
        QMeshFace* rFace_layerEdge_of_eNode = layerEdge_of_eNode->GetRightFace();

        if(!lFace_layerEdge_of_sNode || !rFace_layerEdge_of_sNode || !lFace_layerEdge_of_eNode || !rFace_layerEdge_of_eNode)
            std::cout << "Error: l/r face is NULL, @_record_edgeOfPath_on_which_faceOFLayer" << std::endl;

        //find the same face among above 4 faces
        QMeshFace* sameFace = nullptr;
        if (lFace_layerEdge_of_sNode == lFace_layerEdge_of_eNode || lFace_layerEdge_of_sNode == rFace_layerEdge_of_eNode) {
            sameFace = lFace_layerEdge_of_sNode;
        }
        else if (rFace_layerEdge_of_sNode == lFace_layerEdge_of_eNode || rFace_layerEdge_of_sNode == rFace_layerEdge_of_eNode) {
            sameFace = rFace_layerEdge_of_sNode;
        }

        if (sameFace) {
            Edge->relatedLayerFace = sameFace;
        }
        else {
            std::cout << "Error: there is no relatedLayerFace for the Edge of Path, @_record_edgeOfPath_on_which_faceOFLayer" << std::endl;
        }
    }

}

/*Main functions for the post processing of strip toolpath*/

void toolpathGeneration_stripe::postProcessing_stripPath() {

    for (GLKPOSITION posMesh = m_Waypoints->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* _toolpath = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(posMesh);

        //classify the each toolpath chains // adding the seg_Idx
        this->_classify_chains(_toolpath);

        //mark(remove) small chains
        this->_mark_short_chains(_toolpath, 10.0); //100mm discard threshold

        //connection improvment
        QMeshNode* startNode_4_oneStroke_path = this->_connect_shortGap(_toolpath, 10.0); //shortGap threshold

        //trace the strip path into one-stroke path
        this->_trace_2_oneStroke_path(_toolpath, startNode_4_oneStroke_path);

        //smooth chains
        this->_smooth_chains(_toolpath, 20, 0.5); //smooth time and ratio

        //auxiliary info attachment
        //dist and jump flag
        this->_cal_jumpFlag_4_CCF_toolpath(_toolpath);

        //cut info of nodes
        this->_cal_cutInfo_4_CCF_toolpath(_toolpath);

        //tangent direction info of nodes
        this->_cal_tangentDir_4_CCF_toolpath(_toolpath);
    }
}

void toolpathGeneration_stripe::_classify_chains(QMeshPatch* _toolpath) {

    //pick the endPnts in the _toolpath (open curve)
    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);

        Node->is_endPnt = false; //initialize
        Node->is_Processed = false;
        Node->is_shortChain = false;

        if (Node->GetEdgeNumber() == 1) Node->is_endPnt = true;
        if (Node->GetEdgeNumber() == 0) std::cout << "Error: @_classify_chains" << std::endl;
    }
    //clean the visited Flag for edges
    for (GLKPOSITION Pos = _toolpath->GetEdgeList().GetHeadPosition(); Pos;) {
        QMeshEdge* Edge = (QMeshEdge*)_toolpath->GetEdgeList().GetNext(Pos);
        Edge->visited = false;
        Edge->is_shortChain = false;
    }

    /*trace each open chain*/
    int seg_Idx = 0;
    QMeshNode* start_Node = NULL;
    //iterational tracing open chain
    do {

        start_Node = this->_pickUp_one_startNode(_toolpath, seg_Idx);
        this->_tracing_one_path(_toolpath, start_Node, seg_Idx);
        seg_Idx++;

    } while (!this->_allDetced_4_openChains(_toolpath));

    //iterational tracing close chain
    do {

        start_Node = this->_pickUp_one_startNode_closeLoop(_toolpath, seg_Idx);
        if (start_Node == NULL) break; //no clsed loop is found
        this->_tracing_one_path_closeLoop(_toolpath, start_Node, seg_Idx);
        seg_Idx++;

    } while (!this->_allDetced_4_openChains_closeLoop(_toolpath));

}

QMeshNode* toolpathGeneration_stripe::_pickUp_one_startNode(QMeshPatch* _toolpath, int seg_Idx) {

    QMeshNode* start_Node = NULL;

    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);

        if (Node->is_endPnt && !Node->is_Processed) {
            start_Node = Node;
            start_Node->is_Processed = true;
            start_Node->seg_Idx = seg_Idx;
            break;
        }
    }

    if (start_Node == NULL)
        std::cout << "start node is NULL @ _pickUp_one_startNode" << std::endl;

    return start_Node;
}

void toolpathGeneration_stripe::_tracing_one_path(QMeshPatch* _toolpath, QMeshNode* start_Node, int seg_Idx) {

    QMeshNode* curtNode = start_Node;

    do {

        for (GLKPOSITION Pos = curtNode->GetEdgeList().GetHeadPosition(); Pos;) {
            QMeshEdge* neighEdge = (QMeshEdge*)curtNode->GetEdgeList().GetNext(Pos);

            if (neighEdge->visited) continue;

            QMeshNode* neighNode = neighEdge->GetStartPoint();
            if (neighNode == curtNode) neighNode = neighEdge->GetEndPoint();

            neighEdge->visited = true;
            neighEdge->seg_Idx = seg_Idx;
            curtNode = neighNode;
            curtNode->is_Processed = true;
            curtNode->seg_Idx = seg_Idx;
            break;
        }
    } while (!curtNode->is_endPnt);

}

bool toolpathGeneration_stripe::_allDetced_4_openChains(QMeshPatch* _toolpath) {

    bool allVisited = true;

    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);

        if (Node->is_endPnt) {

            if (!Node->is_Processed) {
                allVisited = false;
                break;
            }
        }
    }

    return allVisited;
}

QMeshNode* toolpathGeneration_stripe::_pickUp_one_startNode_closeLoop(QMeshPatch* _toolpath, int seg_Idx) {

    QMeshNode* start_Node = NULL;

    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);

        if (!Node->is_Processed) {
            start_Node = Node;
            start_Node->is_Processed = true;
            start_Node->seg_Idx = seg_Idx;
            start_Node->is_closeLoop = true;
            break;
        }
    }

    if (start_Node == NULL)
        std::cout << "start node is NULL @ _pickUp_one_startNode_closeLoop, maybe because no clsed loop is found, ignore" << std::endl;

    return start_Node;
}

void toolpathGeneration_stripe::_tracing_one_path_closeLoop(QMeshPatch* _toolpath, QMeshNode* start_Node, int seg_Idx) {

    QMeshNode* curtNode = start_Node;

    do {

        for (GLKPOSITION Pos = curtNode->GetEdgeList().GetHeadPosition(); Pos;) {
            QMeshEdge* neighEdge = (QMeshEdge*)curtNode->GetEdgeList().GetNext(Pos);

            if (neighEdge->visited) continue;

            QMeshNode* neighNode = neighEdge->GetStartPoint();
            if (neighNode == curtNode) neighNode = neighEdge->GetEndPoint();

            neighEdge->visited = true;
            neighEdge->seg_Idx = seg_Idx;
            neighEdge->is_closeLoop = true;
            curtNode = neighNode;
            curtNode->is_Processed = true;
            curtNode->seg_Idx = seg_Idx;
            curtNode->is_closeLoop = true;
            break;
        }
    } while (curtNode != start_Node);

}

bool toolpathGeneration_stripe::_allDetced_4_openChains_closeLoop(QMeshPatch* _toolpath) {

    bool allVisited = true;

    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);

        if (!Node->is_Processed) {
            allVisited = false;
            break;
        }
    }

    return allVisited;
}

void toolpathGeneration_stripe::_mark_short_chains(QMeshPatch* _toolpath, double cutLength) {

    for (GLKPOSITION Pos = _toolpath->GetEdgeList().GetHeadPosition(); Pos;) {
        QMeshEdge* Edge = (QMeshEdge*)_toolpath->GetEdgeList().GetNext(Pos);

        if (Edge->seg_Idx < 0)
            std::cout << "Error:  @ _mark_short_chains, there are still some unProcessed Edge." << std::endl;
    }

    //record the chains length
    std::vector<double> data;
    for (GLKPOSITION Pos = _toolpath->GetEdgeList().GetHeadPosition(); Pos;) {
        QMeshEdge* Edge = (QMeshEdge*)_toolpath->GetEdgeList().GetNext(Pos);

        int index = Edge->seg_Idx;

        if (index >= data.size()) {
            // If index is out of range, resize the vector and set the value
            data.resize(index + 1, 0.0);
            data[index] = 0.0;
        }
        else {
            // If index is within range, update the value
            data[index] += Edge->CalLength();
        }
    }

    //print out and collect shortChain idx
    std::vector<int> shortChain_idx_Set;
    for (size_t i = 0; i < data.size(); i++) {
        if (data[i] != 0.0 && data[i] < cutLength) {
            //std::cout << "Index " << i << ": " << data[i] << std::endl;
            shortChain_idx_Set.push_back(i);
        }
    }

    std::cout << "The number of rest strip path is " << data.size() - shortChain_idx_Set.size() << std::endl;


    for (GLKPOSITION Pos = _toolpath->GetEdgeList().GetHeadPosition(); Pos;) {
        QMeshEdge* Edge = (QMeshEdge*)_toolpath->GetEdgeList().GetNext(Pos);

        int index = Edge->seg_Idx;
        auto it = std::find(shortChain_idx_Set.begin(), shortChain_idx_Set.end(), index);

        if (it != shortChain_idx_Set.end()) {
            Edge->is_shortChain = true;
            Edge->GetStartPoint()->is_shortChain = true;
            Edge->GetEndPoint()->is_shortChain = true;
        }
    }

}

void toolpathGeneration_stripe::_smooth_chains(QMeshPatch* _toolpath, int loop, double ratio) {

    for (size_t i = 0; i < loop; ++i) {
        for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);

            if (Node->is_shortChain || Node->is_endPnt) continue;

            Eigen::Vector3d _pm;
            Node->GetCoord3D(_pm[0], _pm[1], _pm[2]);

            //only smooth mid points
            if (Node->GetEdgeNumber() != 2) {
                std::cout << "Error: @_smooth_chains 1" << std::endl;
                //Node->is_special_SinglePt = true;
                return;
            }

            std::vector<QMeshNode*> neighNodes;
            for (GLKPOSITION Pos = Node->GetEdgeList().GetHeadPosition(); Pos;) {
                QMeshEdge* neighEdge = (QMeshEdge*)Node->GetEdgeList().GetNext(Pos);

                QMeshNode* neighNode = neighEdge->GetStartPoint();
                if (neighNode == Node) neighNode = neighEdge->GetEndPoint();

                neighNodes.push_back(neighNode);
            }

            if (neighNodes.size() != 2) std::cout << "Error: @_smooth_chains 2" << std::endl;

            Eigen::Vector3d _p1, _p2;
            neighNodes[0]->GetCoord3D(_p1[0], _p1[1], _p1[2]);
            neighNodes[1]->GetCoord3D(_p2[0], _p2[1], _p2[2]);

            Eigen::Vector3d pp;
            pp = ratio * (0.5 * (_p1 + _p2)) + (1 - ratio) * _pm;
            Node->SetCoord3D(pp[0], pp[1], pp[2]);

        }
    }
}

QMeshNode* toolpathGeneration_stripe::_connect_shortGap(QMeshPatch* _toolpath, double gap_dist_threshold) {

    std::vector<std::pair<double, std::vector<QMeshNode*>>> endPnts_linkages;

    /*Open chain part*/
    //build a table contains the endPnt of opened longChain
    std::vector<QMeshNode*> endPnts;
    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);

        if (!Node->is_closeLoop && !Node->is_shortChain && Node->is_endPnt)
            endPnts.push_back(Node);
    }

    for (size_t i = 0; i < endPnts.size(); ++i) {

        for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);
            Node->visited = false; //clean flag
        }

        std::vector<QMeshNode*> one_solusion;
        QMeshNode* startNode = endPnts[i];	startNode->visited = true;
        one_solusion.push_back(startNode);

        do {

            QMeshNode* another_endPnt_sameChain = this->_get_another_endPnt_sameChain(startNode, endPnts);
            another_endPnt_sameChain->visited = true;
            one_solusion.push_back(another_endPnt_sameChain);

            if (this->_all_visited(endPnts)) break;

            QMeshNode* next_nearest_startNode = this->_find_nearest_or_second_nearest_startNode(
                another_endPnt_sameChain, endPnts, true);
            // QMeshNode* next_nearest_startNode = this->_find_next_nearest_startNode(another_endPnt_sameChain, endPnts);
            next_nearest_startNode->visited = true;
            one_solusion.push_back(next_nearest_startNode);

            startNode = next_nearest_startNode;

        } while (!this->_all_visited(endPnts));

        //calculate the gap length
        double gap_dist = 0.0;
        int before_seg_Idx = one_solusion[0]->seg_Idx;
        for (int j = 1; j < one_solusion.size(); j++) {

            if (before_seg_Idx != one_solusion[j]->seg_Idx) {

                Eigen::Vector3d p1, p2;
                one_solusion[j]->GetCoord3D(p1[0], p1[1], p1[2]);
                one_solusion[j - 1]->GetCoord3D(p2[0], p2[1], p2[2]);

                gap_dist += (p1 - p2).norm();
            }
            before_seg_Idx = one_solusion[j]->seg_Idx;
        }

        //build one row of the container 
        std::pair<double, std::vector<QMeshNode*>> _element = std::make_pair(gap_dist, one_solusion);
        endPnts_linkages.push_back(_element);
    }

    //for (const auto& pair : endPnts_linkages) {
    //	std::cout << "Distance: " << pair.first << std::endl;
    //	std::cout << "QMeshNode pointers:" << std::endl;
    //	for (QMeshNode* node : pair.second) {
    //		std::cout << node->seg_Idx << " "; // Print the pointer value
    //	}
    //	std::cout << std::endl;
    //}

    // Find the pair with the smallest pair.first
    double minDist = std::numeric_limits<double>::max();
    std::vector<QMeshNode*> minPairSecond;
    for (const auto& pair : endPnts_linkages) {
        if (pair.first < minDist) {
            minDist = pair.first;
            minPairSecond = pair.second;
        }
    }

    //// Output the pair with the smallest pair.first
    //std::cout << "Smallest pair.first: " << minDist << std::endl;
    //std::cout << "Corresponding pair.second:" << std::endl;
    //for (QMeshNode* node : minPairSecond) {
    //	std::cout << node->GetIndexNo()<< " "; // Print the pointer value
    //}std::cout << std::endl;

    //build edges to link the chains (within gap_dist_threshold)
    int before_seg_Idx = minPairSecond[0]->seg_Idx;


    std::cout << "before_seg_Idx " << before_seg_Idx << std::endl;
    for (int j = 1; j < minPairSecond.size(); j++) {

        if (before_seg_Idx != minPairSecond[j]->seg_Idx) {

            Eigen::Vector3d p1, p2;
            minPairSecond[j]->GetCoord3D(p1[0], p1[1], p1[2]);
            minPairSecond[j - 1]->GetCoord3D(p2[0], p2[1], p2[2]);
            double _gap_dist = (p1 - p2).norm();

            this->buildNewEdgetoQMeshPatch(_toolpath,
                minPairSecond[j],
                minPairSecond[j - 1]);

            if (_gap_dist < gap_dist_threshold) { //only smooth shortLinkage
                minPairSecond[j]->is_endPnt = false;
                minPairSecond[j - 1]->is_endPnt = false;
            }

        }

        before_seg_Idx = minPairSecond[j]->seg_Idx;
    }

    /*Close chain part*/
    //get the first pnt of close Chain
    QMeshNode* endNode_closeLoop = minPairSecond.back();

    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);
        Node->visited = false; //clean flag
    }
    for (GLKPOSITION Pos = _toolpath->GetEdgeList().GetHeadPosition(); Pos;) {
        QMeshEdge* Edge = (QMeshEdge*)_toolpath->GetEdgeList().GetNext(Pos);
        Edge->visited = false; //clean the visited Flag for edges
    }

    do {
        QMeshNode* startNode_closeLoop = this->_find_nearest_closeLoopNode(_toolpath, endNode_closeLoop);


        if (startNode_closeLoop == NULL) {
            // all of the unshortChain nodes are handled (open chain OR transfered open chain from losed chaine)
            endNode_closeLoop->is_endPnt = true;
            break;
        }


        startNode_closeLoop->visited = true;
        QMeshEdge* transfer_edge = this->buildNewEdgetoQMeshPatch(_toolpath, endNode_closeLoop, startNode_closeLoop);
        transfer_edge->visited = true;

        //handle the is_endPnt flag
        if (transfer_edge->CalLength() < gap_dist_threshold) { //only smooth shortLinkage
            endNode_closeLoop->is_endPnt = false;
            startNode_closeLoop->is_endPnt = false;
        }
        else {
            endNode_closeLoop->is_endPnt = true;
            startNode_closeLoop->is_endPnt = true;
        }

        endNode_closeLoop = this->_open_closed_loop(_toolpath, startNode_closeLoop);

    } while (!this->_all_closeLoop_opened(_toolpath));

    endNode_closeLoop->is_endPnt = true; // the last node of opened closeLoop should be endPnt

    return minPairSecond[0];
}

QMeshNode* toolpathGeneration_stripe::_get_another_endPnt_sameChain(
    QMeshNode* startNode, std::vector<QMeshNode*> endPnts) {

    int startNode_seg_Idx = startNode->seg_Idx;

    QMeshNode* another_endPnt_sameChain = NULL;

    for (size_t i = 0; i < endPnts.size(); ++i) {

        if (endPnts[i]->visited) continue; //skip self and other visited pnts

        if (endPnts[i]->seg_Idx == startNode_seg_Idx) {
            another_endPnt_sameChain = endPnts[i];
            break;
        }

    }

    if (another_endPnt_sameChain == NULL)
        std::cout << "Error: @_get_another_endPnt_sameChain" << std::endl;

    return another_endPnt_sameChain;
}

QMeshNode* toolpathGeneration_stripe::_find_next_nearest_startNode(
    QMeshNode* another_endPnt_sameChain, std::vector<QMeshNode*> endPnts) {

    Eigen::Vector3d p1;
    another_endPnt_sameChain->GetCoord3D(p1[0], p1[1], p1[2]);

    QMeshNode* next_nearest_startNode = NULL;
    double min_dist = 1e10;

    for (size_t i = 0; i < endPnts.size(); ++i) {

        if (endPnts[i]->visited) continue; //skip self and other visited pnts

        Eigen::Vector3d p2;
        endPnts[i]->GetCoord3D(p2[0], p2[1], p2[2]);
        double temp_dist = (p1 - p2).norm();

        if (temp_dist < min_dist) {
            min_dist = temp_dist;
            next_nearest_startNode = endPnts[i];
        }
    }

    if (next_nearest_startNode == NULL)
        std::cout << "Error: @_find_next_nearest_startNode" << std::endl;

    return next_nearest_startNode;

}

QMeshNode* toolpathGeneration_stripe::_find_nearest_or_second_nearest_startNode(
    QMeshNode* another_endPnt_sameChain, std::vector<QMeshNode*> endPnts, bool returnSecondNearest) {

    // Early check for empty input
    if (endPnts.empty()) {
        std::cout << "Error: @_find_nearest_or_second_nearest_startNode - endPnts is empty" << std::endl;
        return nullptr;
    }

    Eigen::Vector3d p1;
    another_endPnt_sameChain->GetCoord3D(p1[0], p1[1], p1[2]);

    QMeshNode* nearest_startNode = nullptr;
    QMeshNode* second_nearest_startNode = nullptr;
    double min_dist = std::numeric_limits<double>::max();
    double second_min_dist = std::numeric_limits<double>::max();

    for (auto& node : endPnts) {
        if (node->visited) continue; // skip visited nodes

        Eigen::Vector3d p2;
        node->GetCoord3D(p2[0], p2[1], p2[2]);
        double temp_dist = (p1 - p2).norm();

        if (temp_dist < min_dist) {
            // Update the second nearest node to the previous nearest node
            second_min_dist = min_dist;
            second_nearest_startNode = nearest_startNode;

            // Update the nearest node
            min_dist = temp_dist;
            nearest_startNode = node;
        }
        else if (temp_dist < second_min_dist) {
            // Update the second nearest node
            second_min_dist = temp_dist;
            second_nearest_startNode = node;
        }
    }

    if (nearest_startNode == nullptr) {
        std::cout << "Error: @_find_nearest_or_second_nearest_startNode - No valid start node found" << std::endl;
    }

    // Return either the nearest or the second nearest node based on the boolean flag
    return returnSecondNearest ? second_nearest_startNode : nearest_startNode;
}


bool toolpathGeneration_stripe::_all_visited(std::vector<QMeshNode*>endPnts) {

    bool _visited_all = true;

    for (size_t i = 0; i < endPnts.size(); ++i) {

        if (!endPnts[i]->visited)
            _visited_all = false;
    }

    return _visited_all;
}

void toolpathGeneration_stripe::_trace_2_oneStroke_path(QMeshPatch* _toolpath, QMeshNode* sNode) {

    //pick the Pnts in the _toolpath (open curve)
    std::vector<QMeshNode*> toolpath_NodeSet;

    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);
        Node->connectTPathProcessed = false; //initialize
        if (Node->is_shortChain || Node->is_closeLoop)
            Node->connectTPathProcessed = true;
    }

    QMeshNode* startNode = sNode;
    toolpath_NodeSet.push_back(startNode);
    startNode->connectTPathProcessed = true;

    QMeshNode* nextNode = NULL;
    do {

        for (GLKPOSITION Pos = startNode->GetEdgeList().GetHeadPosition(); Pos;) {
            QMeshEdge* neighEdge = (QMeshEdge*)startNode->GetEdgeList().GetNext(Pos);

            nextNode = neighEdge->GetStartPoint();
            if (nextNode == startNode) nextNode = neighEdge->GetEndPoint();

            if (nextNode->connectTPathProcessed) continue;
            toolpath_NodeSet.push_back(nextNode);
            nextNode->connectTPathProcessed = true;
        }

        startNode = toolpath_NodeSet.back();

        //std::cout << startNode->GetIndexNo() << std::endl;

    } while (!this->_detectAll_openChain_Processed(_toolpath));


    //rebuild the node list
    _toolpath->GetNodeList().RemoveAll();
    for (int i = 0; i < toolpath_NodeSet.size(); i++) {

        toolpath_NodeSet[i]->SetIndexNo(_toolpath->GetNodeList().GetCount()); //index should start from 0
        _toolpath->GetNodeList().AddTail(toolpath_NodeSet[i]);

        toolpath_NodeSet[i]->GetEdgeList().RemoveAll();
    }

    //rebuild the edge list
    _toolpath->GetEdgeList().RemoveAll();
    for (int i = 0; i < (toolpath_NodeSet.size() - 1); i++) {
        QMeshEdge* _edge = this->buildNewEdgetoQMeshPatch(_toolpath, toolpath_NodeSet[i], toolpath_NodeSet[i + 1]);
        _edge->seg_Idx = 0;
    }

}

bool toolpathGeneration_stripe::_detectAll_openChain_Processed(QMeshPatch* singlePath) {
    //if all the node being processed return true, else return false.
    for (GLKPOSITION Pos = singlePath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* thisNode = (QMeshNode*)singlePath->GetNodeList().GetNext(Pos);

        if (thisNode->is_shortChain || thisNode->is_closeLoop) continue;

        if (thisNode->connectTPathProcessed == false) return false;
    }
    return true;
}

QMeshNode* toolpathGeneration_stripe::_find_nearest_closeLoopNode(QMeshPatch* _toolpath, QMeshNode* endNode_closeLoop) {

    QMeshNode* nearest_startNode_closeLoop = NULL;
    double min_dist = 1e10;
    Eigen::Vector3d p1; endNode_closeLoop->GetCoord3D(p1[0], p1[1], p1[2]);

    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);

        if (!Node->is_closeLoop || Node->is_shortChain || Node->visited) continue;

        Eigen::Vector3d p2; Node->GetCoord3D(p2[0], p2[1], p2[2]);
        double _dist = (p1 - p2).norm();
        if (_dist < min_dist) {
            min_dist = _dist;
            nearest_startNode_closeLoop = Node;
        }

    }

    if (nearest_startNode_closeLoop == NULL) {
        std::cout << "Error: @ _find_nearest_closeLoopNode, may becase there is no close loop, ignore " << std::endl;
    }

    return nearest_startNode_closeLoop;
}

QMeshNode* toolpathGeneration_stripe::_open_closed_loop(QMeshPatch* _toolpath, QMeshNode* startNode_closeLoop) {

    //find the neighbor node
    QMeshNode* picked_neighNode = NULL;
    QMeshEdge* picked_neighEdge = NULL;

    for (GLKPOSITION Pos = startNode_closeLoop->GetEdgeList().GetHeadPosition(); Pos;) {
        QMeshEdge* neighEdge = (QMeshEdge*)startNode_closeLoop->GetEdgeList().GetNext(Pos);

        if (neighEdge->visited) continue;

        picked_neighNode = neighEdge->GetStartPoint();
        if (picked_neighNode == startNode_closeLoop) picked_neighNode = neighEdge->GetEndPoint();

        picked_neighEdge = neighEdge;
        if (picked_neighNode != NULL) break;
    }

    if (picked_neighNode == NULL) {
        std::cout << "Error: @ _open_closed_loop " << std::endl;
    }

    //delete the picked_neighEdge
    _toolpath->GetEdgeList().Remove(picked_neighEdge);
    startNode_closeLoop->GetEdgeList().Remove(picked_neighEdge);
    picked_neighNode->GetEdgeList().Remove(picked_neighNode);

    //mark the visited flag for one_close_loop
    int seg_idx = startNode_closeLoop->seg_Idx;
    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);

        if (Node->seg_Idx == seg_idx) {
            Node->visited = true;
            Node->is_closeLoop = false;
        }
    }

    return picked_neighNode;
}

bool toolpathGeneration_stripe::_all_closeLoop_opened(QMeshPatch* _toolpath) {

    for (GLKPOSITION Pos = _toolpath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* thisNode = (QMeshNode*)_toolpath->GetNodeList().GetNext(Pos);

        if (thisNode->is_shortChain || !thisNode->is_closeLoop) continue;

        if (thisNode->visited == false) return false;
    }
    return true;
}

QMeshEdge* toolpathGeneration_stripe::buildNewEdgetoQMeshPatch(QMeshPatch* patch, QMeshNode* startNode, QMeshNode* endNode) {

    QMeshEdge* isoEdge = new QMeshEdge;

    //std::cout << "start Node " << startNode->GetIndexNo() << " end Node " << endNode->GetIndexNo() << std::endl;

    isoEdge->SetStartPoint(startNode);
    isoEdge->SetEndPoint(endNode);

    isoEdge->SetMeshPatchPtr(patch);
    //isoEdge->SetIndexNo(patch->GetEdgeList().GetCount() + 1);
    isoEdge->SetIndexNo(patch->GetEdgeList().GetCount());

    (startNode->GetEdgeList()).AddTail(isoEdge);
    (endNode->GetEdgeList()).AddTail(isoEdge);
    patch->GetEdgeList().AddTail(isoEdge);

    return isoEdge;
}

//dist and jump flag
void toolpathGeneration_stripe::_cal_jumpFlag_4_CCF_toolpath(QMeshPatch* ccfPath) {

    for (GLKPOSITION nodePos = ccfPath->GetNodeList().GetHeadPosition(); nodePos;) {
        QMeshNode* Node = (QMeshNode*)ccfPath->GetNodeList().GetNext(nodePos);

        Eigen::Vector3d pp, pp_nxt;
        Node->GetCoord3D(pp[0], pp[1], pp[2]);

        double D = 0.0;
        int lines = ccfPath->GetNodeNumber();
        if (Node->GetIndexNo() == (lines - 1)) { D = 1.0; }
        else {

            GLKPOSITION nextPos = ccfPath->GetNodeList().Find(Node)->next;
            QMeshNode* nextNode = (QMeshNode*)ccfPath->GetNodeList().GetAt(nextPos);
            nextNode->GetCoord3D(pp_nxt[0], pp_nxt[1], pp_nxt[2]);

            D = (pp - pp_nxt).norm();
            if (D > m_jump_detection_threshold) {
                D = 0.0;								// inject the D to the Node/startPnt of Edge
                Node->Jump_preSecEnd = true;			// end of prev section
                nextNode->Jump_nextSecStart = true;		// start of next section
            }
        }
    }
}

//cut info of nodes
void toolpathGeneration_stripe::_cal_cutInfo_4_CCF_toolpath(QMeshPatch* ccfPath) {

    //soft-segment by jump and large turning angle flags
    this->_do_soft_segment_Jump(ccfPath);

    //build table for each segment
    std::vector<std::vector<QMeshNode*>> seg_Set;
    std::vector<QMeshNode*> onePath;
    int global_idx = 0;
    for (GLKPOSITION Pos = ccfPath->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)ccfPath->GetNodeList().GetNext(Pos);

        if (Node->seg_Idx != global_idx) {

            std::vector<QMeshNode*> copy_onePath;
            copy_onePath.assign(onePath.begin(), onePath.end());
            seg_Set.push_back(copy_onePath);
            global_idx = Node->seg_Idx;
            onePath.clear();
        }

        onePath.push_back(Node);
    }

    // Check if onePath still contains nodes and add it to seg_Set
    if (!onePath.empty()) {
        std::vector<QMeshNode*> copy_onePath;
        copy_onePath.assign(onePath.begin(), onePath.end());
        seg_Set.push_back(copy_onePath);
    }


    /*for (int i = 0; i < seg_Set.size(); i++) {
        std::cout << i << " seg_Set[i].size() = " << seg_Set[i].size() << std::endl;
    }*/

    //mark the ccf_cut flag
    for (int i = 0; i < seg_Set.size(); i++) {

        //get the length of each seg
        double _length = 0;
        for (int j = 0; j < (seg_Set[i].size() - 1); j++) {

            Eigen::Vector3d _a, _b;
            seg_Set[i][j]->GetCoord3D(_a[0], _a[1], _a[2]);
            seg_Set[i][j + 1]->GetCoord3D(_b[0], _b[1], _b[2]);
            _length += (_a - _b).norm();

        }
        // discard the short path
        if (_length < min_ccf_length) {

            for (int j = 0; j < seg_Set[i].size(); j++) {
                seg_Set[i][j]->seg_Idx = -1;
            }

            continue;
        }
        //std::cout << i << " seg_Set[i] length = " << _length << std::endl;

        //trace back to min_ccf_length and mark cut info
        _length = 0;
        for (int k = seg_Set[i].size() - 1; k >= 1; --k) {

            Eigen::Vector3d _a, _b;
            seg_Set[i][k]->GetCoord3D(_a[0], _a[1], _a[2]);
            seg_Set[i][k - 1]->GetCoord3D(_b[0], _b[1], _b[2]);
            _length += (_a - _b).norm();

            if (_length > min_ccf_length) {
                seg_Set[i][k - 1]->cut_info = true;
                break;
            }
        }
    }
}

void toolpathGeneration_stripe::_do_soft_segment_Jump(QMeshPatch* patch) {

    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);
        Node->seg_Idx = -1;
    }

    //show diff section by color (each section)
    int global_seg_id = 0;

    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        Node->seg_Idx = global_seg_id;

        if (Node->Jump_preSecEnd)
            global_seg_id++;

        //std::cout << "seg_Idx " << Node->seg_Idx << std::endl;
    }
    std::cout << "global_seg_id " << global_seg_id << std::endl;
}

//tangent dir of nodes
void toolpathGeneration_stripe::_cal_tangentDir_4_CCF_toolpath(QMeshPatch* ccfPath) {

    Eigen::Vector3d last_tangent_dir = Eigen::Vector3d::Zero();
    for (GLKPOSITION nodePos = ccfPath->GetNodeList().GetHeadPosition(); nodePos;) {
        QMeshNode* Node = (QMeshNode*)ccfPath->GetNodeList().GetNext(nodePos);

        Eigen::Vector3d pp, pp_n;
        Node->GetCoord3D(pp[0], pp[1], pp[2]);

        double D = 0.0;
        int lines = ccfPath->GetNodeNumber();
        if (Node->GetIndexNo() == (lines - 1)) { Node->m_printTan = last_tangent_dir; }
        else {

            GLKPOSITION nextPos = ccfPath->GetNodeList().Find(Node)->next;
            QMeshNode* nextNode = (QMeshNode*)ccfPath->GetNodeList().GetAt(nextPos);
            nextNode->GetCoord3D(pp_n[0], pp_n[1], pp_n[2]);

            D = (pp_n - pp).norm();

            if (D > 0.01 && !Node->Jump_preSecEnd) {
                Eigen::Vector3d tangent_dir = (pp_n - pp).normalized();

                last_tangent_dir = tangent_dir;
            }
            else {
                std::cout << "The distance between current node to the next node is Zero" << std::endl;
            }
        }
        Node->m_printTan = last_tangent_dir;
    }
}