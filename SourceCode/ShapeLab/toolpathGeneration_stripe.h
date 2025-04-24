#pragma once
#include "../QMeshLib/PolygenMesh.h"

class toolpathGeneration_stripe {

public:
	toolpathGeneration_stripe(PolygenMesh* isoLayerSet, PolygenMesh* toolpathSet) {
		m_Slices = isoLayerSet;
		m_Waypoints = toolpathSet;
	};
	~toolpathGeneration_stripe() {};

	void generate_toolPath_marchingSquare();
	void postProcessing_stripPath();

private:

	void _initial_index(QMeshPatch* patch);
	size_t _getOrderedEdgeIndexOfFace(QMeshFace* Face, size_t _id_edge);
	std::vector<double> _calculateIsoPoints(double s0, double s1);
	std::vector<Eigen::Vector2i> _matchCrossings(const std::vector<std::vector<size_t>>& crossings);
	void _buildPath(const std::vector<Eigen::Vector3d>& vertices,
		const std::vector<Eigen::Vector2i>& edges, QMeshPatch* singlePath);
	void _mark_boundary_face(QMeshPatch* layer);
	void _build_connection_between_isoNode_with_relatedEdge(
		QMeshPatch* singlePath, QMeshPatch* layer, const std::vector<std::vector<size_t>>& polylineIndices);
	void _build_edges_in_singular_region(QMeshPatch* singlePath, QMeshPatch* layer);
	void _record_edgeOfPath_on_which_faceOFLayer(QMeshPatch* singlePath);

	void _classify_chains(QMeshPatch* _toolpath);
	//open loop
	QMeshNode* _pickUp_one_startNode(QMeshPatch* _toolpath, int seg_Idx);
	void _tracing_one_path(QMeshPatch* _toolpath, QMeshNode* start_Node, int seg_Idx);
	bool _allDetced_4_openChains(QMeshPatch* _toolpath);
	//close loop
	QMeshNode* _pickUp_one_startNode_closeLoop(QMeshPatch* _toolpath, int seg_Idx);
	void _tracing_one_path_closeLoop(QMeshPatch* _toolpath, QMeshNode* start_Node, int seg_Idx);
	bool _allDetced_4_openChains_closeLoop(QMeshPatch* _toolpath);
	void _mark_short_chains(QMeshPatch* _toolpath, double cutLength);

	void _smooth_chains(QMeshPatch* _toolpath, int loop, double ratio);
	QMeshNode* _connect_shortGap(QMeshPatch* _toolpath, double gap_dist_threshold);
	QMeshNode* _get_another_endPnt_sameChain(QMeshNode* startNode, std::vector<QMeshNode*> endPnts);
	QMeshNode* _find_next_nearest_startNode(QMeshNode* another_endPnt_sameChain, std::vector<QMeshNode*> endPnts);
	//updated version
	QMeshNode* _find_nearest_or_second_nearest_startNode(
		QMeshNode* another_endPnt_sameChain, std::vector<QMeshNode*> endPnts, bool returnSecondNearest = false);

	bool _all_visited(std::vector<QMeshNode*> endPnts);
	QMeshNode* _find_nearest_closeLoopNode(QMeshPatch* _toolpath, QMeshNode* endNode_closeLoop);
	QMeshNode* _open_closed_loop(QMeshPatch* _toolpath, QMeshNode* startNode_closeLoop);
	bool _all_closeLoop_opened(QMeshPatch* _toolpath);

	void _trace_2_oneStroke_path(QMeshPatch* _toolpath, QMeshNode* sNode);
	bool _detectAll_openChain_Processed(QMeshPatch* singlePath);

	QMeshEdge* buildNewEdgetoQMeshPatch(QMeshPatch* patch, QMeshNode* startNode, QMeshNode* endNode);

	void _cal_jumpFlag_4_CCF_toolpath(QMeshPatch* ccfPath);
	void _cal_cutInfo_4_CCF_toolpath(QMeshPatch* ccfPath);
	void _do_soft_segment_Jump(QMeshPatch* patch);

	void _cal_tangentDir_4_CCF_toolpath(QMeshPatch* ccfPath);

private:
	PolygenMesh* m_Slices = NULL;
	PolygenMesh* m_Waypoints = NULL;
	double m_jump_detection_threshold = 6.0;
	double min_ccf_length = 100.0;

};