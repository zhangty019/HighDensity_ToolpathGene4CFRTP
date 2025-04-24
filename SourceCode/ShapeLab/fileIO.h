#pragma once

#include "../QMeshLib/PolygenMesh.h"
#include <unordered_map>
#include <vector>
#include <string>

class fileIO
{
public:
    fileIO() {}
    ~fileIO() {}

    int read_layer_files(PolygenMesh* Slices, std::string Dir);
    void input_Projeted_stressField(PolygenMesh* isoLayerSet, std::string dir);
    void outputDirectionField_onNodes(PolygenMesh* isoLayerSet, std::string dir);
    int get_layer_files_name_set(
        std::string Dir, std::vector<std::string>& layer_name_set);
    int get_dField_files_name_set(
        std::string Dir, std::vector<std::string>& dField_name_set);
    bool _find_matched_dField(const std::string& layer_name,
        const std::vector<std::string>& dField_name_set);
    bool read_dField_onNode(const std::string& file_path,
        std::vector<std::vector<double>>& matrix);
    void outputIsoSurface(PolygenMesh* isoSurface, std::string path);
    void output_rawStripToolpath(PolygenMesh* toolpathSet, std::string path);
    void input_raw_stripPath(PolygenMesh* isoLayerSet, PolygenMesh* toolpathSet, std::string path);
    void output_StripToolpath(PolygenMesh* toolpathSet, std::string path);

private:
    std::vector<std::string> _getFilesInDirectory(
        const std::string& dir, const std::string& extension);
    void _readSliceData(PolygenMesh* Slices, std::string folderName,
        std::vector<std::string> layer_name_set);
    void _sortFilesByLeadingNumber(std::vector<std::string>& fileNames);
    bool _matchFileName(const std::vector<std::string>& fileNames,
        const std::string& inputString);
    std::unordered_map<int, std::tuple<double, double, double>> _parseFile(
        const std::string& filePath);
    int _remove_allFile_in_Dir(const std::string& dirPath);
    void _output_OneSurfaceMesh(QMeshPatch* eachLayer, std::string path);
    void _output_One_rawStripToolpath(QMeshPatch* eachPath, std::string path);

    void _read_one_pathObjFile(PolygenMesh* toolpathSet, std::string path, std::string name, int givenId);
    void _assign_waypointNormal_from_layerFacet(
        PolygenMesh* isoLayerSet, PolygenMesh* toolpathSet, std::string name);
    Eigen::Vector3d _get_faceNormal_with_idx(const std::vector<QMeshFace*>& faceSet_of_layer, int idx);

    void _output_One_stripToolpath(QMeshPatch* eachPath, std::string path);
};
