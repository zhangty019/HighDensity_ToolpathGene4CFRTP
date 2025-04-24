#include "fileIO.h"
#include <fstream>
#include <io.h>
#include <iostream>
#include <algorithm>
#include <regex>
#include <unordered_map>
#include <sstream>
#include "dirent.h"

int fileIO::read_layer_files(PolygenMesh* Slices, std::string Dir) {
    std::vector<std::string> layer_name_set = this->_getFilesInDirectory(Dir, ".obj");
    this->_sortFilesByLeadingNumber(layer_name_set);
    this->_readSliceData(Slices, Dir, layer_name_set);
    return Slices->GetMeshList().GetCount();
}

void fileIO::_readSliceData(PolygenMesh* Slices, std::string folderName, std::vector<std::string> layer_name_set) {
    //read slice files and build mesh_patches
    char filename[1024];
    for (size_t i = 0; i < layer_name_set.size(); i++) {
        sprintf_s(filename, sizeof(filename), "%s/%s", folderName.c_str(), layer_name_set[i].c_str());

        QMeshPatch* layers = new QMeshPatch;
        layers->SetIndexNo(Slices->GetMeshList().GetCount()); //index begin from 0
        Slices->GetMeshList().AddTail(layers);
        layers->patchName = layer_name_set[i].c_str();

        layers->inputOBJFile(filename);
    }
    std::cout << "--> Finish inputting files (curved_layers) from: " << folderName << std::endl;
}

std::vector<std::string> fileIO::_getFilesInDirectory(const std::string& dir, const std::string& extension) {
    std::vector<std::string> fileNames;
    DIR* dp;
    struct dirent* entry;

    if ((dp = opendir(dir.c_str())) == NULL) {
        std::cerr << "Cannot open directory: " << dir << std::endl;
        return fileNames;
    }

    while ((entry = readdir(dp)) != NULL) {
        std::string filename = entry->d_name;
        if (filename.find(extension, (filename.length() - extension.length())) != std::string::npos) {
            fileNames.push_back(filename);
        }
    }

    closedir(dp);
    return fileNames;
}

void fileIO::_sortFilesByLeadingNumber(std::vector<std::string>& fileNames) {
    std::sort(fileNames.begin(), fileNames.end(), [](const std::string& a, const std::string& b) {
        std::regex re(R"((\d+))");
        std::smatch matchA, matchB;

        if (std::regex_search(a, matchA, re) && std::regex_search(b, matchB, re)) {
            return std::stoi(matchA.str(1)) < std::stoi(matchB.str(1));
        }
        return a < b; // Fallback to alphabetical order if no leading number is found
        });
}

void fileIO::input_Projeted_stressField(PolygenMesh* isoLayerSet, std::string dir) {
    //get the stress file name set
    std::vector<std::string> files = this->_getFilesInDirectory(dir, ".txt");

    for (GLKPOSITION posMesh = isoLayerSet->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* each_layer = (QMeshPatch*)isoLayerSet->GetMeshList().GetNext(posMesh);

        if (this->_matchFileName(files, each_layer->patchName)) {
            each_layer->isStressLayer = true;
        }
        else { each_layer->isStressLayer = false; }

        if (!each_layer->isStressLayer) continue;

        std::string filePath = dir + "/" + each_layer->patchName + "_stress.txt";

        auto dataMap = this->_parseFile(filePath);

        for (GLKPOSITION posFace = each_layer->GetFaceList().GetHeadPosition(); posFace != nullptr;) {
            QMeshFace* face = (QMeshFace*)each_layer->GetFaceList().GetNext(posFace);

            int inputKey = face->GetIndexNo();

            if (dataMap.find(inputKey) != dataMap.end()) {
                auto tuple = dataMap[inputKey];
                double x = std::get<0>(tuple);
                double y = std::get<1>(tuple);
                double z = std::get<2>(tuple);
                face->stessField_vector << x, y, z;
                face->stessField_flag = true;

                /* Code Continuity*/
                Eigen::Vector3d inVec = face->stessField_vector;
                double stressMag = inVec.stableNorm();
                inVec.stableNormalize();
                face->stressMag = stressMag;
                face->isConstraint = true;
                face->vectorDir = inVec;
                face->stressField = inVec;
            }
            else {
                face->stessField_flag = false;
            }
        }
    }
    std::cout << "--> Finish inputting files (projected_stressField) from: " << dir << std::endl;
}

bool fileIO::_matchFileName(const std::vector<std::string>& fileNames, const std::string& inputString) {
    std::regex pattern(R"((.*\.obj)_stress\.txt)");

    for (const auto& fileName : fileNames) {
        std::smatch match;
        if (std::regex_match(fileName, match, pattern) && match.size() > 1) {
            if (match[1].str() == inputString) {
                return true;
            }
        }
    }

    return false;
}

std::unordered_map<int, std::tuple<double, double, double>> fileIO::_parseFile(const std::string& filePath) {
    std::unordered_map<int, std::tuple<double, double, double>> dataMap;
    std::ifstream file(filePath);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return dataMap;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string identifier;
        double x, y, z;

        if (iss >> identifier >> x >> y >> z) {
            if (identifier[0] == 'f') {
                int key = std::stoi(identifier.substr(1));
                dataMap[key] = std::make_tuple(x, y, z);
            }
        }
    }

    file.close();
    return dataMap;
}

void fileIO::outputDirectionField_onNodes(PolygenMesh* isoLayerSet, std::string dir) {
    // Clean the directory before writing new files
    _remove_allFile_in_Dir(dir);

    for (GLKPOSITION PatchPos = isoLayerSet->GetMeshList().GetHeadPosition(); PatchPos;) {
        QMeshPatch* Patch = (QMeshPatch*)isoLayerSet->GetMeshList().GetNext(PatchPos);

        if (!Patch->isStressLayer) continue;

        std::string patchName = Patch->patchName;
        size_t last_dot = patchName.find_last_of('.');
        if (last_dot != std::string::npos) {
            patchName = patchName.substr(0, last_dot) + ".txt";
        }
        std::string filename = dir + "/" + patchName;

        std::ofstream oFile(filename, std::ofstream::out);
        if (!oFile.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            continue;
        }

        int index = 0;
        for (GLKPOSITION nodePos = Patch->GetNodeList().GetHeadPosition(); nodePos;) {
            QMeshNode* node = (QMeshNode*)Patch->GetNodeList().GetNext(nodePos);

            oFile << index
                << " " << node->directionField[0]
                << " " << node->directionField[1]
                << " " << node->directionField[2]
                << "\n";
            index++;
        }
        oFile.close();
    }

    //output the large area Layers index
    double areaThreshold = 3e3;
    std::cout << "--> Output the layer index with large area: ( " << areaThreshold << " )" << std::endl;
    for (GLKPOSITION PatchPos = isoLayerSet->GetMeshList().GetHeadPosition(); PatchPos;) {
        QMeshPatch* Patch = (QMeshPatch*)isoLayerSet->GetMeshList().GetNext(PatchPos);

        double patchArea = 0;
        for (GLKPOSITION facePos = Patch->GetFaceList().GetHeadPosition(); facePos;) {
            QMeshFace* face = (QMeshFace*)Patch->GetFaceList().GetNext(facePos);

            patchArea += face->CalArea();
        }
        if(patchArea > areaThreshold)
            std::cout << Patch->GetIndexNo() << " ,";
    }
    std::cout << std::endl;

    std::cout << "--> Finish outputting files (directionField_on_Node) into: " << dir << std::endl;
}

int fileIO::_remove_allFile_in_Dir(const std::string& dirPath) {
    DIR* dp;
    struct dirent* entry;

    if ((dp = opendir(dirPath.c_str())) == NULL) {
        std::cerr << "Cannot open directory: " << dirPath << std::endl;
        return -1;
    }

    while ((entry = readdir(dp)) != NULL) {
        std::string filename = entry->d_name;
        // skip "." and ".."
        if (filename == "." || filename == "..") {
            continue;
        }
        std::string filepath = dirPath + "/" + filename;
        if (remove(filepath.c_str()) != 0) {
            std::cerr << "Error deleting file: " << filepath << std::endl;
        }
    }

    closedir(dp);
    return 0;
}

int fileIO::get_layer_files_name_set(std::string Dir, std::vector<std::string>& layer_name_set) {
    layer_name_set = this->_getFilesInDirectory(Dir, ".obj");
    this->_sortFilesByLeadingNumber(layer_name_set);
    return layer_name_set.size();
}

int fileIO::get_dField_files_name_set(std::string Dir, std::vector<std::string>& dField_name_set) {
    dField_name_set = this->_getFilesInDirectory(Dir, ".txt");
    this->_sortFilesByLeadingNumber(dField_name_set);
    return dField_name_set.size();
}

bool fileIO::_find_matched_dField(const std::string& layer_name, const std::vector<std::string>& dField_name_set) {
    std::string dField_name = layer_name.substr(0, layer_name.find_last_of('.')) + ".txt";
    return std::find(dField_name_set.begin(), dField_name_set.end(), dField_name) != dField_name_set.end();
}

bool fileIO::read_dField_onNode(const std::string& file_path, std::vector<std::vector<double>>& matrix) {
    std::ifstream file(file_path);

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::vector<double> row;
            std::istringstream iss(line);
            double num;
            while (iss >> num) {
                row.push_back(num);
            }
            matrix.push_back(row);
        }
        file.close();

        for (auto& row : matrix) {
            if (row.size() > 1) {
                row.erase(row.begin());
            }
        }
        return true;
    }
    else {
        std::cout << "Unable to open file: " << file_path << std::endl;
        return false;
    }
}

void fileIO::outputIsoSurface(PolygenMesh* isoSurface, std::string path) {

    this->_remove_allFile_in_Dir(path);

    std::string LAYER_dir;

    for (GLKPOSITION posMesh = isoSurface->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* each_layer = (QMeshPatch*)isoSurface->GetMeshList().GetNext(posMesh);

        LAYER_dir = path + "/" + each_layer->patchName;

        this->_output_OneSurfaceMesh(each_layer, LAYER_dir);
    }
}

void fileIO::_output_OneSurfaceMesh(QMeshPatch* eachLayer, std::string path) {

    double pp[3];
    std::ofstream nodeSelection(path);

    int index = 0;
    for (GLKPOSITION posNode = eachLayer->GetNodeList().GetHeadPosition(); posNode != nullptr;) {
        QMeshNode* node = (QMeshNode*)eachLayer->GetNodeList().GetNext(posNode);
        node->GetCoord3D(pp[0], pp[1], pp[2]);
        nodeSelection << "v " << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;
        index++; node->SetIndexNo(index);

    }
    for (GLKPOSITION posFace = eachLayer->GetFaceList().GetHeadPosition(); posFace != nullptr;) {
        QMeshFace* face = (QMeshFace*)eachLayer->GetFaceList().GetNext(posFace);
        nodeSelection << "f " << face->GetNodeRecordPtr(0)->GetIndexNo()
            << " " << face->GetNodeRecordPtr(1)->GetIndexNo()
            << " " << face->GetNodeRecordPtr(2)->GetIndexNo() << std::endl;
    }
    nodeSelection.close();
}

void fileIO::output_rawStripToolpath(PolygenMesh* toolpathSet, std::string path) {

    this->_remove_allFile_in_Dir(path);

    std::string Path_dir;

    for (GLKPOSITION posMesh = toolpathSet->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* each_path = (QMeshPatch*)toolpathSet->GetMeshList().GetNext(posMesh);

        Path_dir = path + "/" + each_path->patchName;

        this->_output_One_rawStripToolpath(each_path, Path_dir);
    }
}

void fileIO::_output_One_rawStripToolpath(QMeshPatch* eachPath, std::string path) {

    double pp[3];
    std::ofstream file_(path);

    int index = 0;
    for (GLKPOSITION posNode = eachPath->GetNodeList().GetHeadPosition(); posNode != nullptr;) {
        QMeshNode* node = (QMeshNode*)eachPath->GetNodeList().GetNext(posNode);
        node->GetCoord3D(pp[0], pp[1], pp[2]);
        file_ << "v " << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;
        index++; node->SetIndexNo(index);

    }
    for (GLKPOSITION posEdge = eachPath->GetEdgeList().GetHeadPosition(); posEdge != nullptr;) {
        QMeshEdge* edge = (QMeshEdge*)eachPath->GetEdgeList().GetNext(posEdge);
        file_ << "l " << edge->GetStartPoint()->GetIndexNo()
            << " " << edge->GetEndPoint()->GetIndexNo() << std::endl;
    }

    for (GLKPOSITION posEdge = eachPath->GetEdgeList().GetHeadPosition(); posEdge != nullptr;) {
        QMeshEdge* edge = (QMeshEdge*)eachPath->GetEdgeList().GetNext(posEdge);
        //+1 is compatible with the below function: _get_faceNormal_with_idx 
        file_ << "#li " << (edge->relatedLayerFace->GetIndexNo() + 1) << std::endl;
    }

    file_.close();
}


void fileIO::input_raw_stripPath(PolygenMesh* isoLayerSet, PolygenMesh* toolpathSet, std::string path) {

    // collect the names of the isoLayerset
    std::vector<std::string> layer_list;
    for (GLKPOSITION posMesh = isoLayerSet->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* each_layer = (QMeshPatch*)isoLayerSet->GetMeshList().GetNext(posMesh);

        //layer index is already start from 0;
        layer_list.push_back(each_layer->patchName);

        //calculate the faceNormal once to speed up the Function: _get_faceNormal_with_idx
        for (GLKPOSITION posFace = each_layer->GetFaceList().GetHeadPosition(); posFace != nullptr;) {
            QMeshFace* face = (QMeshFace*)each_layer->GetFaceList().GetNext(posFace);

            face->CalPlaneEquation();
        }
    }

    // output name list
    for (const auto& layer : layer_list) {
        std::cout << "layer: " << layer << std::endl;
    }

    // collect the name of the path folder
    std::vector<std::string> stripPath_list;
    DIR* dir;
    struct dirent* ent;
    if ((dir = opendir(path.c_str())) != nullptr) {
        while ((ent = readdir(dir)) != nullptr) {
            std::string filename(ent->d_name);
            // Check if the file is a .obj file
            if (filename.size() > 4 && filename.substr(filename.size() - 4) == ".obj") {
                // Add the filename (with extension) to the path_list
                stripPath_list.push_back(filename);
            }
        }
        closedir(dir);
    }
    else {
        // Could not open directory
        perror("Could not open directory");
        return;
    }

    // output name list
    for (const auto& strip_path : stripPath_list) {
        std::cout << "strip_path: " << strip_path << std::endl;
    }

    // find the same name of the stripPath file by iterating through each name in layer_list
    bool nameFound = false;
    for (size_t i = 0; i < layer_list.size(); ++i) {
        const std::string& name = layer_list[i];
        // Check if the name exists in stripPath_list
        auto it = std::find(stripPath_list.begin(), stripPath_list.end(), name);
        // If name is found in stripPath_list, set the flag to true and break out of the loop
        if (it != stripPath_list.end()) {
            nameFound = true;
            // Construct the filename based on the name
            std::cout << "Found matching name in stripPath_list. Filename: " << name << " at index: " << i << std::endl;
            this->_read_one_pathObjFile(toolpathSet, path, name, i);
            this->_assign_waypointNormal_from_layerFacet(isoLayerSet, toolpathSet, name);
        }
        else {
            std::cout << "No file from layer_list is found in stripPath_list. Filename: " << name << " at index: " << i << std::endl;
        }
    }

}


void fileIO::_read_one_pathObjFile(PolygenMesh* toolpathSet, std::string path, std::string name, int givenId) {

    char filename[1024];
    sprintf(filename, "%s%s%s", path.c_str(), "/", name.c_str());

    QMeshPatch* waypoint = new QMeshPatch;
    //waypoint->SetIndexNo(toolpathSet->GetMeshList().GetCount()); //index begin from 0
    waypoint->SetIndexNo(givenId);
    toolpathSet->GetMeshList().AddTail(waypoint);
    waypoint->patchName = name;

    waypoint->inputStripOBJFile(filename);
}

void fileIO::_assign_waypointNormal_from_layerFacet(
    PolygenMesh* isoLayerSet, PolygenMesh* toolpathSet, std::string name) {

    //find the layer and tp patch with the same name
    QMeshPatch* _layer = NULL;
    QMeshPatch* _path = NULL;

    for (GLKPOSITION posMesh = isoLayerSet->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* each_layer = (QMeshPatch*)isoLayerSet->GetMeshList().GetNext(posMesh);

        if (each_layer->patchName == name) {
            _layer = each_layer;
            break;
        }
    }

    for (GLKPOSITION posMesh = toolpathSet->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* each_path = (QMeshPatch*)toolpathSet->GetMeshList().GetNext(posMesh);

        if (each_path->patchName == name) {
            _path = each_path;
            break;
        }
    }

    if (_layer == NULL || _path == NULL) std::cout << "Error: @ _assign_waypointNormal_from_layerFacet" << std::endl;
    _path->attached_Layer = _layer;

    //build a vector of QMeshFace to speed up the calculation of _get_faceNormal_with_idx
    std::vector<QMeshFace*> faceSet_of_layer;
    for (GLKPOSITION Pos = _layer->GetFaceList().GetHeadPosition(); Pos;) {
        QMeshFace* Face = (QMeshFace*)_layer->GetFaceList().GetNext(Pos);

        faceSet_of_layer.push_back(Face);
    }

    for (GLKPOSITION Pos = _path->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* curtNode = (QMeshNode*)_path->GetNodeList().GetNext(Pos);

        std::vector<int> edge_facet_idx;
        for (GLKPOSITION Pos = curtNode->GetEdgeList().GetHeadPosition(); Pos;) {
            QMeshEdge* neighEdge = (QMeshEdge*)curtNode->GetEdgeList().GetNext(Pos);

            edge_facet_idx.push_back(neighEdge->attached_face_idx);
        }

        Eigen::Vector3d curtNode_n = Eigen::Vector3d::Zero();
        for (int idx : edge_facet_idx) {

            //std::cout << idx << std::endl;

            Eigen::Vector3d n = this->_get_faceNormal_with_idx(faceSet_of_layer, idx);

            curtNode_n += n;
        }
        curtNode_n.normalize();
        curtNode->SetNormal(curtNode_n[0], curtNode_n[1], curtNode_n[2]);
    }

}

Eigen::Vector3d fileIO::_get_faceNormal_with_idx(const std::vector<QMeshFace*> &faceSet_of_layer, int idx) {

    QMeshFace* attachedFace = faceSet_of_layer[idx - 1];
    if (attachedFace == NULL) 
        std::cout << "Error: @ _get_faceNormal_with_idx, reason: face id start from 1" << std::endl;

    Eigen::Vector3d n;
    //attachedFace->CalPlaneEquation();
    attachedFace->GetNormal(n[0], n[1], n[2]);
    n.normalize();
    return -n; //the face normal is downward, reverse it
}

void fileIO::output_StripToolpath(PolygenMesh* toolpathSet, std::string path) {

    this->_remove_allFile_in_Dir(path);

    std::string Path_dir;

    for (GLKPOSITION posMesh = toolpathSet->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* each_path = (QMeshPatch*)toolpathSet->GetMeshList().GetNext(posMesh);

        std::string sPath_name = each_path->patchName.substr(0, each_path->patchName.find_last_of('.')) + ".txt";

        Path_dir = path + "/" + sPath_name;

        this->_output_One_stripToolpath(each_path, Path_dir);
        std::cout << "* " << std::endl;
    }

}

void fileIO::_output_One_stripToolpath(QMeshPatch* eachPath, std::string path) {

    Eigen::Vector3d offset = { 70.0,0.0,0.0 }; // only for ncc3_half

    double pp[3],pn[3];
    std::ofstream file_(path);

    for (GLKPOSITION posNode = eachPath->GetNodeList().GetHeadPosition(); posNode != nullptr;) {
        QMeshNode* node = (QMeshNode*)eachPath->GetNodeList().GetNext(posNode);

        node->GetCoord3D(pp[0], pp[1], pp[2]);
        node->GetNormal(pn[0], pn[1], pn[2]);

        file_
            << pp[0] + offset[0] << " " << pp[1] + offset[1] << " " << pp[2] + offset[2] << " "
            << pn[0] << " " << pn[1] << " " << pn[2] <<
            " " << node->m_printTan[0] << " " << node->m_printTan[1] << " " << node->m_printTan[2] <<
            " " << 0.0547 <<
            std::endl;
     
    }
    file_.close();

    std::cout << "the offset of position is " << offset.transpose() << std::endl;
}