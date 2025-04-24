#include "directionFieldOpt.h"
#include "meshOperation.h"


void directionFieldOpt::initialize(PolygenMesh* isoLayerSet, std::string model_name) {

	m_isoLayerSet = isoLayerSet;

	meshOperation* meshOperator = new meshOperation();
	meshOperator->initialiseIndices(m_isoLayerSet);
	delete meshOperator;

	//the a vector which shows the constraints number of each layer
    for (GLKPOSITION PatchPos = m_isoLayerSet->GetMeshList().GetHeadPosition(); PatchPos;) {
        QMeshPatch* Patch = (QMeshPatch*)m_isoLayerSet->GetMeshList().GetNext(PatchPos);

        int _constraintNumber = 0;

        for (GLKPOSITION facePos = Patch->GetFaceList().GetHeadPosition(); facePos;) {
            QMeshFace* thisFace = (QMeshFace*)Patch->GetFaceList().GetNext(facePos);

            if (thisFace->isConstraint) {

                _constraintNumber++;
            }
        }
		Patch->constraintNumber = _constraintNumber;
    }

	this->_setOptParameters(model_name);
}

void directionFieldOpt::_setOptParameters(std::string model_name) {

	//Weighting of stress = stressWeight_ratio * pow(10,[0,-stressWeight_range])

	if (model_name == "TshapeBracketNew") {
		
		b_Weighting = { 5.0 * pow(10,0), 1.0 * pow(10,-3) };
		stressWeight_range = 2.0;
		stressWeight_ratio = 0.1;
		/*b_Weighting = { 5.0 * pow(10,-15), 1.0 * pow(10,-16) };
		stressWeight_range = 2.0;
		stressWeight_ratio = 0.1;*/
	}
	if (model_name == "bladeSmall") {
		b_Weighting = { 5.0 * pow(10,-6), 1.0 * pow(10,-8) };
		stressWeight_range = 2.0;
		stressWeight_ratio = 0.0005;
	}
	if (model_name == "ncc2") {

		//b_Weighting = { 8.0 * pow(10,0), 1.0 * pow(10,-2) };
		b_Weighting = { 8.0 * pow(10,-10), 1.0 * pow(10,-20) };
		stressWeight_range = 4.0;
		stressWeight_ratio = 0.1;
	}
	if (model_name == "GEBracketModified"){
		b_Weighting = { 8.0 * pow(10,1), 1.0 * pow(10,-2) };
		stressWeight_range = 4.0;
		stressWeight_ratio = 0.0005;
	}
	if (model_name == "bridge") {
		b_Weighting = { 1.0 * pow(10,-2), 1.0 * pow(10,-3) };
		stressWeight_range = 2.0;
		stressWeight_ratio = 0.0005;
	}
	if (model_name == "ncc3_half") {

		b_Weighting = { 8.0 * pow(10,-1), 1.0 * pow(10,-2) };
		stressWeight_range = 4.0;
		stressWeight_ratio = 0.001;
	}
	if (model_name == "frameBA" || model_name == "curvedFrame") {

		b_Weighting = { 8.0 * pow(10, 0), 8.0 * pow(10, 0) };
		stressWeight_range = 4.0;
		stressWeight_ratio = 0.01;
	}
	if (model_name == "shelf_ty") {

		b_Weighting = { 8.0 * pow(10, -5), 8.0 * pow(10, -5) };
		stressWeight_range = 1.0;
		stressWeight_ratio = 0.01;
	}

}

void directionFieldOpt::run() {

	for (GLKPOSITION PatchPos = m_isoLayerSet->GetMeshList().GetHeadPosition(); PatchPos;) {
		QMeshPatch* Patch = (QMeshPatch*)m_isoLayerSet->GetMeshList().GetNext(PatchPos);

		this->_installCoordinateSystem(Patch);
		this->_installSlopes(Patch);
		this->_classifyBoundary(Patch);
		this->_smoothVectorField(Patch);
		this->_resolveSingularity(Patch);
		this->_calAngleDiff_dirctionField_vs_stressField(Patch);
		Patch->draw_directionField_onFace = true; //draw the directionField on face
		this->_getNodeDirectionField(Patch);
		Patch->draw_Field_onFace = false;//draw the directionField on Node
	}

}

//prev version
//void directionFieldOpt::_installCoordinateSystem(QMeshPatch* materialMesh){
//
//	for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
//		QMeshFace* thisFace = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);
//		thisFace->CalPlaneEquation();
//
//		double nx, ny, nz;
//		thisFace->GetNormal(nx, ny, nz);
//
//		Eigen::Vector3d targetVector(nx, ny, nz);
//		Eigen::Vector3d xDir(1, 0, 0);
//		Eigen::Vector3d yDir(0, 1, 0);
//		Eigen::Vector3d zDir(0, 0, 1);
//
//		Eigen::Quaterniond rotat;
//		rotat.setFromTwoVectors(zDir, targetVector);
//		zDir = rotat * zDir;
//		assert((zDir[2] - targetVector[2]) < 1e-4);
//		xDir = rotat * xDir;
//		yDir = rotat * yDir;
//
//		thisFace->xDir = xDir;
//		thisFace->yDir = yDir;
//		thisFace->zDir = zDir;
//
//		//std::cout << thisFace->xDir << std::endl;
//	}
//
//}

void directionFieldOpt::_installCoordinateSystem(QMeshPatch* materialMesh)
{

	Eigen::Vector3d xDir_base(1, 0, 0);
	Eigen::Vector3d yDir_base(0, 1, 0);
	Eigen::Vector3d zDir_base(0, 0, 0);
	Eigen::Vector3d zDir_orig(0, 0, 1);

	int count = 0;
	for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
		QMeshFace* thisFace = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);
		thisFace->CalPlaneEquation();

		double nx, ny, nz;
		thisFace->GetNormal(nx, ny, nz);

		Eigen::Vector3d targetVector(nx, ny, nz);
		zDir_base = zDir_base + targetVector;
		zDir_base.stableNormalize();
		count++;
	}

	zDir_base.stableNormalize();

	Eigen::Quaterniond rotat_base;
	rotat_base.setFromTwoVectors(zDir_orig, zDir_base);
	zDir_orig = rotat_base * zDir_orig;
	assert((zDir_base[2] - zDir_orig[2]) < 1e-4);
	xDir_base = rotat_base * xDir_base;
	yDir_base = rotat_base * yDir_base;

	for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
		QMeshFace* thisFace = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);
		thisFace->CalPlaneEquation();

		double nx, ny, nz;
		thisFace->GetNormal(nx, ny, nz);

		Eigen::Vector3d targetVector(nx, ny, nz);
		/*Eigen::Vector3d xDir(1, 0, 0);
		Eigen::Vector3d yDir(0, 1, 0);
		Eigen::Vector3d zDir(0, 0, 1);*/
		Eigen::Vector3d xDir = xDir_base;
		Eigen::Vector3d yDir = yDir_base;
		Eigen::Vector3d zDir = zDir_base;


		Eigen::Quaterniond rotat;
		rotat.setFromTwoVectors(zDir, targetVector);
		zDir = rotat * zDir;
		assert((zDir[2] - targetVector[2]) < 1e-4);
		xDir = rotat * xDir;
		yDir = rotat * yDir;

		thisFace->xDir = xDir;
		thisFace->yDir = yDir;
		thisFace->zDir = zDir;
	}

}

void directionFieldOpt::_installSlopes(QMeshPatch* materialMesh){

	for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
		QMeshFace* thisFace = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);

		if (!thisFace->isConstraint) continue;

		Eigen::Vector3d dirVec(thisFace->vectorDir[0], thisFace->vectorDir[1], thisFace->vectorDir[2]);
		thisFace->xDir.stableNormalize();
		thisFace->yDir.stableNormalize();
		dirVec.stableNormalize();
		double cos_angle = thisFace->xDir.dot(dirVec);
		double sin_angle = thisFace->yDir.dot(dirVec);
		double slope;
		double angle = atan2(sin_angle, cos_angle);
		slope = cos(2 * angle);
		thisFace->slope1 = slope;
		slope = sin(2 * angle);
		thisFace->slope2 = slope;
	}
}

void directionFieldOpt::_classifyBoundary(QMeshPatch* materialMesh) {
	std::cout << "Classifying boundary....\n";
	std::set<QMeshEdge*> boundarySet;

	for (GLKPOSITION edgePos = materialMesh->GetEdgeList().GetHeadPosition(); edgePos; ) {
		QMeshEdge* thisEdge = (QMeshEdge*)materialMesh->GetEdgeList().GetNext(edgePos);

		if (!thisEdge->boundary) continue;
		boundarySet.insert(thisEdge);
	}

	if (boundarySet.empty()) {
		std::cout << "No boundary edges found.\n";
		return;
	}

	QMeshEdge* startEdge = *boundarySet.begin();
	double angleChange = 0;
	std::set<QMeshEdge*> tempSet;
	tempSet.insert(startEdge);
	QMeshEdge* nextEdge = nullptr;

	while (!boundarySet.empty()) {
		boundarySet.erase(startEdge);

		QMeshNode* checkNode = startEdge->GetEndPoint();
		for (int i = 0; i < checkNode->GetEdgeNumber(); i++) {
			nextEdge = checkNode->GetEdgeRecordPtr(i + 1);
			if (nextEdge->boundary && nextEdge != startEdge) {
				break;
			}
		}

		if (nextEdge == nullptr) {
			std::cout << "Next edge not found.\n";
			return;
		}

		bool isIn = boundarySet.find(nextEdge) != boundarySet.end();
		tempSet.insert(nextEdge);
		angleChange += calculateAngleChange(startEdge, nextEdge);

		if (isIn) {
			startEdge = nextEdge;
			startEdge->visited = true;
		}
		else {
			classifyEdges(tempSet, angleChange);
			tempSet.clear();
			angleChange = 0;
			if (!boundarySet.empty()) {
				startEdge = *boundarySet.begin();
			}
		}
	}
}

double directionFieldOpt::calculateAngleChange(QMeshEdge* startEdge, QMeshEdge* nextEdge) {
	double xs, ys, zs, xe, ye, ze;
	startEdge->GetStartPoint()->GetCoord3D(xs, ys, zs);
	startEdge->GetEndPoint()->GetCoord3D(xe, ye, ze);
	Eigen::Vector3d startVector(xe - xs, ye - ys, ze - zs);

	nextEdge->GetStartPoint()->GetCoord3D(xs, ys, zs);
	nextEdge->GetEndPoint()->GetCoord3D(xe, ye, ze);
	Eigen::Vector3d nextVector(xe - xs, ye - ys, ze - zs);

	QMeshFace* thisFace = (QMeshFace*)startEdge->GetLeftFace();
	thisFace->GetNormal(xs, ys, zs);
	Eigen::Vector3d normalVec(xs, ys, zs);
	normalVec.stableNormalize();

	nextVector = nextVector - nextVector.dot(normalVec) * normalVec;

	Eigen::Vector3d crossProd = startVector.cross(nextVector);
	crossProd = crossProd / startVector.stableNorm();
	crossProd = crossProd / nextVector.stableNorm();

	double sinA = normalVec.dot(crossProd);
	double cosA = startVector.dot(nextVector) / (startVector.stableNorm() * nextVector.stableNorm());
	return atan2(sinA, cosA);
}

void directionFieldOpt::classifyEdges(std::set<QMeshEdge*>& edgeSet, double angleChange) {
	for (auto edge : edgeSet) {
		if (angleChange > 0) {
			edge->isOutBound = true;
			edge->GetLeftFace()->outBFace = true;
		}
		else {
			edge->isInBound = true;
			edge->GetLeftFace()->inBFace = true;
		}
	}
}


void directionFieldOpt::_smoothVectorField(QMeshPatch* materialMesh){

	int num = materialMesh->GetFaceNumber();

	std::cout << materialMesh->constraintNumber << std::endl;

	Eigen::VectorXd m10 = Eigen::VectorXd::Zero(num + materialMesh->constraintNumber);
	Eigen::VectorXd m20 = Eigen::VectorXd::Zero(num + materialMesh->constraintNumber);

	Eigen::SparseMatrix<double> M;
	this->_buildLaplacianMatrix(materialMesh, M);

	this->_createSolutionVector(materialMesh, 1, m10);
	this->_createSolutionVector(materialMesh, 2, m20);

	Eigen::VectorXd m1_sol(num);
	Eigen::VectorXd m2_sol(num);
	//std::cout << "Check2\n";


	std::cout << "Initialising Sparse Solver for vector smoothing..\n";
	Eigen::SparseLU <Eigen::SparseMatrix<double>> SolverN;// (PardisoLU/SparseLU)
	Eigen::SparseMatrix<double> S = M.transpose() * M;
	SolverN.analyzePattern(S);
	SolverN.factorize(S);
	if (SolverN.info() != Eigen::Success)
		std::cout << "error here: factorize fail!" << std::endl;
	SolverN.compute(S);
	//std::cout << "Smoothing 1\n";
	m10 = M.transpose() * m10;
	m1_sol = SolverN.solve(m10);

	//std::cout << "Smoothing 2\n";
	m20 = M.transpose() * m20;
	m2_sol = SolverN.solve(m20);


	for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
		QMeshFace* face = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);
		int index = face->FieldIndexNumber;


		Eigen::Vector2d repVec(m1_sol(index), m2_sol(index));
		repVec.stableNormalize();

		//std::cout << repVec.transpose() << std::endl;

		double theta_ = atan2(repVec.y(), repVec.x());
		theta_ = 0.5 * theta_;
		double c1 = cos(theta_);
		double s1 = sin(theta_);

		if ((abs(theta_ - 1.5708) < 1e-2) || (abs(-theta_ - 1.5708) < 1e-2)) {
			face->singularWarning = true;

			//mark the singularWarning flag for nodes
			for (int _i = 0; _i < 3; _i++) {
				face->GetNodeRecordPtr(_i)->singularWarning = true;
			}//end

			//std::cout << "singularity....\n";
		}

		Eigen::Vector3d newNormal = s1 * face->yDir + c1 * face->xDir;
		newNormal.stableNormalize();
		//newNormal.stableNormalize();


		face->vectorDir[0] = newNormal.x();
		face->vectorDir[1] = newNormal.y();
		face->vectorDir[2] = newNormal.z();

		//std::cout << face->vectorDir[0] << std::endl;

	}

}

void directionFieldOpt::_buildLaplacianMatrix(QMeshPatch* materialMesh, 
	Eigen::SparseMatrix<double>& LaplacianMatrix) {

	this->_normalizeWeights(materialMesh);
	int num = materialMesh->GetFaceNumber();

	Eigen::SparseMatrix<double> diagonalMultiplier;

	LaplacianMatrix.resize(num + materialMesh->constraintNumber, num);
	LaplacianMatrix.reserve(Eigen::VectorXd::Constant(num + materialMesh->constraintNumber, 10));

	diagonalMultiplier.resize(num + materialMesh->constraintNumber, num + materialMesh->constraintNumber);
	diagonalMultiplier.setIdentity();

	//std::cout << "I am here 1\n";

	for (GLKPOSITION edgePos = materialMesh->GetEdgeList().GetHeadPosition(); edgePos;) {
		QMeshEdge* edge = (QMeshEdge*)materialMesh->GetEdgeList().GetNext(edgePos);

		QMeshFace* f1 = edge->GetRightFace();
		QMeshFace* f2 = edge->GetLeftFace();

		if (f1 == nullptr || f2 == nullptr) continue;

		int i = f1->FieldIndexNumber;
		int j = f2->FieldIndexNumber;

		LaplacianMatrix.insert(i, j) = -1;
		LaplacianMatrix.insert(j, i) = -1;

		if (LaplacianMatrix.coeff(i, i) == 0) {
			LaplacianMatrix.insert(i, i) = 1;
		}
		else {
			LaplacianMatrix.coeffRef(i, i)++;
		}

		if (LaplacianMatrix.coeff(j, j) == 0) {
			LaplacianMatrix.insert(j, j) = 1;
		}
		else {
			LaplacianMatrix.coeffRef(j, j)++;
		}
	}

	//std::cout << "I am here 3\n";
	int pivotInsertCount = 0;
	for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
		QMeshFace* face = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);
		if (!face->isConstraint) continue;

		int i = face->FieldIndexNumber;
		int indexx = pivotInsertCount + num;
		LaplacianMatrix.insert(indexx, i) = 0.001;

		// the weighting strategy
		if (face->isBoundaryConstraint) {
			if (face->outBFace) {
				LaplacianMatrix.coeffRef(indexx, i) = b_Weighting[0]; // 0
			}
			else {
				LaplacianMatrix.coeffRef(indexx, i) = b_Weighting[1]; // -2
			}
		}
		else {
			double weight_ = std::pow(10, face->stressMag);
			LaplacianMatrix.coeffRef(indexx, i) = weight_;
		}
		pivotInsertCount++;
	}

	for (int l = 0; l < num; l++) {
		if (LaplacianMatrix.coeff(l, l) < 1) std::cout << "DANGER!!!!!\n";
		diagonalMultiplier.coeffRef(l, l) = 1 / LaplacianMatrix.coeff(l, l);
	}

	LaplacianMatrix = diagonalMultiplier * LaplacianMatrix;
}

void directionFieldOpt::_normalizeWeights(QMeshPatch* materialMesh){

	std::cout << "adjusting weights...\n";
	double minStress = 9e10;
	double maxStress = -9e10;
	for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
		QMeshFace* face = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);
		if (!face->isConstraint) continue;
		if (face->isBoundaryConstraint) continue;
		double thisStress = face->stressMag;

		if (thisStress < minStress) minStress = thisStress;
		if (thisStress > maxStress) maxStress = thisStress;
	}

	//std::cout << "Stress\n";
	//std::cout << maxStress << std::endl;
	//std::cout << minStress << std::endl;

	
	if (maxStress != minStress) {

		for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
			QMeshFace* face = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);
			if (!face->isConstraint) continue;
			if (face->isBoundaryConstraint) continue;
			double thisStress = face->stressMag;

			face->stressMag =
				stressWeight_ratio * (-stressWeight_range + stressWeight_range * (thisStress - minStress) / (maxStress - minStress));
		}
	}
	else {

		for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
			QMeshFace* face = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);
			if (!face->isConstraint) continue;
			if (face->isBoundaryConstraint) continue;
			double thisStress = face->stressMag;

			face->stressMag = 0;
		}
	}

}

void directionFieldOpt::_createSolutionVector(QMeshPatch* materialMesh, int column, Eigen::VectorXd& u)
{
	int num = materialMesh->GetFaceNumber();

	std::cout << "Creating vector " << column << std::endl;
	int pivotInsertCount = 0;
	for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
		QMeshFace* face = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);
		if (!face->isConstraint) continue;
		int index = face->FieldIndexNumber;
		//if (face->slope < 0) std::cout << "Again error\n";

		//1 -> Xlocal ; 2 -> Ylocal
		//the weighting strategy
		if (column == 1) {
			if (face->isBoundaryConstraint) {
				if (face->outBFace) {
					u(num + pivotInsertCount) = b_Weighting[0] * face->slope1; //0
				}
				else {
					u(num + pivotInsertCount) = b_Weighting[1] * face->slope1; //-2
				}
			}
			else {
				double weight_ = pow(10, face->stressMag);
				u(num + pivotInsertCount) = weight_ * face->slope1;
			}
		}
		else {
			if (face->isBoundaryConstraint) {
				if (face->outBFace) {
					u(num + pivotInsertCount) = b_Weighting[0] * face->slope2; //0
				}
				else {
					u(num + pivotInsertCount) = b_Weighting[1] * face->slope2; //-2
				}
			}
			else {
				double weight_ = pow(10, face->stressMag);
				u(num + pivotInsertCount) = weight_ * face->slope2;
			}
		}
		pivotInsertCount++;

	}


}

void directionFieldOpt::_resolveSingularity(QMeshPatch* materialMesh){

	for (int k = 0; k < 30; k++) {
		for (GLKPOSITION edgePos = materialMesh->GetEdgeList().GetHeadPosition(); edgePos;) {
			QMeshEdge* edge = (QMeshEdge*)materialMesh->GetEdgeList().GetNext(edgePos);

			QMeshFace* LeftFace = edge->GetLeftFace();
			QMeshFace* RightFace = edge->GetRightFace();

			if (!LeftFace || !RightFace) {
				continue;
			}


			if (LeftFace->singularWarning) {

				if (!RightFace->singularWarning) {
					//std::cout << "r...\n";
					Eigen::Vector3d leftVec(LeftFace->vectorDir[0], LeftFace->vectorDir[1], LeftFace->vectorDir[2]);
					Eigen::Vector3d rightVec(RightFace->vectorDir[0], RightFace->vectorDir[1], RightFace->vectorDir[2]);

					if (leftVec.dot(rightVec) < -0.95) {
						//std::cout << "resolving...\n";
						LeftFace->vectorDir[0] = -LeftFace->vectorDir[0];
						LeftFace->vectorDir[1] = -LeftFace->vectorDir[1];
						LeftFace->vectorDir[2] = -LeftFace->vectorDir[2];

					}
					//std::cout << leftVec.dot(rightVec) << std::endl;
					LeftFace->singularWarning = false;
				}
			}

			if (RightFace->singularWarning) {
				if (!LeftFace->singularWarning) {
					//std::cout << "r...\n";
					Eigen::Vector3d leftVec(LeftFace->vectorDir[0], LeftFace->vectorDir[1], LeftFace->vectorDir[2]);
					Eigen::Vector3d rightVec(RightFace->vectorDir[0], RightFace->vectorDir[1], RightFace->vectorDir[2]);

					if (rightVec.dot(leftVec) < -0.95) {
						//std::cout << "resolving...\n";
						RightFace->vectorDir[0] = -RightFace->vectorDir[0];
						RightFace->vectorDir[1] = -RightFace->vectorDir[1];
						RightFace->vectorDir[2] = -RightFace->vectorDir[2];

					}
					//std::cout << leftVec.dot(rightVec) << std::endl;
					RightFace->singularWarning = false;
				}
			}
		}
	}
}

void directionFieldOpt::_getNodeDirectionField(QMeshPatch* materialMesh) {

	for (GLKPOSITION nodePos = materialMesh->GetNodeList().GetHeadPosition(); nodePos;) {
		QMeshNode* node = (QMeshNode*)materialMesh->GetNodeList().GetNext(nodePos);

		std::vector<Eigen::Vector3d> neiFace_DirectionField;

		for (GLKPOSITION neiFacePos = node->GetFaceList().GetHeadPosition(); neiFacePos;) {
			QMeshFace* neiFace = (QMeshFace*)node->GetFaceList().GetNext(neiFacePos);

			neiFace_DirectionField.push_back(Eigen::Vector3d(neiFace->vectorDir[0], neiFace->vectorDir[1], neiFace->vectorDir[2]));
		}

		/*std::cout << "output a neiFace_DirectionField: " << std::endl;
		for (const auto& direction : neiFace_DirectionField) {
			std::cout << "[" << direction[0] << ", " << direction[1] << ", " << direction[2] << "]" << std::endl;
		}*/

		// Step 1: Calculate angles between all pairs of vectors
		int n = neiFace_DirectionField.size();
		if (n <= 1) continue; // No need to process if there's only one or no direction

		std::vector<double> maxAngle(n, 0.0); // Store the maximum angle for each vector

		for (int i = 0; i < n; ++i) {
			for (int j = i + 1; j < n; ++j) {
				double dotProduct = neiFace_DirectionField[i].dot(neiFace_DirectionField[j]);
				double angle = acos(dotProduct / (neiFace_DirectionField[i].norm() * neiFace_DirectionField[j].norm()));
				maxAngle[i] = std::max(maxAngle[i], angle);
				maxAngle[j] = std::max(maxAngle[j], angle);
			}
		}

		// Step 2: Find the vector with the maximum deviation
		int maxIndex = std::distance(maxAngle.begin(), std::max_element(maxAngle.begin(), maxAngle.end()));
		neiFace_DirectionField.erase(neiFace_DirectionField.begin() + maxIndex);

		// Step 3: Compute the average direction
		Eigen::Vector3d averageDirection(0, 0, 0);
		for (const auto& vec : neiFace_DirectionField) {
			averageDirection += vec;
		}
		averageDirection.normalize();

		//std::cout << "Average direction: [" << averageDirection[0] << ", " << averageDirection[1] << ", " << averageDirection[2] << "]" << std::endl;
		//for singularWarning nodes, directely give the first vector of neiFace_DirectionField
		if (node->singularWarning) {
			node->directionField = neiFace_DirectionField[0];
			//node->SetNormal(neiFace_DirectionField[0][0], neiFace_DirectionField[0][1], neiFace_DirectionField[0][2]);
		}
		else {
			node->directionField = averageDirection;
			//node->SetNormal(averageDirection[0], averageDirection[1], averageDirection[2]);
		}

	}
}

//void directionFieldOpt::_calAngleDiff_dirctionField_vs_stressField(QMeshPatch* materialMesh) {
//
//	if (!materialMesh->isStressLayer) return; // If it's not a stress layer, exit the function.
//
//	double avgDiff = 0.0;  // To accumulate angle differences
//	double maxDiff = 0.0;  // To record the maximum angle difference
//	int _num = 0;  // Counter for the number of faces with stress field data
//	std::vector<double> angleDiffs; // To store all angle differences
//
//	for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
//		QMeshFace* face = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);
//
//		if (face->stessField_flag) { // Check if the face has stress field data
//
//			Eigen::Vector3d face_DirectionField = face->vectorDir;  // Direction field vector
//			Eigen::Vector3d face_stressField = face->stessField_vector;  // Stress field vector
//
//			double dotProduct = face_DirectionField.dot(face_stressField);  // Dot product of the two vectors
//
//			double normA = face_DirectionField.norm();  // Magnitude of direction field vector
//			double normB = face_stressField.norm();  // Magnitude of stress field vector
//
//			double cosTheta = dotProduct / (normA * normB);  // Cosine of the angle between the vectors
//
//			double angleRad = acos(cosTheta);  // Convert cosine to angle in radians
//			double angleDeg = angleRad * 180.0 / 3.141592654;  // Convert angle to degrees
//
//			if (angleDeg > 90.0) angleDeg = 180.0 - angleDeg;  // Ensure angle is in [0, 90] degrees range
//
//			avgDiff += angleDeg;  // Accumulate angle differences
//			angleDiffs.push_back(angleDeg);  // Store the current angle difference
//
//			if (angleDeg > maxDiff) maxDiff = angleDeg;  // Update maximum angle difference
//			_num++;  // Increment the counter
//		}
//	}
//
//	// Compute variance
//	double variance = 0.0;
//	double meanDiff = avgDiff / _num;  // Calculate the average angle difference
//	for (double diff : angleDiffs) {
//		variance += (diff - meanDiff) * (diff - meanDiff);  // Sum of squared differences
//	}
//	variance /= _num;  // Final variance calculation
//
//	// Output the results: average angle difference, maximum angle difference, and variance
//	std::cout << "The AngleDiff between direction and stress Field of " << materialMesh->patchName
//		<< " is " << meanDiff << ", maxDiff: " << maxDiff << ", variance: " << variance << std::endl;
//}

void directionFieldOpt::_calAngleDiff_dirctionField_vs_stressField(QMeshPatch* materialMesh) {

	if (!materialMesh->isStressLayer) return;

	double avgDiff = 0.0;
	double maxDiff = 0.0;
	int _num = 0;
	std::vector<double> angleDiffs;

	// Create a unique filename using patch name or a unique identifier
	std::stringstream ss;
	//ss << "../DataSet/FigureData/angle_differences_planar_" << materialMesh->patchName << ".txt";
	ss << "../DataSet/FigureData/angle_differences" << materialMesh->patchName << ".txt";
	std::string filename = ss.str();

	// Open the file to store the angle differences
	std::ofstream outFile(filename);

	// Check if the file opened successfully
	if (!outFile.is_open()) {
		std::cerr << "Failed to open the file: " << filename << std::endl;
		return;
	}

	for (GLKPOSITION facePos = materialMesh->GetFaceList().GetHeadPosition(); facePos;) {
		QMeshFace* face = (QMeshFace*)materialMesh->GetFaceList().GetNext(facePos);

		if (face->stessField_flag) {

			Eigen::Vector3d face_DirectionField = face->vectorDir;
			//face_DirectionField = { 1,0,0 };
			Eigen::Vector3d face_stressField = face->stessField_vector;

			double dotProduct = face_DirectionField.dot(face_stressField);

			double normA = face_DirectionField.norm();
			double normB = face_stressField.norm();

			double cosTheta = dotProduct / (normA * normB);

			double angleRad = acos(cosTheta);
			double angleDeg = angleRad * 180.0 / 3.141592654;

			if (angleDeg > 90.0) angleDeg = 180.0 - angleDeg;

			avgDiff += angleDeg;
			angleDiffs.push_back(angleDeg);

			// Write each angle difference to the file
			outFile << angleDeg << std::endl;

			if (angleDeg > maxDiff) maxDiff = angleDeg;
			_num++;
		}
	}

	// Close the file after writing
	outFile.close();

	// Calculate variance
	double variance = 0.0;
	double meanDiff = avgDiff / _num;
	for (double diff : angleDiffs) {
		variance += (diff - meanDiff) * (diff - meanDiff);
	}
	variance /= _num;

	// Output the average, max, and variance
	std::cout << "The AngleDiff between direction and stress Field of " << materialMesh->patchName
		<< " is " << meanDiff << ", maxDiff: " << maxDiff << ", variance: " << variance << std::endl;
}
