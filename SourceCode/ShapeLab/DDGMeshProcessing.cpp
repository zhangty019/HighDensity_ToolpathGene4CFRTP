#include "DDGMeshProcessing.h"
#include "Types.h"

void DDGMeshProcessing::meshProcessing()
{
	//mesh.computeSmoothestSection();
	mesh.parameterize();
	mesh.glueParameterization();

	double meanEdgeLength = 0.;
	for (auto e = mesh.edges.begin(); e != mesh.edges.end(); e++)
	{
		meanEdgeLength += e->length();
	}
	meanEdgeLength /= (double)mesh.edges.size();

	double k = mesh.fieldDegree;
	for (auto& vit : mesh.vertices) {
		DDG::Vector Z = vit.fieldVector(k, 1).unit() * meanEdgeLength;
		std::vector<double> gradient = { Z.x, Z.y, Z.z };
		perVerticeGradient.push_back(gradient);
	}
}

QMeshPatch* DDGMeshProcessing::toQMeshPatch() {

	int nv = mesh.vertices.size();
	int nf = mesh.faces.size();

	float* v = new float[mesh.vertices.size() * 3];
	uint* f = new uint[mesh.faces.size() * 3];

	std::vector<float> fn;
	std::vector<float> ft;
	std::vector<double> stripeIndices;
	std::vector<double> fieldIndices;

	for (int i = 0; i < nv; i++) {
		v[i * 3] = (float)mesh.vertices[i].position.x;
		v[i * 3 + 1] = (float)mesh.vertices[i].position.y;
		v[i * 3 + 2] = (float)mesh.vertices[i].position.z;
	}

	for (int i = 0; i < mesh.faces.size(); i++) {
		auto face = mesh.faces[i];
		f[i * 3] = (uint)face.he->vertex->index;
		f[i * 3 + 1] = (uint)face.he->next->vertex->index;
		f[i * 3 + 2] = (uint)face.he->next->next->vertex->index;

		fn.push_back(face.normal().x);
		fn.push_back(face.normal().y);
		fn.push_back(face.normal().z);

		//record the paramIndex(stripeIndices)
		stripeIndices.push_back(face.paramIndex[0]);
		stripeIndices.push_back(face.paramIndex[1]);
		//record the fieldIndex
		fieldIndices.push_back(face.fieldIndex(2.));

		DDG::HalfEdgeIter he = face.he;

		for (int j = 0; j < 3; j++)
		{
			DDG::Complex z = he->texcoord / (2. * M_PI) + DDG::Complex(.5, .5);
			ft.push_back(z.re);
			ft.push_back(z.im);
			he = he->next;
		}
	}

	//build the new mesh
	QMeshPatch* newMesh = new QMeshPatch;
	newMesh->constructionFromVerFaceTable(nv, v, nf, f);

	int faceIndex = 0;
	GLKPOSITION Pos;
	for (Pos = newMesh->GetFaceList().GetHeadPosition(); Pos != NULL; faceIndex++) {
		QMeshFace* face = (QMeshFace*)(newMesh->GetFaceList().GetNext(Pos));
		for (int i = 0; i < 3; i++) {
			face->SetVerticeNormal(i, fn[3 * faceIndex], fn[3 * faceIndex + 1], fn[3 * faceIndex + 2]);
			face->SetTextureCoord(i, ft[6 * faceIndex + 2 * i], ft[6 * faceIndex + 2 * i + 1]);
		}
		// install the stripeIndices and fieldIndex
		face->stripeIndices[0] = stripeIndices[2 * faceIndex];
		face->stripeIndices[1] = stripeIndices[2 * faceIndex + 1];
		face->fieldIndices = fieldIndices[faceIndex];

		/*std::cout << "Face " << faceIndex << " stripeIndices: ["
			<< face->stripeIndices[0] << ", " << face->stripeIndices[1] << "]"
			<< ", fieldIndex: " << face->fieldIndex << std::endl;*/

			/*std::cout << "Face " << faceIndex << " boundary: " << face->isBoundaryFace() << std::endl;*/
	}
	return newMesh;
}


void DDGMeshProcessing::setInitField(std::vector<std::vector<double>>& initField) {
	for (auto& vit : mesh.vertices) {
		int idx = vit.index;
		DDG::HalfEdgeIter he = vit.he;
		std::vector<std::pair<int, double>> vv_idx_angle;

		DDG::Vector initFieldVector(initField[idx][0], initField[idx][1], initField[idx][2]);
		if (initFieldVector.norm2() < 1e-8) {
			vit.directionField = DDG::Complex(0, 0);
			continue;
		}

		//vector_in_face <==> check Scalar Triple Product < eps && \theta in 
		// first project to plane
		// check in triangle
		if (vit.onBoundary()) {
			// check the first he
			while (!he->flip->onBoundary) {
				he = he->flip->next;
			}
			vit.he = he;
		}
		do
		{
			int idx = he->flip->vertex->index;
			double angle = he->angularCoordinate;
			vv_idx_angle.push_back(std::make_pair(idx, angle));
			he = he->next->next->flip;
		} while (he != vit.he && (!he->onBoundary));

		if (he->onBoundary) {
			int idx = he->flip->vertex->index;
			double angle = he->angularCoordinate;
			vv_idx_angle.push_back(std::make_pair(idx, angle));
		}
		double angleSacle = 2 * M_PI / vit.angleSum();
		double angle = 0;
		bool findit = false;
		if (vit.onBoundary()) {
			for (int i = 0; i < vv_idx_angle.size() - 1; i++) {
				//https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
				//https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
				auto P1 = vit.position;
				auto P2 = mesh.vertices[vv_idx_angle[i].first].position;
				auto P3 = mesh.vertices[vv_idx_angle[i + 1].first].position;

				DDG::Vector v0 = P2 - P1;
				DDG::Vector v1 = P3 - P1;
				DDG::Vector v2 = initFieldVector * 1000;

				double d00 = DDG::dot(v0, v0);
				double d01 = DDG::dot(v0, v1);
				double d11 = DDG::dot(v1, v1);
				double d20 = DDG::dot(v2, v0);
				double d21 = DDG::dot(v2, v1);

				double denom = d00 * d11 - d01 * d01;
				double v = (d11 * d20 - d01 * d21) / denom;
				double w = (d00 * d21 - d01 * d20) / denom;
				double u = 1.0f - v - w;

				if (v >= 0 && w >= 0 && u <= 1) {
					DDG::Vector P = u * P1 + v * P2 + w * P3;
					angle = DDG::angle(P2 - P1, P - P1) + vv_idx_angle[i].second;
					findit = true;
				}
			}
			if (findit == false) {
				for (int i = 0; i < vv_idx_angle.size() - 1; i++) {
					//https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
					//https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
					auto P1 = vit.position;
					auto P2 = mesh.vertices[vv_idx_angle[i].first].position;
					auto P3 = mesh.vertices[vv_idx_angle[i + 1].first].position;

					DDG::Vector v0 = P2 - P1;
					DDG::Vector v1 = P3 - P1;
					DDG::Vector v2 = -initFieldVector * 1000;

					double d00 = DDG::dot(v0, v0);
					double d01 = DDG::dot(v0, v1);
					double d11 = DDG::dot(v1, v1);
					double d20 = DDG::dot(v2, v0);
					double d21 = DDG::dot(v2, v1);

					double denom = d00 * d11 - d01 * d01;
					double v = (d11 * d20 - d01 * d21) / denom;
					double w = (d00 * d21 - d01 * d20) / denom;
					double u = 1.0f - v - w;

					if (v >= 0 && w >= 0 && u <= 1) {
						DDG::Vector P = u * P1 + v * P2 + w * P3;
						angle = DDG::angle(P2 - P1, P - P1) + vv_idx_angle[i].second;
						findit = true;
					}
				}
			}
		}
		// loop
		else {
			for (int i = 0; i < vv_idx_angle.size(); i++) {
				DDG::Vector P1, P2, P3;
				if (i == vv_idx_angle.size() - 1) {
					P1 = vit.position;
					P2 = mesh.vertices[vv_idx_angle[i].first].position;
					P3 = mesh.vertices[vv_idx_angle[0].first].position;
				}
				else {
					P1 = vit.position;
					P2 = mesh.vertices[vv_idx_angle[i].first].position;
					P3 = mesh.vertices[vv_idx_angle[i + 1].first].position;
				}

				DDG::Vector v0 = P2 - P1;
				DDG::Vector v1 = P3 - P1;
				DDG::Vector v2 = initFieldVector;

				double d00 = DDG::dot(v0, v0);
				double d01 = DDG::dot(v0, v1);
				double d11 = DDG::dot(v1, v1);
				double d20 = DDG::dot(v2, v0);
				double d21 = DDG::dot(v2, v1);

				double denom = d00 * d11 - d01 * d01;
				double v = (d11 * d20 - d01 * d21) / denom;
				double w = (d00 * d21 - d01 * d20) / denom;
				double u = 1.0f - v - w;

				if (v >= -1e-3 && w >= -1e-3 && u <= 1 + 1e-3) {
					DDG::Vector P = u * P1 + v * P2 + w * P3;
					angle = DDG::angle(P2 - P1, P - P1) + vv_idx_angle[i].second;
					findit = true;
				}

			}
		}

		if (!(findit || vit.onBoundary())) {
			int jiji = 0;
			//std::cout << "can not find vertex index v_idx " << vit.index << std::endl;
		}
		double phi = 2 * angle + M_PI;
		vit.directionField = DDG::Complex(cos(phi), sin(phi));
	}

}
