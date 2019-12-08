#include "Scene.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Camera.h"
#include "Color.h"
#include "Helpers.h"
#include "Matrix4.h"
#include "Model.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"
#include "tinyxml2.h"

using namespace tinyxml2;
using namespace std;

void copy(double dest[4][4], double src[4][4]) {
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) dest[i][j] = src[i][j];
}

Matrix4 getTranslationTransformation(Translation *translation) {
	double mat[4][4] = {{1.0, .0, .0, translation->tx},
	                    {.0, 1.0, .0, translation->ty},
	                    {.0, .0, 1.0, translation->tz},
	                    {.0, .0, .0, 1.0}};

	return Matrix4(mat);
}

Matrix4 getScalingTransformation(Scaling *scaling) {
	double mat[4][4] = {{scaling->sx, .0, .0, .0},
	                    {.0, scaling->sy, .0, .0},
	                    {.0, .0, scaling->sz, .0},
	                    {.0, .0, .0, 1.0}};

	return Matrix4(mat);
}

Matrix4 getRotationTransformation(Rotation *rotation) {
	Vec3 u(rotation->ux, rotation->uy, rotation->uz, -1);
	u = normalizeVec3(u);

	char min_component = 'x';
	if (abs(u.y) < abs(u.x)) min_component = 'y';
	if (abs(u.z) < abs(u.y) && abs(u.z) < abs(u.x)) min_component = 'z';

	Vec3 v;
	switch (min_component) {
		case 'x':
			v = Vec3(.0, -u.z, u.y, -1);
			break;
		case 'y':
			v = Vec3(-u.z, .0, u.x, -1);
			break;
		case 'z':
			v = Vec3(-u.y, u.x, .0, -1);
			break;
	}
	v = normalizeVec3(v);

	Vec3 w = crossProductVec3(u, v);

	double mat0[4][4] = {{u.x, u.y, u.z, .0},
	                     {v.x, v.y, v.z, .0},
	                     {w.x, w.y, w.z, .0},
	                     {.0, .0, .0, 1.0}};
	Matrix4 M(mat0);

	double mat1[4][4] = {{u.x, v.x, w.x, .0},
	                     {u.y, v.y, w.y, .0},
	                     {u.z, v.z, w.z, .0},
	                     {.0, .0, .0, 1.0}};
	Matrix4 M_inverse(mat1);

	double theta = (rotation->angle * M_PI) / 180.0;
	double mat2[4][4] = {{1.0, .0, .0, .0},
	                     {.0, cos(theta), -sin(theta), .0},
	                     {.0, sin(theta), cos(theta), .0},
	                     {.0, .0, .0, 1.0}};
	Matrix4 Rx_theta(mat2);

	return multiplyMatrixWithMatrix(M_inverse,
	                                multiplyMatrixWithMatrix(Rx_theta, M));
}

Matrix4 Scene::getModelTransformation(int id, char type) {
	switch (type) {
		case 't':
			return getTranslationTransformation(translations[id]);
		case 's':
			return getScalingTransformation(scalings[id]);
		case 'r':
			return getRotationTransformation(rotations[id]);
	}

	return Matrix4();  // ?
}

Matrix4 getCameraTransformation(Camera *c) {
	double mat[4][4] = {{c->u.x, c->u.y, c->u.z, -dotProductVec3(c->u, c->pos)},
	                    {c->v.x, c->v.y, c->v.z, -dotProductVec3(c->v, c->pos)},
	                    {c->w.x, c->w.y, c->w.z, -dotProductVec3(c->w, c->pos)},
	                    {.0, .0, .0, 1.0}};

	return Matrix4(mat);
}

Matrix4 getOrthographicTransformation(Camera *c) {
	double n = c->near, f = c->far, l = c->left, r = c->right, b = c->bottom,
	       t = c->top;
	double mat[4][4] = {{2.0 / (r - l), .0, .0, -(r + l) / (r - l)},
	                    {.0, 2.0 / (t - b), .0, -(t + b) / (t - b)},
	                    {.0, .0, -2.0 / (f - n), -(f + n) / (f - n)},
	                    {.0, .0, .0, 1.0}};

	return Matrix4(mat);
}

Matrix4 getPerspectiveTransformation(Camera *c) {
	double n = c->near, f = c->far, l = c->left, r = c->right, b = c->bottom,
	       t = c->top;
	double mat[4][4] = {{2.0 * n / (r - l), .0, (r + l) / (r - l), .0},
	                    {.0, 2.0 * n / (t - b), (t + b) / (t - b), .0},
	                    {.0, .0, -(f + n) / (f - n), -2.0 * f * n / (f - n)},
	                    {.0, .0, -1.0, .0}};

	return Matrix4(mat);
}

Matrix4 getProjectionTransformation(Camera *c, int projectionType) {
	return (projectionType == 0) ? getOrthographicTransformation(c)
	                             : getPerspectiveTransformation(c);
}

Matrix4 getViewportTransformation(Camera *c) {
	double mat[4][4] = {{c->horRes / 2.0, .0, .0, (c->horRes - 1.0) / 2.0},
	                    {.0, c->verRes / 2.0, .0, (c->verRes - 1.0) / 2.0},
	                    {.0, .0, .5, .5},
	                    {0, 0, 0, 0}};

	return Matrix4(mat);
}

bool Scene::backFaceCulling(Vec3 *v0, Vec3 *v1, Vec3 *v2, Vec3 v) {
	Vec3 norm =
	    crossProductVec3(subtractVec3(*v1, *v0), subtractVec3(*v2, *v0));
	Vec3 eye_to_point = subtractVec3(*v0, v);

	if (dotProductVec3(norm, eye_to_point) > 0)
		return true;
	else
		return false;
}

// bool visible(double num, double den, double *t_e, double *t_l) {
// 	if (den > .0) {
// 		double t = num / den;
// 		if (t > *t_l) return false;
// 		if (t > *t_e) *t_e = t;
// 	} else if (den < .0) {
// 		double t = num / den;
// 		if (t < *t_e) return false;
// 		if (t < *t_l) *t_l = t;
// 	} else if (num > .0)
// 		return false;

// 	return true;
// }

// void clipping(Camera *camera, Vec3 v0, Vec3 v1) {
// 	double t_e = .0, t_l = 1.0;
// 	double dx = v1.x - v0.x, dy = v1.y - v0.y, dz = v1.z - v0.z;
// 	// bool isVisible = false;

// 	if (visible(dx, camera->left - v0.x, *t_e, *t_l))
// 		if (visible(-dx, v0.x - camera->right, *t_e, *t_l))
// 			if (visible(dy, camera->bottom - v0.y, *t_e, *t_l))
// 				if (visible(-dy, v0.y - camera->top, *t_e, *t_l))
// 					if (visible(dz, camera->near - v0.z, *t_e, *t_l))
// 						if (visible(-dz, v0.z - camera->far, *t_e, *t_l)) {
// 							// isVisible = true;
// 							if (t_l < 1.0) {
// 								Vec3 dt = multiplyVec3WithScalar(
// 								    Vec3(dx, dy, dz, -1), t_l);
// 								v1 = addVec3(v0, dt);
// 							} else if (t_e > .0) {
// 								Vec3 dt = multiplyVec3WithScalar(
// 								    Vec3(dx, dy, dz, -1), t_e);
// 								v0 = addVec3(v0, dt);
// 							}
// 						}
// }

Model *Scene::doTransformations(Camera *camera, Model *model) {
	// Transformations
	Model *model_t = new Model();
	model_t->modelId = model->modelId;
	model_t->type = model->type;
	model_t->numberOfTriangles = model->numberOfTriangles;

	Matrix4 transformations = getIdentityMatrix();
	for (int i = 0; i < model->transformationIds.size(); i++) {
		int id = model->transformationIds[i] - 1;
		char type = model->transformationTypes[i];

		Matrix4 transformation = getModelTransformation(id, type);

		Matrix4 tmp = multiplyMatrixWithMatrix(transformation, transformations);
		copy(transformations.val, tmp.val);
	}

	Matrix4 camera_t = getCameraTransformation(camera);
	Matrix4 projection_t = getProjectionTransformation(camera, projectionType);
	Matrix4 viewport_t = getViewportTransformation(camera);

	Matrix4 tmp = multiplyMatrixWithMatrix(camera_t, transformations);
	transformations = multiplyMatrixWithMatrix(projection_t, tmp);

	for (int i = 0; i < model->numberOfTriangles; i++) {
		int *vertexIds = model->triangles[i].vertexIds;
		// vector< Triangle* > triangles_t;
		Triangle triangle_t;
		vector<Vec3 *> vertices_tmp;
		vector<vector<double> > mat(3);

		for (int j = 0; j < 3; j++) {
			Vec4 tmp(vertices[*(vertexIds + j) - 1]->x,
			         vertices[*(vertexIds + j) - 1]->y,
			         vertices[*(vertexIds + j) - 1]->z, 1.0,
			         vertices[*(vertexIds + j) - 1]->colorId);

			tmp = multiplyMatrixWithVec4(transformations, tmp);

			double _mat[4] = {tmp.x, tmp.y, tmp.z, tmp.t};
			mat[j].assign(_mat, _mat + 4);

			Vec3 *vertex_t = new Vec3();
			vertex_t->x = mat[j][0];
			vertex_t->y = mat[j][1];
			vertex_t->z = mat[j][2];
			vertex_t->colorId = tmp.colorId - 1;

			vertices_tmp.push_back(vertex_t);
		}

		bool isVisible = true;
		if (cullingEnabled)
			isVisible = backFaceCulling(vertices_tmp[0], vertices_tmp[1],
			                            vertices_tmp[2], camera->v);

		// CLIPPING STAGE
		// vector< Vec3* > polygon;
		// Vec3 _v0(vertices_tmp[0]);
		// Vec3 _v1(vertices_tmp[1]);
		// clipping(camera, _v0, _v1);
		// polygon.push_back(_v0, _v1);

		// Vec3 __v1(vertices_tmp[1]);
		// Vec3 __v2(vertices_tmp[2]);
		// clipping(camera, __v1, __v2);

		// Vec3 ___v2(vertices_tmp[2]);
		// Vec3 ___v0(vertices_tmp[0]);
		// clipping(camera, ___v2, ___v0);

		for (int j = 0; j < 3 && isVisible; j++) {
			for (int _i = 0; _i < 3; _i++) {
				double total = 0;
				for (int _j = 0; _j < 4; _j++)
					total += viewport_t.val[_i][_j] * (mat[j][_j] / mat[j][3]);

				switch (_i) {
					case 0:
						vertices_tmp[j]->x = total;
						break;
					case 1:
						vertices_tmp[j]->y = total;
						break;
					case 2:
						vertices_tmp[j]->z = total;
						break;
				}
			}

			vertices_t.push_back(vertices_tmp[j]);
			triangle_t.vertexIds[j] = vertices_t.size() - 1;
		}

		if (isVisible)
			model_t->triangles.push_back(triangle_t);
		else
			model_t->numberOfTriangles--;
	}

	// models_t.push_back(model_t);
	return model_t;
}

void Scene::midpoint(int x0, int y0, int x1, int y1, Color c0, Color c1) {
	// Bresenham
	const bool steep = abs(y1 - y0) > abs(x1 - x0);
	if (steep) {
		swap(x0, y0);
		swap(x1, y1);
	}

	if (x0 > x1) {
		swap(x0, x1);
		swap(y0, y1);
		swap(c0, c1);
	}

	const int dx = x1 - x0, dy = abs(y1 - y0);

	float error = dx / 2.0f;
	const int ystep = (y0 < y1) ? 1 : -1;
	int y = y0;
	Color c = c0, dc;
	dc.r = (c1.r - c0.r) / dx;
	dc.g = (c1.g - c0.g) / dx;
	dc.b = (c1.b - c0.b) / dx;

	for (int x = x0; x <= x1; x++) {
		if (steep)
			image[y][x] = Color(c);
		else
			image[x][y] = Color(c);

		error -= dy;
		if (error < 0) {
			y += ystep;
			error += dx;
		}

		c.r += dc.r;
		c.g += dc.g;
		c.b += dc.b;
	}
}

void Scene::midpoint(Model *model) {
	if (model->type == 1) return;

	for (int i = 0; i < model->numberOfTriangles; i++) {
		int *vertexIds = model->triangles[i].vertexIds;

		for (int j = 0; j < 3; j++) {
			int x0 = vertices_t[*(vertexIds + j)]->x,
			    y0 = vertices_t[*(vertexIds + j)]->y,
			    x1 = vertices_t[*(vertexIds + ((j + 1) % 3))]->x,
			    y1 = vertices_t[*(vertexIds + ((j + 1) % 3))]->y;
			Color c0 = *colorsOfVertices[vertices_t[*(vertexIds + j)]->colorId],
			      c1 =
			          *colorsOfVertices[vertices_t[*(vertexIds + ((j + 1) % 3))]
			                                ->colorId];

			midpoint(x0, y0, x1, y1, c0, c1);
		}
	}
}

double f01(int x, int y, int x0, int y0, int x1, int y1) {
	double f = x * (y0 - y1) - y * (x0 - x1) + x0 * y1 - y0 * x1;

	return f;
}

double f12(int x, int y, int x1, int y1, int x2, int y2) {
	double f = x * (y1 - y2) - y * (x1 - x2) + x1 * y2 - y1 * x2;

	return f;
}

double f20(int x, int y, int x2, int y2, int x0, int y0) {
	double f = x * (y2 - y0) - y * (x2 - x0) + x2 * y0 - y2 * x0;

	return f;
}

void Scene::triangleRasterization(Model *model) {
	if (model->type == 0) return;

	for (int i = 0; i < model->numberOfTriangles; i++) {
		int *vertexIds = model->triangles[i].vertexIds;

		int x0 = vertices_t[*(vertexIds)]->x, y0 = vertices_t[*(vertexIds)]->y,
		    x1 = vertices_t[*(vertexIds + 1)]->x,
		    y1 = vertices_t[*(vertexIds + 1)]->y,
		    x2 = vertices_t[*(vertexIds + 2)]->x,
		    y2 = vertices_t[*(vertexIds + 2)]->y;
		Color *c0 = colorsOfVertices[vertices_t[*(vertexIds)]->colorId],
		      *c1 = colorsOfVertices[vertices_t[*(vertexIds + 1)]->colorId],
		      *c2 = colorsOfVertices[vertices_t[*(vertexIds + 2)]->colorId];

		int x_min = min(min(x0, x1), x2), x_max = max(max(x0, x1), x2);
		int y_min = min(min(y0, y1), y2), y_max = max(max(y0, y1), y2);

		for (int y = y_min; y <= y_max; y++) {
			if (y >= ((int)image[0].size() - 1) || y < 0)
				continue;  // remove after clipping impl
			for (int x = x_min; x <= x_max; x++) {
				if (x >= ((int)image.size() - 1) || x < 0)
					continue;  // remove after clipping impl
				double alpha = f12(x, y, x1, y1, x2, y2) /
				               f12(x0, y0, x1, y1, x2, y2),
				       beta = f20(x, y, x2, y2, x0, y0) /
				              f20(x1, y1, x2, y2, x0, y0),
				       gamma = f01(x, y, x0, y0, x1, y1) /
				               f01(x2, y2, x0, y0, x1, y1);

				if (alpha >= 0 && beta >= 0 && gamma >= 0) {
					image[x][y].r =
					    alpha * c0->r + beta * c1->r + gamma * c2->r;
					image[x][y].g =
					    alpha * c0->g + beta * c1->g + gamma * c2->g;
					image[x][y].b =
					    alpha * c0->b + beta * c1->b + gamma * c2->b;
				}
			}
		}
	}
}

/*
        Transformations, clipping, culling, rasterization are done here.
        You can define helper functions inside Scene class implementation.
*/
void Scene::forwardRenderingPipeline(Camera *camera) {
	// TODO
	for (Model *model : models) {
		Model *model_t = doTransformations(camera, model);
		midpoint(model_t);
		triangleRasterization(model_t);

		vertices_t.clear();
	}
}

/*
        Parses XML file
*/
Scene::Scene(const char *xmlPath) {
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g,
	       &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL) pElement->QueryBoolText(&cullingEnabled);

	// read projection type
	pElement = pRoot->FirstChildElement("ProjectionType");
	if (pElement != NULL) pElement->QueryIntText(&projectionType);

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL) {
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d", &cam->left, &cam->right,
		       &cam->bottom, &cam->top, &cam->near, &cam->far, &cam->horRes,
		       &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL) {
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL) {
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty,
		       &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL) {
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL) {
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux,
		       &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read models
	pElement = pRoot->FirstChildElement("Models");

	XMLElement *pModel = pElement->FirstChildElement("Model");
	XMLElement *modelElement;
	while (pModel != NULL) {
		Model *model = new Model();

		pModel->QueryIntAttribute("id", &model->modelId);
		pModel->QueryIntAttribute("type", &model->type);

		// read model transformations
		XMLElement *pTransformations =
		    pModel->FirstChildElement("Transformations");
		XMLElement *pTransformation =
		    pTransformations->FirstChildElement("Transformation");

		pTransformations->QueryIntAttribute("count",
		                                    &model->numberOfTransformations);

		while (pTransformation != NULL) {
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			model->transformationTypes.push_back(transformationType);
			model->transformationIds.push_back(transformationId);

			pTransformation =
			    pTransformation->NextSiblingElement("Transformation");
		}

		// read model triangles
		XMLElement *pTriangles = pModel->FirstChildElement("Triangles");
		XMLElement *pTriangle = pTriangles->FirstChildElement("Triangle");

		pTriangles->QueryIntAttribute("count", &model->numberOfTriangles);

		while (pTriangle != NULL) {
			int v1, v2, v3;

			str = pTriangle->GetText();
			sscanf(str, "%d %d %d", &v1, &v2, &v3);

			model->triangles.push_back(Triangle(v1, v2, v3));

			pTriangle = pTriangle->NextSiblingElement("Triangle");
		}

		models.push_back(model);

		pModel = pModel->NextSiblingElement("Model");
	}
}

/*
        Initializes image with background color
*/
void Scene::initializeImage(Camera *camera) {
	if (this->image.empty()) {
		for (int i = 0; i < camera->horRes; i++) {
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++) {
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	// if image is filled before, just change color rgb values with the
	// background color
	else {
		for (int i = 0; i < camera->horRes; i++) {
			for (int j = 0; j < camera->verRes; j++) {
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
        If given value is less than 0, converts value to 0.
        If given value is more than 255, converts value to 255.
        Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value) {
	if (value >= 255.0) return 255;
	if (value <= 0.0) return 0;
	return (int)(value);
}

/*
        Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera) {
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--) {
		for (int i = 0; i < camera->horRes; i++) {
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
			     << makeBetweenZeroAnd255(this->image[i][j].g) << " "
			     << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
        Converts PPM image in given path to PNG file, by calling ImageMagick's
   'convert' command. os_type == 1 		-> Ubuntu os_type == 2 ->
   Windows os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType) {
	string command;

	// call command on Ubuntu
	if (osType == 1) {
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2) {
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else {
	}
}