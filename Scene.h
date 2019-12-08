#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Matrix4.h"
#include "Model.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"

using namespace std;

class Scene {
   public:
	Color backgroundColor;
	bool cullingEnabled;
	int projectionType;

	vector<vector<Color> > image;
	vector<Camera*> cameras;
	vector<Vec3*> vertices;
	vector<Color*> colorsOfVertices;
	vector<Scaling*> scalings;
	vector<Rotation*> rotations;
	vector<Translation*> translations;
	vector<Model*> models;

	Scene(const char* xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

   private:
	vector<Vec3*> vertices_t;
	Matrix4 getModelTransformation(int id, char type);
	bool backFaceCulling(Vec3* v0, Vec3* v1, Vec3* v2, Vec3 v);
	Model* doTransformations(Camera* camera, Model* model);
	void midpoint(int x0, int y0, int x1, int y1, Color c0, Color c1);
	void midpoint(Model* model);
	void triangleRasterization(Model* model);
};

#endif
