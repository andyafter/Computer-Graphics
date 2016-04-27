#pragma once

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include "glm/glm.hpp"
#include "utilities.h"
#include "sceneStructs.h"
#include "image.h"

using namespace std;

class Scene {
private:
    ifstream fp_in;
    int loadMaterial(string materialid);
    int loadGeom(string objectid);
    int loadCamera();
public:
    Scene(string filename);
    ~Scene();

	int lightIdx = 0;//??? multilight
	std::vector<int> lightIdxs;
    std::vector<Geom> geoms;
    std::vector<Material> materials;
	std::vector<image> textures;
    RenderState state;
};