#ifndef PPM_H_
#define PPM_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <list>
#include "omp.h"
#include "vecmath.h"

#include "scene_parser.hpp"
#include "image.hpp"
#include "utils.hpp"
#include "camera.hpp"
#include "group.hpp"
#include "sphere.hpp"
#include "material.hpp"

Image photonMap(SceneParser &parser, int samples);

#endif // PPM_H_