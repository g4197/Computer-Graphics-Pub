#ifndef PT_H_
#define PT_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <omp.h>
#include "vecmath.h"

#include "utils.hpp"
#include "scene_parser.hpp"
#include "image.hpp"
#include "utils.hpp"
#include "camera.hpp"
#include "group.hpp"
#include "sphere.hpp"
#include "material.hpp"

Image pathTrace(SceneParser &parser, int samples);

#endif //PT_H_