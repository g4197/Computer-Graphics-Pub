#ifndef MESH_H
#define MESH_H

#include <vector>
#include "object3d.hpp"
#include "triangle.hpp"
#include "Vector2f.h"
#include "Vector3f.h"
#include "kdtree.hpp"
#include "tiny_obj_loader.hpp"


class Mesh : public Object3D {

public:
    Mesh(const char *filename, Material *m);
    bool intersect(const Ray &r, Hit &h, double tmin) override;
    void print() override {
        printf("A mesh with %lu triangles\n", triangle_vec_.size());
    }
private:
    void buildTree();
    KDTree *tree_;
    vector<Object3D *> triangle_vec_;
};

#endif
