#ifndef KDTREE_H_
#define KDTREE_H_

#include <algorithm>
#include <iostream>
#include "object3d.hpp"
using namespace std;

class KDTNode : public Object3D {
public:
    KDTNode(Vector3f l, Vector3f r, Object3D *lc, Object3D *rc);
    ~KDTNode();
    bool intersect(const Ray &r, Hit &h, double tmin) override;
private:
    Object3D *lc_, *rc_;
};

class KDTree {
public:
    KDTree(std::vector<Object3D *> objects);
    ~KDTree();
    Object3D *build(int layer, int l, int r, std::vector<Object3D *> &objects);
    bool intersect(const Ray &r, Hit &h, double tmin);
private:
    Object3D *root_;
};

#endif // KDTREE_H_