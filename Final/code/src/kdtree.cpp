#include "kdtree.hpp"
#include "triangle.hpp"

KDTNode::KDTNode(Vector3f l, Vector3f r, Object3D *lc, Object3D *rc) :
    Object3D(nullptr, l, r), lc_(lc), rc_(rc) {}

KDTNode::~KDTNode() {
    // only delete KDTNodes
    if (dynamic_cast<KDTNode *>(lc_) != nullptr) {
        delete lc_;
    }
    if (dynamic_cast<KDTNode *>(rc_) != nullptr) {
        delete rc_;
    }
}

bool KDTNode::intersect(const Ray &r, Hit &h, double tmin) {
    if (!aabb_.intersect(r)) {
        return false;
    } else {
        bool ret = false;
        ret |= lc_->intersect(r, h, tmin);
        ret |= rc_->intersect(r, h, tmin);
        return ret;
    }
}

KDTree::KDTree(std::vector<Object3D *> objects) {
    this->root_ = build(0, 0, objects.size() - 1, objects);
}

Object3D *KDTree::build(int layer, int l, int r, std::vector<Object3D *> &objects) {
    if (r == l) {
        return objects[l];
    }
    int m = (l + r) / 2, layer_mod_3 = layer % 3;
    nth_element(
        objects.begin() + l, objects.begin() + m, objects.begin() + r, 
        [layer_mod_3] (Object3D *x, Object3D *y) {
        return x->aabb().l()[layer_mod_3] < y->aabb().l()[layer_mod_3];
    });
    Object3D *lc = build(layer + 1, l, m, objects);
    Object3D *rc = build(layer + 1, m + 1, r, objects);
    return new KDTNode(minVec(lc->aabb().l(), rc->aabb().l()), 
                       maxVec(lc->aabb().r(), rc->aabb().r()), lc, rc);
}

KDTree::~KDTree() {
    delete root_;
}

bool KDTree::intersect(const Ray &r, Hit &h, double tmin) {
    return root_->intersect(r, h, tmin);
}
