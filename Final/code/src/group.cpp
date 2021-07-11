#include "group.hpp"
Group::Group() : tree_(nullptr) {
    container_ = new std::vector<Object3D *>;
}

Group::Group (int num_objects) : tree_(nullptr) {
    container_ = new std::vector<Object3D *>(num_objects);
}

Group::~Group() {
    delete container_;
}

Object3D *Group::operator[](int idx) const {
    return (*container_)[idx];
}

bool Group::intersect(const Ray &r, Hit &h, double tmin) {
    // if (!tree_) {
    //    bool ret = false;
    //    for (int i = 0; i < container_->size(); ++i) {
    //        ret |= (*container_)[i]->intersect(r, h, tmin);
    //    }
    //    return ret;
    // } else {
    bool ret = tree_->intersect(r, h, tmin);
    // }
    return ret;
}

void Group::addObject(int index, Object3D *obj) {
    assert(container_->size() > index);
    (*container_)[index] = obj;
}

void Group::buildTree() {
    tree_ = new KDTree(*container_);
}

int Group::size() {
    return container_->size();
}