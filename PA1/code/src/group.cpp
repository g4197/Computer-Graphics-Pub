#include "group.hpp"
Group::Group() {
    container_ = new std::vector<Object3D *>;
}

Group::Group (int num_objects) {
    container_ = new std::vector<Object3D *>(num_objects);
}

Group::~Group() {
    delete container_;
}

bool Group::intersect(const Ray &r, Hit &h, float tmin) {
    bool ret = false;
    for (int i = 0; i < container_->size(); ++i) {
        ret |= (*container_)[i]->intersect(r, h, tmin);
    }
    return ret;
}

void Group::addObject(int index, Object3D *obj) {
    assert(container_->size() > index);
    (*container_)[index] = obj;
}

int Group::getGroupSize() {
    return container_->size();
}