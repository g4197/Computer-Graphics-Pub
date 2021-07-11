#ifndef GROUP_H
#define GROUP_H


#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <iostream>
#include <vector>


// TODO (PA2): Implement Group - copy from PA1
class Group : public Object3D {

public:

    Group() {
        container_ = new std::vector<Object3D *>;
    }

    explicit Group (int num_objects) {
        container_ = new std::vector<Object3D *>(num_objects);
    }

    ~Group() override {
        delete container_;
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
		bool ret = false;
        for (int i = 0; i < container_->size(); ++i) {
            ret |= (*container_)[i]->intersect(r, h, tmin);
        }
        return ret;
    }

    void drawGL() override {
        for (size_t i = 0; i < container_->size(); ++i) {
            (*container_)[i]->drawGL();
        }
    }

    void addObject(int index, Object3D *obj) {
        assert(container_->size() > index);
        (*container_)[index] = obj;
    }

    int getGroupSize() {
        return container_->size();
    }

private:
    std::vector<Object3D *> *container_;
};

#endif
	
