#ifndef GROUP_H
#define GROUP_H


#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <iostream>
#include <vector>


// TODO: Implement Group - add data structure to store a list of Object*
class Group : public Object3D {

public:

    Group();

    explicit Group (int num_objects);

    ~Group() override;

    bool intersect(const Ray &r, Hit &h, float tmin) override;

    void addObject(int index, Object3D *obj);

    int getGroupSize();

    void print() override {
        for (Object3D *o : *container_) {
            o->print();
        }
    }

private:
    std::vector<Object3D *> *container_;
};

#endif
	
