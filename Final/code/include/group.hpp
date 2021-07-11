#ifndef GROUP_H
#define GROUP_H


#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include "kdtree.hpp"
#include <iostream>
#include <vector>

class Group : public Object3D {

public:

    Group();

    explicit Group (int num_objects);

    ~Group() override;

    Object3D *operator[](int idx) const;

    bool intersect(const Ray &r, Hit &h, double tmin) override;

    void addObject(int index, Object3D *obj);

    void buildTree();

    int size();

    void print() override {
        for (Object3D *o : *container_) {
            o->print();
        }
    }

private:
    std::vector<Object3D *> *container_;
    KDTree *tree_;
};

#endif
	
