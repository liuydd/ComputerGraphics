//PA1框架
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

    Group() {

    }

    explicit Group (int num_objects): object_list(num_objects) {

    }

    ~Group() override {

    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        bool flag = false;
        for(auto obj: object_list){
            if(obj){
                flag = obj->intersect(r, h, tmin) || flag;
            }
        }
        return flag;
    }

    void addObject(int index, Object3D *obj) {
        object_list.insert(object_list.begin() + index, obj);
    }

    int getGroupSize() {
        return object_list.size();
    }

private:
    std::vector<Object3D*> object_list;
};

#endif
	
