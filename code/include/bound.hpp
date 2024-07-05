//参考已有代码：https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project
#ifndef BOUND_H
#define BOUND_H 

#include <vecmath.h>
#include <vector>
#include "object3d.hpp"
#include <cmath>
#include "constants.hpp"
using namespace std;

class Bound{
public:
    void cmin(Vector3f& a, const Vector3f& b){
        if(a.x() > b.x()) a.x() = b.x();
        if(a.y() > b.y()) a.y() = b.y();
        if(a.z() > b.z()) a.z() = b.z();
    }

    void cmax(Vector3f& a, const Vector3f& b){
        if(a.x() < b.x()) a.x() = b.x();
        if(a.y() < b.y()) a.y() = b.y();
        if(a.z() < b.z()) a.z() = b.z();
    }

    Bound(){
        bound[0] = Vector3f(-INF);
        bound[1] = Vector3f(INF);
    }

    Bound(const vector<Vector3f> v, Object3D *obj){
        this->o = obj;
        bound[0] = v[0];
        bound[1] = v[0];
        for(auto i : v){
            cmin(bound[0], i);
            cmax(bound[1], i);
        }
    }

    void set(const Vector3f &lo, const Vector3f &hi){
        bound[0] = lo;
        bound[1] = hi;
    }

    bool intersect(const Ray &r, Hit &h, double tmin){
        Vector3f o(r.getOrigin());
        Vector3f dir(r.getDirection());
        Vector3f invdir(1 / dir);
        vector<int> sgn = { invdir.x() < 0, invdir.y() < 0, invdir.z() < 0 };
        tmin = INF;
        float t_min, t_max, tymin, tymax, tzmin, tzmax;
        t_min = (bound[sgn[0]].x() - o.x()) * invdir.x();
        t_max = (bound[1 - sgn[0]].x() - o.x()) * invdir.x();
        tymin = (bound[sgn[1]].y() - o.y()) * invdir.y();
        tymax = (bound[1 - sgn[1]].y() - o.y()) * invdir.y();
        tzmin = (bound[sgn[2]].z() - o.z()) * invdir.z();
        tzmax = (bound[1 - sgn[2]].z() - o.z()) * invdir.z();
        if ((t_min > tymax) || (tymin > t_max)) return false;
        if (tymin > t_min) t_min = tymin;
        if (tymax < t_max) t_max = tymax;
        if ((t_min > tzmax) || (tzmin > t_max)) return false;
        if (tzmin > t_min) t_min = tzmin;
        if (tzmax < t_max) t_max = tzmax;
        tmin = t_min;
        return true;
    }

protected:
    Object3D *o; //un-transformed object
    Vector3f bound[2];
};


#endif