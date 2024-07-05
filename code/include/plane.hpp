//PA1框架，setHit函数参考已有代码https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project，其余独立实现
#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement Plane representing an infinite plane
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D {
public:
    Plane(): P(Vector3f::UP), d(0) {

    }

    Plane(const Vector3f &normal, float d, Material *m) : Object3D(m), P(normal), d(d) {
        assert(normal.length() == 1);
        N = Vector3f::cross(Vector3f::UP, P);
    }

    ~Plane() override = default;

    void setHit(Hit &h, float t, Vector3f normal, Vector3f dir, Vector3f next_origin, Vector3f N, Material *material){
        Vector2f grad = Vector2f::ZERO;
        Vector3f new_normal = normal;
        float u = Vector3f::dot(P - dir * normal, N);
        float v = next_origin.y();
        float f = material->getBump().getDisturb(u, v, grad);
        if (f > 1e-4 || f < -1e-4) {  //扰动梯度求法向量
            if (N.squaredLength() > 1e-4 || N.squaredLength() < -1e-4){
                new_normal = Vector3f::cross(N + normal * grad[0], Vector3f::UP + normal * grad[1]);
                new_normal.normalized();
            }
        }
        h.set(t, material, new_normal, material->getColor(u, v), next_origin, material->getTexture().isTexture);
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        Vector3f origin(r.getOrigin());
        Vector3f dir(r.getDirection());
        float x = Vector3f::dot(P, dir);
        if(fabs(x) < 1e-6) return false; //光线与平面平行
        float t = (d - Vector3f::dot(P, origin)) / Vector3f::dot(P, dir);
        if (t < tmin || t > h.getT()) return false;

        //梯度扰动
        Vector3f next_origin(origin + dir * t);
        setHit(h, t, P, dir, next_origin, N, material);
        // h.set(t, material, P);
        return true;
    }

protected:
    Vector3f P; //法向量
    Vector3f N; 
    float d;

};

#endif //PLANE_H
		

