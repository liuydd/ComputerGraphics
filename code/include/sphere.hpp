#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement functions and add more fields as necessary

class Sphere : public Object3D {
public:
    Sphere(): radius(0), center(0, 0, 0) {
        // unit ball at the center
    }

    Sphere(const Vector3f &center, float radius, Material *material, const Vector3f &velocity) : Object3D(material), center(center), radius(radius), velocity(velocity) {
        if(velocity != Vector3f::ZERO) isMove = true;
    }

    ~Sphere() override = default;

    //参考已有代码https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project
    void setHit(Hit &h, float t, Vector3f N, Vector3f P, Material *material){
        // 更新 Hit 对象
        float u = 0.5 + atan2(N.x(), N.z()) / (2 * M_PI);
        float v = 0.5 - asin(N.y()) / M_PI;
        Vector2f grad(0, 0); // 存储凹凸贴图的梯度信息
        float f = material->getBump().getDisturb(u, v, grad);
        if (f > 1e-4 || f < -1e-4) {  //扰动梯度求法向量
            float phi = u * 2 * M_PI;
            float theta = M_PI - v * M_PI;
            Vector3f pu(-N.z(), 0, N.x()), pv(N.y() * cos(phi), -radius * sin(theta), N.y() * sin(phi));
            N = Vector3f::cross(pu + N * grad[0] / (2 * M_PI), pv + N * grad[1] / M_PI).normalized();
        }
        h.set(t, material, N, material->getColor(u, v), P, material->getTexture().isTexture);
    }

    //除对运动模糊处理是参考已有代码https://github.com/zrporz/THU-CG-2023外其余独立完成
    bool intersect(const Ray &r, Hit &h, float tmin) override {
        //
        // Vector3f origin(r.getOrigin());
        Vector3f origin;
        if(isMove) origin = r.getOrigin() - RND2 * velocity; //运动造成中心点的随机偏移
        else origin = r.getOrigin();
        Vector3f dir(r.getDirection());
        Vector3f l(center - origin); //由光源O指向球心C的向量
        dir.normalize();
        float OH = Vector3f::dot(l, dir); //OC在射线方向的投影
        float d = sqrt(fabs(l.squaredLength() - OH * OH));
        if (d > radius) return false;
        float tt = sqrt(fabs(radius * radius - d * d));
        float t = OH - tt;
        if (t < tmin || t > h.getT()) return false;
        Vector3f P(origin + dir * t); //与球面的相交点
        Vector3f N((P - center).normalized());

        setHit(h, t, N, P, material); //实现纹理贴图
        // h.set(t, material, N);
        return true;
    }

protected:
    Vector3f center;
    float radius;
    Vector3f velocity;
    bool isMove = false;
};


#endif
