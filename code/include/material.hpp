#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "hit.hpp"
#include "texture.h"
#include "utils.hpp"
#include <iostream>

// TODO: Implement Shade function that computes Phong introduced in class.
class Material {
public:
    //函数部分独立实现，Shade函数参考PA1 readme公式
    explicit Material(const Vector3f &d_color, const Vector3f &s_color = Vector3f::ZERO, float s = 0, const Vector3f &color = Vector3f::ZERO, const Vector3f &e_color = Vector3f::ZERO, float r = 0, const Vector4f &t = Vector4f(1, 0, 0, 0), const char *textureFile = "", const char *bumpFile = "") :
            diffuseColor(d_color), specularColor(s_color), shininess(s), color(color), emission(e_color), refr(r), texture(textureFile), bump(bumpFile) {
                this->type = t / (t.x() + t.y() + t.z() + t.w());
    }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }

    Vector3f getEmissionColor() const {
        return emission;
    }

    Vector4f getType() const{
        return type;
    }

    float getRefr() const{
        return refr;
    }

    Texture getBump() const {
        return bump;
    }

    Texture getTexture() const {
        return texture;
    }

    Vector3f getColor(float u, float v) const{
        if(!texture.isTexture) return this->color;
        else return texture.getColor(u, v);
    }

    float getDisturb(float u, float v, Vector2f& grad) {
        return bump.getDisturb(u, v, grad);
    }


    Vector3f Shade(const Ray &ray, const Hit &hit,
                   const Vector3f &dirToLight, const Vector3f &lightColor) {
        Vector3f shaded = Vector3f::ZERO;
        // 
        Vector3f N = hit.getNormal().normalized();
        Vector3f V = -ray.getDirection().normalized();
        Vector3f L_i = dirToLight.normalized();
        Vector3f R_i = 2 * (Vector3f::dot(N, L_i)) * N - L_i;
        R_i.normalized();
        shaded += specularColor * (pow(relu(Vector3f::dot(V, R_i)), shininess));
        bool isTexture;
        if(!isTexture) shaded += diffuseColor * relu(Vector3f::dot(L_i, N));
        else shaded += hit.getColor(isTexture) * relu(Vector3f::dot(L_i, N));
        shaded = lightColor * shaded;
        return shaded;
    }

    Vector3f sampleGlossy(const Vector3f &wi, const Vector3f &N) {
        float shininess = 50.0f; // Adjust shininess for glossy reflections
        Vector3f reflection = reflect(wi, N);
        Vector3f randomVector = randomUnitVector();
        Vector3f glossyDirection = (reflection + randomVector * (1.0f / shininess)).normalized();
        return glossyDirection;
    }

    Vector3f color; //颜色
    Vector3f diffuseColor;
    Vector3f specularColor; // 镜面反射系数

//参数部分参考已有代码：https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project
protected:
    float shininess; // 高光指数
    Vector3f emission;       // 发光系数
    float refr;              // 折射率
    Vector4f type;           // 种类
    // string type;
    Texture texture;   // 颜色纹理
    Texture bump;    //凹凸纹理

    Vector3f reflect(const Vector3f &I, const Vector3f &N) const {
        return I - 2 * Vector3f::dot(I, N) * N;
    }

    Vector3f randomUnitVector() const {
        float z = 2.0f * RND2 - 1.0f;
        float a = 2.0f * M_PI * RND2;
        float r = sqrt(1.0f - z * z);
        return Vector3f(r * cos(a), r * sin(a), z);
    }

    float relu(float x) { return std::max((float)0, x); }
};


#endif // MATERIAL_H
