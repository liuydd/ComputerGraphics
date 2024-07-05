#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include <vecmath.h>
#include <float.h>
#include <cmath>
#include "utils.hpp"

//PA1框架
class Camera {
public:
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH) {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up).normalized();
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

// TODO: Implement Perspective camera
// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera {

public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
            const Vector3f &up, int imgW, int imgH, float angle) : Camera(center, direction, up, imgW, imgH) {
        // angle is in radian.(弧度制)
        cx = imgW / 2;
        cy = imgH / 2;
        fx = cx / tan(angle / 2);
        fy = cy / tan(angle / 2);
    }

    Ray generateRay(const Vector2f &point) override {
        // 
        Vector3f d_rc = Vector3f((point[0] - cx) / fx, (cy - point[1]) / fy, 1).normalized();
        Matrix3f R(horizontal, -up, direction);
        return Ray(center, R * d_rc);
    }

protected:
    float fx, fy, cx, cy;
};

//参考已有代码：https://github.com/zrporz/THU-CG-2023
class DofCamera : public Camera {
    float fxy;
    float distance; // 焦距
    float radius; // 光圈半径
public:
    DofCamera(const Vector3f& center, const Vector3f& direction,
        const Vector3f& up, int imgW, int imgH, float angle, float distance, float radius)
        : Camera(center, direction, up, imgW, imgH), distance(distance), radius(radius)
    {
        // angle is in radian.(弧度制)
        cx = imgW / 2;
        cy = imgH / 2;
        fx = cx / tan(angle / 2);
        fy = cy / tan(angle / 2);
        fxy = imgH / (2 * tan(angle / 2) * (distance + 1));
    }

    Ray generateRay(const Vector2f& point) override
    {
        float csx = distance * (point.x() - cx) / fx;
        float csy = distance * (point.y() - cy) / fy;
        //在光圈上随机采样一个点
        float theta = 2 * M_PI * RND2;
        float r = RND2 * radius;
        float lensX = r * cos(theta);
        float lensY = r * sin(theta);
        Vector3f point2 = lensY * up + lensX * horizontal;
        Vector3f dir(csx - point2.x(), -csy - point2.y(), distance);
        Vector3f lensPoint = center + horizontal * point2.x() - up * point2.y();
        Vector3f finalDir = Matrix3f(horizontal, -up, direction) * dir.normalized();
        return Ray(lensPoint, finalDir);
    }

protected:
    float fx;
    float fy;
    float cx;
    float cy;
};

#endif //CAMERA_H
