//参考已有代码https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project
#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP

#include <tuple>
#include "curve.hpp"
#include "object3d.hpp"
#include "triangle.hpp"
#include "constants.hpp"
#include "bound.hpp"
#include "utils.hpp"
#include <math.h>


class RevSurface : public Object3D {
    Curve* pCurve;
    Bound bbox;
    // Definition for drawable surface.
    typedef std::tuple<unsigned, unsigned, unsigned> Tup3u;
    std::vector<Triangle> triangles;

public:
    RevSurface(Curve* pCurve, Material* material)
        : pCurve(pCurve)
        , Object3D(material)
    {
        // Check flat.
        for (const auto& cp : pCurve->getControls()) {
            if (cp.z() != 0.0) {
                printf("Profile of revSurface must be flat on xy plane.\n");
                exit(0);
            }
        }
        Vector3f lo(-pCurve->radius, pCurve->y[0] - 3, -pCurve->radius);
        Vector3f hi(pCurve->radius, pCurve->y[1] + 3, pCurve->radius);
        bbox.set(lo, hi);
    }

    ~RevSurface() override { delete pCurve; }

    inline bool intersect(const Ray& r, Hit& h, float tmin) override {
        return newtonIntersect(r, h);
    }

private:
    bool newtonIntersect(const Ray& r, Hit& h) {
        float t, theta, mu;
        // 检测射线r是否与某个包围盒相交
        if (!bbox.intersect(r, h, t) || t > h.getT()) return false;
        getUV(r, t, theta, mu); // 计算射线与曲线的参数值theta和mu。
        Vector3f normal, point;
        // 利用牛顿迭代法求交点和法线
        if (!newton(r, t, theta, mu, normal, point)) return false;
        // 检查参数值的有效性
        if (!isnormal(mu) || !isnormal(theta) || !isnormal(t)) return false;
        if (t < 0 || mu < pCurve->range[0] || mu > pCurve->range[1] || t > h.getT()) return false;
        h.set(t, material, normal.normalized(), material->getColor(theta / (2 * M_PI), mu), point, material->getTexture().isTexture);
        return true;
    }

    bool newton(const Ray& r, float& t, float& theta, float& mu, Vector3f& normal, Vector3f& point) {
     // 牛顿迭代法 theta:位置参数 [0,2\pi] mu:位置比例 [0,1]
        // 求解交点和法线
        Vector3f dmu, dtheta;
        for (int i = 0; i < newton_depth; ++i) {
            normalizeTheta(theta);
            clampMu(mu);
            // 计算交点和切向量
            point = getPoint(theta, mu, dtheta, dmu);
            Vector3f f = r.getOrigin() + r.getDirection() * t - point;
            if (f.squaredLength() < NEWTON_EPS) return true;
            normal = Vector3f::cross(dmu, dtheta);
            // 牛顿迭代更新参数
            float D = Vector3f::dot(r.getDirection(), normal);
            t -= Vector3f::dot(dmu, Vector3f::cross(dtheta, f)) / D;
            mu -= Vector3f::dot(r.getDirection(), Vector3f::cross(dtheta, f)) / D;
            theta += Vector3f::dot(r.getDirection(), Vector3f::cross(dmu, f)) / D;
        }
        return false;
    }

    void getUV(const Ray& r, const float& t, float& theta, float& mu) {
        Vector3f pt(r.getOrigin() + r.getDirection() * t);
        theta = atan2(-pt.z(), pt.x()) + M_PI;
        mu = (pCurve->y[1] - pt.y()) / (pCurve->y[1] - pCurve->y[0]);
    }

    Vector3f getPoint(const float& theta, const float& mu, Vector3f& dtheta, Vector3f& dmu) {
        //根据给定的参数值theta和mu，计算参数曲面上对应点的位置，并计算该点处的切向量关于theta和mu的偏导数。
        Vector3f pt;
        Quat4f rot;
        rot.setAxisAngle(theta, Vector3f::UP);
        Matrix3f rotMat = Matrix3f::rotation(rot);
        CurvePoint cp;
        pCurve->evaluate(mu);
        for (int j = 0; j < pCurve->s.size(); ++j) {
            cp.V += pCurve->controls[pCurve->lsk + j] * pCurve->s[j];
            cp.T += pCurve->controls[pCurve->lsk + j] * pCurve->ds[j];
        }
        pt = rotMat * cp.V;
        dmu = rotMat * cp.T;
        dtheta = Vector3f(-cp.V.x() * sin(theta), 0, -cp.V.x() * cos(theta));
        return pt;
    }

};

#endif // REVSURFACE_HPP

