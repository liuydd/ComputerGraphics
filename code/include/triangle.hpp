//参考已有代码https://github.com/zrporz/THU-CG-2023
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>
using namespace std;

// TODO: implement this class and add more fields as necessary,
class Triangle: public Object3D {

public:
	Triangle() = delete;

    // a b c are three vertex positions of the triangle
	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c, Material* m) : Object3D(m), vertices({a, b, c}) {
		normal = Vector3f::cross(a-b, a-c).normalized();
		Vector3f edge1 = a - c;
        Vector3f edge2 = b - a;
        d = Vector3f::dot(normal, a);
		e1sq = Vector3f::dot(edge1, edge1);
		e1e2 = Vector3f::dot(edge1, edge2);
		e2sq = Vector3f::dot(edge2, edge2);
		inverDeno = e1sq * e2sq - e1e2 * e1e2;
		is_norm = false;
		is_texture = false;
	}


	bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
		Vector3f origin(ray.getOrigin());
        Vector3f dir(ray.getDirection());
		float dir_len = dir.length();
        float x = Vector3f::dot(normal, dir);
        if(fabs(x) < 1e-6) return false; //光线与平面平行
		float d = Vector3f::dot(vertices[0], normal);
		float t = (d - Vector3f::dot(normal, origin)) / Vector3f::dot(normal, dir);
        if (t < tmin || t > hit.getT()) return false;
		Vector3f next_origin(origin + t * dir);
		if(!isinTriangle(next_origin)) return false;
		t = t / dir_len;
		if(inverDeno == 0) return false;
		Vector3f edge1 = vertices[0] - vertices[2];
		Vector3f edge2 = vertices[1] - vertices[0];
		Vector3f ans = origin + t * dir - vertices[0];
		float e1ans = Vector3f::dot(edge1, ans);
		float e2ans = Vector3f::dot(edge2, ans);
		float u = (-e2sq * e1ans + e1e2 * e2ans) / inverDeno;
		if(u < 0 || u > 1) return false;
		float v = (e1sq * e2ans - e1e2 * e1ans) / inverDeno;
		if(v < 0 || v > 1) return false;
		if(u + v > 1) return false;
		float uu = 0, vv = 0;
		getUV(next_origin, uu, vv);
		Vector3f new_normal = getNorm(next_origin);
		hit.set(t, this->material, new_normal, material->getColor(uu, vv), next_origin, is_texture);
		// hit.set(t, material, normal);
        return true;
	}
	Vector3f normal;
	Vector3f vertices[3];
	// Vector2f bound[2];

protected:
	bool isinTriangle(const Vector3f& c){
		Vector3f edge0(vertices[0] - c);
        Vector3f edge1(vertices[1] - c);
        Vector3f edge2(vertices[2] - c);
        Vector3f c0(Vector3f::cross(edge0, edge1));
        Vector3f c1(Vector3f::cross(edge1, edge2));
		Vector3f c2(Vector3f::cross(edge2, edge0));
        return Vector3f::dot(c0, normal) >= -1e-6 && Vector3f::dot(c1, normal) >= -1e-6 && Vector3f::dot(c2, normal) >= -1e-6;
	}
	void getUV(const Vector3f& p, float& u, float& v) {	//计算p点的纹理坐标
        if (!is_texture)
			return;
        Vector3f va = (vertices[0] - p);
		Vector3f vb = (vertices[1] - p);
		Vector3f vc = (vertices[2] - p);
        float ra = Vector3f::cross(vb, vc).length();
    	float rb = Vector3f::cross(vc, va).length();
    	float rc = Vector3f::cross(va, vb).length();
        Vector2f uv = (ra * at + rb * bt + rc * ct) / (ra + rb + rc);
        u = uv.x();
        v = uv.y();
    }

	Vector3f getNorm(const Vector3f& p) {	//获取p点处的法线向
		if(!is_norm)
			return this->normal;
		Vector3f va = (vertices[0] - p);
		Vector3f vb = (vertices[1] - p);
		Vector3f vc = (vertices[2] - p);
		float ra = Vector3f::cross(vb, vc).length();
    	float rb = Vector3f::cross(vc, va).length();
    	float rc = Vector3f::cross(va, vb).length();
		return (ra * an + rb * bn + rc * cn).normalized();
	}

	float d;
	float e1sq;
	float e1e2;
	float e2sq;
	float inverDeno;	// 用于计算相交参数的倒数
	bool is_texture;
	bool is_norm;
	Vector2f at, bt, ct;	// 纹理坐标的权重
    Vector3f an, bn, cn;
};

#endif //TRIANGLE_H
