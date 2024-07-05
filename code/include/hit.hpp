//独立实现
#ifndef HIT_H
#define HIT_H

#include <vecmath.h>
#include "ray.hpp"

class Material;

class Hit {
public:
    Vector3f color;
    Vector3f p;
    bool isTexture;

    // constructors
    Hit() {
        material = nullptr;
        t = 1e38;
    }

    Hit(float _t, Material *m, const Vector3f &n) {
        t = _t;
        material = m;
        normal = n;
    }

    Hit(const Hit &h) {
        t = h.t;
        material = h.material;
        normal = h.normal;
    }

    // destructor
    ~Hit() = default;

    float getT() const {
        return t;
    }

    Material *getMaterial() const {
        return material;
    }

    const Vector3f &getNormal() const {
        return normal;
    }

    Vector3f getColor(bool& _isTexture) const{
        _isTexture = isTexture;
        return color;
    }

    void set(float _t, Material *m, const Vector3f &n, const Vector3f &c = Vector3f::ZERO, const Vector3f &_p = Vector3f::ZERO, bool _isTexture = false) {
        t = _t;
        material = m;
        normal = n;
        color = c;
        p = _p;
        isTexture = _isTexture;
    }

private:
    float t;
    Material *material;
    Vector3f normal;

};

inline std::ostream &operator<<(std::ostream &os, const Hit &h) {
    os << "Hit <" << h.getT() << ", " << h.getNormal() << ">";
    return os;
}

#endif // HIT_H
