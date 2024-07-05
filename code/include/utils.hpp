//参考smallpt与已有代码https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project
#pragma once
#include <cmath>
#include <cstdlib>
#include <random>

static std::mt19937 mersenneTwister;
static std::uniform_real_distribution<float> uniform;
#define RND1 (2.0 * uniform(mersenneTwister) - 1.0) //[-1, 1]均匀分布
#define RND2 (uniform(mersenneTwister)) //[0,1]均匀分布
#define PI (3.1415926536)
#define FLT_EPSILON	__FLT_EPSILON__

inline float ReLU(float x) { return x < 0 ? 0 : x > 1 ? 1 : x; }

inline float toFloat(float x) { return float(int(pow(ReLU(x), 1 / 2.2) * 255 + .5)) / 255.0; }

inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }

inline void normalizeTheta(float& theta) { //规范角度
        if (theta < 0.0) {
            theta += 2 * M_PI;
        } else if (theta >= 2 * M_PI) {
            theta = fmod(theta, 2 * M_PI);
        }
    }

inline void clampMu(float& mu) { //限制参数范围
    if (mu >= 1) {
        mu = 1.0 - FLT_EPSILON;
    } else if (mu <= 0) {
        mu = FLT_EPSILON;
    }
}