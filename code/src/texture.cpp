//参考已有代码https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project
#define STB_IMAGE_IMPLEMENTATION
#include <vecmath.h>
#include <iostream>
#include "stb_image.h"
#include <cstring>
#include "texture.h"
using namespace std;

float compare_zero(float u){
    if(u < 0) return 0;
    else return u;
}

int Texture::getIndex(float u, float v) const {
    int x = int(u * width + width) % width;
    int y = int(v * height + height) % height;
    // 将像素坐标限制在图像边界内
    int index = (y * width + x) * channel;
    return index;
}

float Texture::getGray(int index) const { return (picture[index] / 255. - 0.5) * 2; }

Texture::Texture(const char *textureFile) {
    if (strlen(textureFile) > 0) {
        picture = stbi_load(textureFile, &width, &height, &channel, 0);
        isTexture = true;
    } else {
        isTexture = false;
    }
}

Vector3f Texture::getColor(int u, int v) const {
    u = compare_zero(u);
    // cout<<"u "<<u<<endl;
    // cout<<"width "<<width<<endl;
    u = u > width - 1 ? width - 1 : u;
    // cout<<"u "<<u<<endl;
    v = compare_zero(v);
    // cout<<"v "<<v<<endl;
    // cout<<"height "<<height<<endl;
    v = v > height - 1 ? height - 1 : v;
    // cout<<"v "<<v<<endl;
    int index = (v * width + u) * channel;
    return Vector3f(picture[index + 0], picture[index + 1], picture[index + 2]) / 255.0;
}

Vector3f Texture::getColor(float u, float v) const {
     // 将小数部分限制在 [0, 1) 范围内
    u -= int(u);
    v -= int(v);
    u = u < 0 ? 1 + u : u;
    v = v < 0 ? 1 + v : v;
    // 使用双线性插值计算纹理颜色
    u = u * width;
    v = height * (1 - v);
    int iu = (int)u, iv = (int)v;
    Vector3f ret = Vector3f::ZERO;
    float alpha = u - iu;
    float beta = v - iv;
    // 使用双线性插值计算纹理颜色
    ret += (1 - alpha) * (1 - beta) * getColor(iu, iv) + alpha * (1 - beta) * getColor(iu + 1, iv) + (1 - alpha) * beta * getColor(iu, iv + 1) + alpha * beta * getColor(iu + 1, iv + 1);
    return ret;
}

float Texture::getDisturb(float u, float v, Vector2f &grad) const {
    if(!isTexture) return 0.0;
    float disturb = getGray(getIndex(u, v));
    float du = 1.0 / width, dv = 1.0 / height;
    //计算梯度
    grad[0] = width * (getGray(getIndex(u + du, v)) - getGray(getIndex(u - du, v))) * 0.5;
    grad[1] = height * (getGray(getIndex(u, v + dv)) - getGray(getIndex(u, v - dv))) * 0.5;
    return disturb;
}
