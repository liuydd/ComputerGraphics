//参考已有代码https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project
#ifndef TEXTURE_H
#define TEXTURE_H
#include <string>

#include <vecmath.h>
using std::string;
struct Texture {  // 纹理
    unsigned char *picture;
    int width, height, channel;
    bool isTexture;
    Texture(const char *textureFile);
    Vector3f getColor(float u, float v) const;
    Vector3f getColor(int u, int v) const;
    float getDisturb(float u, float v, Vector2f &grad) const;
    inline int getIndex(float u, float v) const;
    inline float getGray(int index) const;
};

#endif  // !TEXTURE_H