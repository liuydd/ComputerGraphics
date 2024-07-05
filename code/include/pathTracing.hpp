//参考smallpt和已有代码https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project
#include<iostream>
#include "camera.hpp"
#include "constants.hpp"
#include "group.hpp"
#include "hit.hpp"
#include "image.hpp"
#include "light.hpp"
#include "ray.hpp"
#include "scene_parser.hpp"
#include "utils.hpp"
#include "omp.h"
using namespace std;

inline float Trowbridge_Reitz_GGX_D(Vector3f normal, Vector3f halfVector,float a)
{
    float a2=a*a;
    float NdotH =std::max(Vector3f::dot(normal,halfVector),0.0f);
    float NdotH2 =NdotH*NdotH;
    float nom=a2;
    float denom=(NdotH2*(a2-1.0f)+1.0f);
    denom=M_PI*denom*denom;
    return nom/std::max(denom,0.00001f);
}
inline float Schick_GGX(float NdotV,float k)
{
    float nom=NdotV;
    float denom=NdotV*(1.0f-k)+k;
    return nom/std::max(denom,0.00001f);
}
inline float Schick_GGXSmith_G(Vector3f N,Vector3f V,Vector3f L,float k)
{
    k = std::pow(k+1.0f,2.0f)/8.0f;
    float NdotV=std::max(Vector3f::dot(N,V),0.0f);
    float NdotL=std::max(Vector3f::dot(N,L),0.0f);
    float ggx1=Schick_GGX(NdotV,k);
    float ggx2=Schick_GGX(NdotL,k);
    return ggx1*ggx2;
}
    

float Schick_Fresnel_F(float cosTheta,float F0)
{
    return F0+(1.0-F0)*std::pow(1.0-cosTheta,5.0f);
}

class PathTracer{
public:
    const SceneParser& scene;
    int samps; //采样数
    const char* fout; //输出
    PathTracer(const SceneParser& scene, int samps, const char* fout);
    void render();
    Vector3f radiance(Ray ray, Group* group, int depth);
};

PathTracer::PathTracer(const SceneParser& scene, int samps, const char* fout): scene(scene), samps(samps), fout(fout){}

Vector3f PathTracer::radiance(Ray ray, Group* group, int depth){
    Vector3f color = Vector3f::ZERO;
    Vector3f cf = Vector3f::ONE;
    while(true){ //采用迭代的方式
        if(cf.max() < 1e-3) return color;
        // Russian Roulette终止策略
        float p = max(cf.x(), max(cf.y(), cf.z()));
        if (depth > 5) {
            if (RND2 > p) return color;
            cf = cf / p;
        }

        Hit hit;
        if(!group->intersect(ray, hit, 1e-5)){ //如果不相交 //这里tmin不能为0，否则基本是全黑，因为会导致光线在表面上反复hit！
            return color;
        }
        else{
            // 获取相交点的材质信息
            Material* material = hit.getMaterial();
            bool is_texture;
            Vector3f hit_color = hit.getColor(is_texture);
            if (!is_texture)
                hit_color = material->color;
            // 获取材质的自发光颜色、法线和下一条光线的起点
            color += cf * material->getEmissionColor();
            cf = cf * hit_color;
            // cout<<"cf"<<cf.x()<<" "<<cf.y()<<" "<<cf.z()<<"refcolr"<<hit_color.x()<<" "<<hit_color.y()<<" "<<hit_color.z()<<"depth "<<depth<<endl;
            Vector3f hit_emission = material->getEmissionColor();
            Vector3f N = hit.getNormal().normalized();
            Vector3f next_origin = ray.getOrigin() + hit.getT() * ray.getDirection();
            Vector3f ray_direction = ray.getDirection().normalized();

            // 根据随机数生成的类型决策来确定下一条光线的方向和类型
            float random_type = RND2;
            float b = Vector3f::dot(ray_direction, N);
            Vector3f next_direction;

            // NEE（Next Event Estimation）
            // 直接采样光源位置并估计其对当前交点的贡献
            // for (int i = 0; i < scene.getNumLights(); ++i) {
            //     Light* light = scene.getLight(i);
            //     Vector3f light_dir, light_col;
            //     light->getIllumination(next_origin, light_dir, light_col);
            //     float distance_squared = light_dir.squaredLength();
            //     light_dir.normalized();
            //     Ray shadow_ray(next_origin, light_dir);
            //     Hit shadow_hit;
            //     if (!group->intersect(shadow_ray, shadow_hit, 1e-5)) {
            //         color += cf * light_col * max(0.0f, Vector3f::dot(N, light_dir)) / distance_squared;
            //     }
            // }

            //根据类型计算下一条光线的方向
            if(random_type < material->getType().x()){ //漫反射
                Vector3f z_ = Vector3f::cross(ray_direction, N);
                Vector3f x_ = Vector3f::cross(z_, N);
                z_.normalize();
                x_.normalize();
                if (b < 0) next_direction = RND1 * z_ + RND1 * x_ + RND2 * N;
                else next_direction = RND1 * z_ + RND1 * x_ - RND2 * N;
                next_direction.normalize();
                ray = Ray(next_origin, next_direction);
            }
            else if(random_type < material->getType().x() + material->getType().y()){ //镜面反射
                next_direction = ray_direction - N * (2 * b); //光线根据法向做轴变换
                next_direction.normalized();
                ray = Ray(next_origin, next_direction);
            }
            else if(random_type < material->getType().x() + material->getType().y() + material->getType().z()){ //折射
                float n = material->getRefr(); // 折射率
                float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n)); // 被反射光线比例
                if (b > 0) N.negate();
                else n = 1.0 / n;
                
                float cos1 = -Vector3f::dot(N, ray_direction); // 入射角theta余弦值
                float cos2 = 1.0 - n * n * (1.0 - cos1 * cos1); // 出射角theta2 cos(theta2)的平方
                Vector3f reflect = (ray_direction + N * (cos1 * 2));
                if (cos2 < 0) { // 全反射
                    ray = Ray(next_origin, reflect);
                } else {
                    // 用菲涅尔项计算光线在界面上发生反射和折射的概率。
                    float Rprob = R0 + (1.0 - R0) * pow(1.0 - cos1, 5.0);
                    Vector3f refrac = (ray_direction * n) + (N * (n * cos1 - sqrt(cos2)));
                    refrac.normalized();
                    if (cos2 > 0 && RND2 > Rprob) ray = Ray(next_origin, refrac);
                    else ray = Ray(next_origin, reflect);
                }
            }
            else{ //glossy
                // glossy部分与王伊荷和陈畅讨论过解决方法，以及参考相关公式
                Vector3f wo = ray.getDirection();
                double r1 = 2*M_PI*RND2, r2 = RND2, r2s = sqrt(r2);
                Vector3f nl = Vector3f::dot(N, ray.getDirection()) < 0? N: N*(-1);
                Vector3f w = nl;
                Vector3f mm = Vector3f::cross(Vector3f(0, 1, 0), w);
                Vector3f mn = Vector3f::cross(Vector3f(1, 0, 0), w);
                Vector3f u = (fabs(w.x())>0.1? mm : mn).normalized();
                Vector3f v = Vector3f::cross(w, u);
                Vector3f next_direction = (u * cos(r1)*r2s + v*sin(r1)*r2s+w*sqrt(1-r2)).normalized();
                ray = Ray(next_origin, next_direction);
                float alpha = Vector3f::dot(N, wo);
                float roughness=0.02f;
                Vector3f V = next_direction;
                Vector3f L = -wo;
                Vector3f H =(V+L).normalized();
                float D =Trowbridge_Reitz_GGX_D(N,H,roughness);
                float G =Schick_GGXSmith_G(N,V,L,roughness);
                float F = Schick_Fresnel_F(alpha,0.50f);
                float diffsue=Vector3f::dot(N,V)+0.5f;
                float divide=1.0f/(4*std::max(Vector3f::dot(N,L),0.0001f)*std::max(Vector3f::dot(N,V),0.0001f));
                float Specular = D*F*G*divide;
                hit.getMaterial()->diffuseColor = (0.9, 0.9, 0.9);
                Vector3f new_color = diffsue*(hit.getMaterial()->diffuseColor) / M_PI + (hit.getMaterial()->specularColor)*Specular;
                color += new_color;
            }
            depth ++;
        }
    }
    return color;
}

void PathTracer::render(){
    Camera* camera = scene.getCamera();
    Image outImg(camera->getWidth(), camera->getHeight());
    Vector3f background = scene.getBackgroundColor();
    Group* group = scene.getGroup();
    #pragma omp parallel for schedule(dynamic, 1)  // OpenMP
    time_t start = time(NULL);
    for (int y = 0; y < scene.getCamera()->getHeight(); ++y) {
        float elapsed = (time(NULL) - start), progress = (1. + y) / camera->getHeight();
        fprintf(stderr, "\rRendering (%d spp) %5.2f%% Time: %.2f/%.2f sec",
                    samps, progress * 100., elapsed, elapsed / progress);
        for (int x = 0; x < scene.getCamera()->getWidth(); ++x) {
            Vector3f color(0, 0, 0);
            for (int s = 0; s < samps; ++s) {
                // double u = (x + drand48()) / scene.getCamera()->getWidth();
                // double v = (y + drand48()) / scene.getCamera()->getHeight();
                Ray r = camera->generateRay(Vector2f(x + RND1 , y + RND1 )); //使用画布网格随机扰动的操作，实现抗锯齿效果
                color += radiance(r, group, 0);
            }
            color = color / float(samps);
            // cout<<color[0]<<" "<<color[1]<<" "<<color[2]<<endl;
            outImg.SetPixel(x, y, Vector3f(toFloat(color.x()), toFloat(color.y()), toFloat(color.z())));
        }
    }
    outImg.SaveBMP(fout);
}

