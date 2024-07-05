//独立实现
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>

#include "scene_parser.hpp"
#include "image.hpp"
#include "camera.hpp"
#include "group.hpp"
#include "light.hpp"
#include "pathTracing.hpp"

#include <string>
#include <chrono>

using namespace std;

int main(int argc, char *argv[]) {
    for (int argNum = 1; argNum < argc; ++argNum) {
        std::cout << "Argument " << argNum << " is: " << argv[argNum] << std::endl;
    }

    if (argc != 5) {
        cout << "Usage: ./bin/PA1 <input scene file> <output bmp file>" << endl;
        return 1;
    }
    string inputFile = argv[1];
    string outputFile = argv[2];  // only bmp is allowed.
    int samples = atoi(argv[3]); //采样数
    string method = argv[4]; //方式

    auto start = std::chrono::high_resolution_clock::now();

    // TODO: Main RayCasting Logic
    // First, parse the scene using SceneParser.
    SceneParser sceneParser(inputFile.c_str());
    if(method == "pt"){
        printf("Method: path tracing\n");
        PathTracer pt(sceneParser, samples, outputFile.c_str());
        pt.render();
    }
    else{
        Camera *camera = sceneParser.getCamera();
        Image img(camera->getWidth(), camera->getHeight());
        // Then loop over each pixel in the image, shooting a ray
        // through that pixel and finding its intersection with
        // the scene.  Write the color at the intersection to that
        // pixel in your output image.
        for(int x = 0; x < camera->getWidth(); ++x){
            for(int y = 0; y< camera->getHeight(); ++y){
                Ray camRay = camera->generateRay(Vector2f(x, y));
                Group * baseGroup = sceneParser.getGroup();
                Hit hit;
                bool isIntersect = baseGroup->intersect(camRay, hit, 0);
                if(isIntersect){
                    Vector3f finalColor = Vector3f::ZERO;
                    for(int li = 0; li < sceneParser.getNumLights(); ++li){
                        Light *light = sceneParser.getLight(li);
                        Vector3f L, lightColor;
                        light->getIllumination(camRay.pointAtParameter(hit.getT()), L, lightColor);
                        finalColor += hit.getMaterial()->Shade(camRay, hit, L, lightColor);
                    }
                    img.SetPixel(x, y, finalColor);
                }
                else{
                    img.SetPixel(x, y, sceneParser.getBackgroundColor());
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    cout << "运行时间: " << duration << " 秒" << endl;
    // cout << "Hello! Computer Graphics!" << endl;
    // img.SaveBMP(outputFile.c_str());
    return 0;
}

