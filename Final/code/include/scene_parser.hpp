#ifndef SCENE_PARSER_H
#define SCENE_PARSER_H

#include <cassert>
#include <vector>
#include <vecmath.h>
#include <sphere.hpp>

class Camera;
class Material;
class Object3D;
class Group;
class Sphere;
class Plane;
class Triangle;
class Transform;
class Mesh;
class Curve;
class RevSurface;

using std::vector;

#define MAX_PARSER_TOKEN_LENGTH 1024

class SceneParser {
public:

    SceneParser() = delete;
    SceneParser(const char *filename);

    ~SceneParser();

    Camera *getCamera() const {
        return camera;
    }

    Vector3f getBackgroundColor() const {
        return background_color;
    }

    int getNumMaterials() const {
        return num_materials;
    }

    Material *getMaterial(int i) const {
        assert(i >= 0 && i < num_materials);
        return materials[i];
    }

    int getNumLights() const {
        return num_lights;
    }

    Sphere *getLight(int i) const {
        assert(i >= 0 && i < num_lights);
        return lights[i];
    }

    Group *getGroup() const {
        return group;
    }

private:

    void parseFile();
    void parsePerspectiveCamera();
    void parseBackground();
    void parseLights();
    void parseMaterials();
    Material *parseMaterial();
    Object3D *parseObject(char token[MAX_PARSER_TOKEN_LENGTH]);
    Group *parseGroup();
    Sphere *parseSphere();
    Plane *parsePlane();
    Triangle *parseTriangle();
    Mesh *parseTriangleMesh();
    Transform *parseTransform();
    Curve *parseBezierCurve();
    Curve *parseBsplineCurve();
    RevSurface *parseRevSurface();

    int getToken(char token[MAX_PARSER_TOKEN_LENGTH]);

    Vector3f readVector3f();

    double readFloat();
    int readInt();

    FILE *file;
    Camera *camera;
    Vector3f background_color;
    int num_materials;
    Material **materials;
    Material *current_material;
    Sphere **lights;
    int num_lights;
    Group *group;
};

#endif // SCENE_PARSER_H
