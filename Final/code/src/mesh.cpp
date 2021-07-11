#define TINYOBJLOADER_IMPLEMENTATION

#include "mesh.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <sstream>

bool Mesh::intersect(const Ray &r, Hit &h, double tmin) {
    return tree_->intersect(r, h, tmin);
}

Mesh::Mesh(const char *filename, Material *material) : Object3D(material), tree_(nullptr) {
    tinyobj::ObjReader reader;
    bool ret = reader.ParseFromFile(filename);
    if (!ret) {
        cerr << "Parse obj failed " << filename << endl;
        exit(1);
    }
    auto shapes = reader.GetShapes();
    auto attrib = reader.GetAttrib();
    for (size_t s = 0; s < shapes.size(); ++s) {
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); ++f) {
            size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);
            assert(fv == 3);
            Vector3f vertices[3], normal = Vector3f::ZERO;
            Vector2f texcoord[3];
            for (size_t v = 0; v < 3; ++v) {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
                tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
                tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];
                if (idx.normal_index >= 0) {
                    tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
                    tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
                    tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
                    normal = Vector3f(nx, ny, nz);
                }
                if (idx.texcoord_index >= 0) {
                    tinyobj::real_t tx = attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
                    tinyobj::real_t ty = attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
                    texcoord[v] = Vector2f(tx, ty);
                }
                vertices[v] = Vector3f(vx, vy, vz);
            }
            if (normal != Vector3f::ZERO) {
                triangle_vec_.push_back(new Triangle(vertices[0], vertices[1], vertices[2], normal, material_,
                                                     texcoord[0], texcoord[1], texcoord[2]));
            } else {
                triangle_vec_.push_back(new Triangle(vertices[0], vertices[1], vertices[2], material_,
                                                     texcoord[0], texcoord[1], texcoord[2]));
            }
            
            index_offset += fv;

            // per-face material
            shapes[s].mesh.material_ids[f];
        }
    }
    buildTree();
}

void Mesh::buildTree() {
    tree_ = new KDTree(triangle_vec_);
}
