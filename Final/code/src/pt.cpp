#include "pt.hpp"

Vector3f traceRay(const Ray &r, int depth, Group *group) {
    Hit hit;
    if (!group->intersect(r, hit, kTmin)) return Vector3f();
    double t = hit.t();
    Vector3f intersect_pnt = r.pointAtParameter(t);
    Vector3f normal = hit.normal(), n = hit.into() ? normal : -normal;
    Material *material = hit.material();
    Vector3f f = material->diffuseColor(hit);
    double p = max(f);
    if (++depth > 5) {
        if (erand48() < p) f *= (1 / p);
        else return material->ambientColor(); //R.R.
    }
    if (depth > kMaxDepth) {
        return material->ambientColor();
    }
    if (material->reflectType() == Material::kDiffuse) {
        double phi = 2 * kPi * erand48();
        double sq_sin_theta = erand48(), sin_theta = sqrt(sq_sin_theta);
        Vector3f w = normal, u, v;
        getXYbyZ(u, v, w);
        Vector3f ray_dir = (u * cos(phi) * sin_theta + v * sin(phi) * sin_theta + w * sqrt(1 - sq_sin_theta)).norm();
        return material->ambientColor() + f * traceRay(Ray(intersect_pnt, ray_dir), depth, group);
    } else if (material->reflectType() == Material::kSpecular) {
        Vector3f reflect_dir = r.direction() - n * 2 * n.dot(r.direction());
        return material->ambientColor() + f * traceRay(Ray(intersect_pnt, reflect_dir), depth, group);
    } else /*if (material->reflectType() == Material::kRefractive)*/ {
        Vector3f reflect_dir = r.direction() - n * 2 * n.dot(r.direction());
        Ray reflect_ray(intersect_pnt, reflect_dir);
        bool into = hit.into();
        double n_src = kNAir, n_dest = material->refractiveIndex();
        double n_ratio = into ? n_src / n_dest : n_dest / n_src;
        double sin_src = r.direction().dot(normal);
        double sq_cos_dest = 1 - n_ratio * n_ratio * (1 - sin_src * sin_src);
        if (sq_cos_dest< 0) {
            return material->ambientColor() + f * traceRay(reflect_ray, depth, group);
        }
        Vector3f dest_ray_dir = (r.direction() * n_ratio - n * ((into ? 1 : -1) * (sin_src * n_ratio + sqrt(sq_cos_dest)))).norm();
        double a = n_dest - n_src, b = n_dest + n_src;
        double r_normal = a * a / (b * b);
        double c = 1 - (into ? -sin_src : dest_ray_dir.dot(n));
        double r_dir = r_normal + (1 - r_normal) * c * c * c * c * c;
        double t_r = 1 - r_dir, p_reflect = .25 + .5 * r_dir, r_p = r_dir / p_reflect, t_p = t_r / (1 - p_reflect);
        Ray refract_ray(intersect_pnt, dest_ray_dir);
        return material->ambientColor() + f * (depth > 2 ? (erand48() < p_reflect ? 
            traceRay(reflect_ray, depth, group) * r_p : 
            traceRay(refract_ray, depth, group) * t_p) : 
            traceRay(reflect_ray, depth, group) * r_dir + traceRay(refract_ray, depth, group) * t_r);
    }
}

Vector3f tracePixel(SceneParser *parser, int x, int y, int samples) {
    Vector3f pixel;
    for (int sy = 0; sy < 2; sy++) {
        for (int sx = 0; sx < 2; sx++) {
            Vector3f r;
            for (int s = 0; s < samples; s++) {
                double r1 = 2 * erand48(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                double r2 = 2 * erand48(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                double cur_x = x + sx - 0.5 + dx, cur_y = y + sy - 0.5 + dy;
                r = r + traceRay(parser->getCamera()->generateRay(Vector2f(cur_x, cur_y)), 0, parser->getGroup()) * (1. / samples);
            }
            pixel = pixel + clamp(r) * .25;
        }
    }
    return pixel;
}

Image pathTrace(SceneParser &parser, int samples) {
    Camera *cam = parser.getCamera();
    int w = cam->getWidth(), h = cam->getHeight();
    Image img(w, h);
    Group *group = parser.getGroup();
    auto start = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 1)
    for (int y = 0; y < h; y++) {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%, Time %.2fs", samples * 4, 100. * y / (h - 1), omp_get_wtime() - start);
        for (int x = 0; x < w; x++) {
            img.SetPixel(x, y, tracePixel(&parser, x, y, samples));
        }
    }
    return img;
}