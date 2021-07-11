#include "ppm.hpp"

struct HPoint {
    Vector3f f, pos, normal, flux; 
    double r2; 
    unsigned int n;
    int x;
    int y;
    HPoint(const Vector3f &f_, const Vector3f &pos_, const Vector3f &normal_, int x_, int y_) : 
        f(f_), pos(pos_), normal(normal_), x(x_), y(y_) {}
};

vector<HPoint *> hitpoints;

unsigned int num_hash, num_photon;
double hash_s; 
vector<vector<HPoint *>> hash_grid; 
AABB hp_aabb;

inline unsigned int photonHash(const int ix, const int iy, const int iz) {
    return (unsigned int)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % num_hash;
}

void build_hash_grid(const int w, const int h) {
    hp_aabb.reset();
    for (auto hp : hitpoints) {
        hp_aabb.fit(hp->pos);
    }

    Vector3f ssize = hp_aabb.r() - hp_aabb.l();
    double irad = ((ssize.x() + ssize.y() + ssize.z()) / 3.0) / ((w + h) / 2.0) * 2.0;

    hp_aabb.reset();
    int vphoton = 0; 
    for (auto hp : hitpoints) {
        hp->r2 = irad * irad; 
        hp->n = 0;
        hp->flux = Vector3f();
        vphoton++; 
        hp_aabb.fit(hp->pos - irad); 
        hp_aabb.fit(hp->pos + irad);
    }

    hash_s = 1.0 / (irad * 2.0); 
    num_hash = vphoton; 

    hash_grid.resize(num_hash);
    for (unsigned int i = 0; i < num_hash; i++) hash_grid[i].clear();
    for (auto hp : hitpoints) { 
        Vector3f BMin = ((hp->pos - irad) - hp_aabb.l()) * hash_s;
        Vector3f BMax = ((hp->pos + irad) - hp_aabb.l()) * hash_s;
        for (int iz = abs(int(BMin.z())); iz <= abs(int(BMax.z())); iz++) {
            for (int iy = abs(int(BMin.y())); iy <= abs(int(BMax.y())); iy++) {
                for (int ix = abs(int(BMin.x())); ix <= abs(int(BMax.x())); ix++) {
                    int hv = photonHash(ix, iy, iz); 
                    hash_grid[hv].push_back(hp);
                }
            }
        }
    }
}
Ray genRay(Vector3f* flux, Sphere *light) {
    *flux = light->material()->ambientColor() * 
            light->material()->diffuseColor(Hit()) * 
            (kPi * 4.0 * light->radius() * light->radius()); // colored light
    double p = 2.0 * kPi * erand48(),t = 2.0 * acos(sqrt(1.0 - erand48()));
    double st = sin(t);
    Vector3f d(cos(p) * st,sin(p) * st, cos(t));
    return Ray(light->center() + (light->radius() + 2 * kEps) * d, d);
}

void ppmTraceRay(const Ray &r, int depth, Group *group, bool m, const Vector3f &flux, const Vector3f &flux_decay,
                 int x, int y) {
    Hit hit;
    if (!group->intersect(r, hit, kTmin) || ++depth > kMaxDepth) return;
    double t = hit.t();
    Vector3f intersect_pnt = r.pointAtParameter(t);
    Vector3f normal = hit.normal(), n = hit.into() ? normal : -normal;
    Material *material = hit.material();
    Vector3f f = material->diffuseColor(hit);
    double p = max(f);
    if (material->reflectType() == Material::kDiffuse) {
        double phi = 2 * kPi * erand48();
        double sq_sin_theta = erand48(), sin_theta = sqrt(sq_sin_theta);
        Vector3f w = normal, u, v;
        getXYbyZ(u, v, w);
        Vector3f ray_dir = (u * cos(phi) * sin_theta + v * sin(phi) * sin_theta + w * sqrt(1 - sq_sin_theta)).norm();
        if (m) {
            HPoint *hp = new HPoint(f * flux_decay, intersect_pnt, n, x, y);
            #pragma omp critical
            hitpoints.push_back(hp);
        } else {
            Vector3f hh = (intersect_pnt - hp_aabb.l()) * hash_s;
            int ix = abs(int(hh.x())), iy = abs(int(hh.y())), iz = abs(int(hh.z()));
            for (auto hitpoint : hash_grid[photonHash(ix, iy, iz)]) {
                Vector3f v = hitpoint->pos - intersect_pnt;
                if ((hitpoint->normal.dot(n) > kEps) && (v.dot(v) <= hitpoint->r2)) {
                    double g = (hitpoint->n * kAlpha + kAlpha) / (hitpoint->n * kAlpha + 1.0);
                    hitpoint->r2 = hitpoint->r2 * g; 
                    hitpoint->n++;
                    hitpoint->flux = (hitpoint->flux + hitpoint->f * flux * (1. / kPi)) * g;
                }
            }
            if (erand48() < p) ppmTraceRay(Ray(intersect_pnt, ray_dir), depth, group, m, f * flux * (1.0 / p), flux_decay, x, y);
        }
    } else if (material->reflectType() == Material::kSpecular) {
        Vector3f reflect_dir = r.direction() - n * 2 * n.dot(r.direction());
        ppmTraceRay(Ray(intersect_pnt, reflect_dir), depth, group, m, f * flux, f * flux_decay, x, y);
    } else /*if (material->reflectType() == Material::kRefractive)*/ {
        Vector3f reflect_dir = r.direction() - n * 2 * n.dot(r.direction());
        Ray reflect_ray(intersect_pnt, reflect_dir);
        bool into = hit.into();
        double n_src = kNAir, n_dest = material->refractiveIndex();
        double n_ratio = into ? n_src / n_dest : n_dest / n_src;
        double sin_src = r.direction().dot(normal);
        double sq_cos_dest = 1 - n_ratio * n_ratio * (1 - sin_src * sin_src);
        if (sq_cos_dest < 0) {
            return ppmTraceRay(reflect_ray, depth, group, m, flux, flux_decay, x, y);
        }
        Vector3f dest_ray_dir = (r.direction() * n_ratio - n * ((into ? 1 : -1) * (sin_src * n_ratio + sqrt(sq_cos_dest)))).norm();
        double a = n_dest - n_src, b = n_dest + n_src;
        double r_normal = a * a / (b * b);
        double c = 1 - (into ? -sin_src : dest_ray_dir.dot(n));
        double r_dir = r_normal + (1 - r_normal) * c * c * c * c * c;
        double r_p = r_dir;
        Ray refract_ray(intersect_pnt, dest_ray_dir);
        Vector3f color_flux_decay = f * flux_decay;
        if (m) {
            ppmTraceRay(reflect_ray, depth, group, m, flux, color_flux_decay * r_dir, x, y);
            ppmTraceRay(refract_ray, depth, group, m, flux, color_flux_decay * (1.0 - r_dir), x, y);
        } else {
            erand48() < r_p ? ppmTraceRay(reflect_ray, depth, group, m, flux, color_flux_decay, x, y) :
                              ppmTraceRay(refract_ray, depth, group, m, flux, color_flux_decay, x, y);
        }
    }
}

void init() {
    hp_aabb.reset();
    for (int i = 0; i < hitpoints.size(); ++i) {
        delete hitpoints[i];
    }

    hitpoints.clear();
    hash_grid.clear();
}

void singlePhotonMap(SceneParser &parser, int samples, Vector3f **pix, double start) {
    init();
    Camera *cam = parser.getCamera();
    int w = cam->getWidth(), h = cam->getHeight();
    Group *group = parser.getGroup();
    #pragma omp parallel for schedule(dynamic, 1)
    for (int y = 0; y < h; y++) {
        fprintf(stderr, "\rHitPointPass %5.2f%%, Time %.2fs", 100. * y / (h - 1), omp_get_wtime() - start);
        for (int x = 0; x < w; x++) {
            double r1 = 2 * erand48(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
            double r2 = 2 * erand48(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
            double cur_x = x + dx, cur_y = y + dy;
            ppmTraceRay(parser.getCamera()->generateRay(Vector2f(cur_x, cur_y)), 0, parser.getGroup(), true, Vector3f(), Vector3f(1, 1, 1), x, y);
        }
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Hitpoint Size %ld\n", hitpoints.size());
    build_hash_grid(w, h);
    num_photon = samples; 
    cerr << "Num photon: " << num_photon << endl;
    Vector3f vw = Vector3f(1, 1, 1);
    int num_lights = parser.getNumLights();
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < num_photon; i++) {
        double p = 100. * (i + 1) / num_photon;
        fprintf(stderr, "\rPhotonPass %5.2f%%, Time %.2fs", p, omp_get_wtime() - start); 
        Vector3f f;
        for(int j = 0; j < 1000; j++) {
            Sphere *light = parser.getLight(j % num_lights);
            Ray r = genRay(&f, light);
            ppmTraceRay(r, 0, group, false, f, vw, 0, 0);
        }
    }
    fprintf(stderr, "\n");
    for (auto hp : hitpoints) {
        pix[hp->x][hp->y] += hp->flux * (1.0 / (kPi * hp->r2 * num_photon * 1000.0));
    }
}

Image photonMap(SceneParser &parser, int samples) {
    Camera *cam = parser.getCamera();
    int w = cam->getWidth(), h = cam->getHeight();
    Image img(w, h);
    Vector3f **pix;
    pix = new Vector3f *[w];
    for (int i = 0; i < w; ++i) {
        pix[i] = new Vector3f[h];
    }

    auto start = omp_get_wtime();
    for (int i = 0; i < kPPMTurn; ++i) {
        fprintf(stderr, "PPM turn %d/%d...\n", i + 1, kPPMTurn);
        singlePhotonMap(parser, samples, pix, start);
    }

    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < h; ++j) {
            img.SetPixel(i, j, pix[i][j] / kPPMTurn);
        }
    }
    return img;
}