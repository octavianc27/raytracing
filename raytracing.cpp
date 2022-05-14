#include <iostream>
#include <cmath>
#include <vector>
#include "Vector3f.h"
#include <fstream>
#include <random>
#include <omp.h>
using namespace std;

// 1- ray tracing
// 2- Lambertian shading, point-like source
// 3- antialiasing
// 4- shadows
// 5- OOP for planes and triangles
#pragma region // Hittable
// All objects who interact with light
class Hittable
{
protected:
    Vector3f object_color;
    float diffusion_coeff;
    float specular_coeff;
    float phong_exp;

public:
    Vector3f get_color();
    float get_diffc();
    float get_specc();
    float get_phongexp();
    virtual float distance_from_p(Vector3f, Vector3f) = 0;
    virtual Vector3f get_normal(Vector3f) = 0;
};

// Get an object's color
Vector3f Hittable::get_color()
{
    return object_color;
}

// Get an object's diffuse coefficient
float Hittable::get_diffc()
{
    return diffusion_coeff;
}

// Get specular coefficient
float Hittable::get_specc()
{
    return specular_coeff;
}

// Get phong's exponent
float Hittable::get_phongexp()
{
    return phong_exp;
}
#pragma endregion
#pragma region // Sphere
class Sphere : public Hittable
{
    Vector3f center;
    float radius;

public:
    Sphere(Vector3f, float, Vector3f, float, float, float);
    void set_values(Vector3f, float, Vector3f, float, float, float);
    float distance_from_p(Vector3f, Vector3f);
    Vector3f get_normal(Vector3f);
};

// Initialise sphere
Sphere::Sphere(Vector3f c, float r, Vector3f col, float dif_c, float spec_c, float p)
{
    center = c;
    radius = r;
    object_color = col;
    diffusion_coeff = dif_c;
    specular_coeff = spec_c;
    phong_exp = p;
}

// Set sphere
void Sphere::set_values(Vector3f c, float r, Vector3f col, float dif_c, float spec_c, float p)
{
    center = c;
    radius = r;
    object_color = col;
    diffusion_coeff = dif_c;
    specular_coeff = spec_c;
    phong_exp = p;
}

// Get the normal direction to a point on the sphere
Vector3f Sphere::get_normal(Vector3f point)
{
    return (point - center).normalized();
}

// Return distance from a point to sphere on a given normalised direction, or -1 if no intersection
float Sphere::distance_from_p(Vector3f point, Vector3f dir)
{
    Vector3f q = point - center;
    float a = 1;
    float b = 2 * dir.innerProduct(q);
    float c = q.length_squared() - radius * radius;
    float d = b * b - 4 * a * c;
    // no intersection
    if (d < 0)
    {
        return -1;
    }
    // intersection
    else
    {
        float t1 = (-b + sqrt(d)) / (2 * a);
        float t2 = (-b - sqrt(d)) / (2 * a);
        if (t1 > t2)
            swap(t1, t2);
        if (t1 >= 0)
            return t1;
        return t2;
    }
}
#pragma endregion
#pragma region // Plane
class Plane : public Hittable
{
    Vector3f normal;
    Vector3f point_on_plane;

public:
    Plane(Vector3f, Vector3f, Vector3f, float, float, float);
    void set_values(Vector3f, Vector3f, Vector3f, float, float, float);
    float distance_from_p(Vector3f, Vector3f) override;
    Vector3f get_normal(Vector3f) override;
};
// Initialise plane
Plane::Plane(Vector3f n, Vector3f po, Vector3f col, float dif_c, float spec_c, float p)
{
    normal = n.normalized();
    point_on_plane = po;
    object_color = col;
    diffusion_coeff = dif_c;
    specular_coeff = spec_c;
    phong_exp = p;
}

// Set plane
void Plane::set_values(Vector3f n, Vector3f po, Vector3f col, float dif_c, float spec_c, float p)
{
    normal = n.normalized();
    point_on_plane = po;
    object_color = col;
    diffusion_coeff = dif_c;
    specular_coeff = spec_c;
    phong_exp = p;
}

// Get the normal direction to a plane
Vector3f Plane::get_normal(Vector3f point)
{
    return normal;
}

// Return distance from a point to the plane on a given direction, or -1 if no intersection
float Plane::distance_from_p(Vector3f point, Vector3f dir)
{
    // t=-(n.point)/(n.direction)
    float term1 = normal.innerProduct(point - point_on_plane);
    float term2 = normal.innerProduct(dir.normalized());
    return -term1 / term2;
}
#pragma endregion
#pragma region // Light source
class Source
{
    Vector3f location;
    float intensity;

public:
    Source();
    Source(Vector3f, float);
    Vector3f get_location();
    float get_intensity();
};
// Constructor
Source::Source(Vector3f loc, float intt)
{
    location = loc;
    intensity = intt;
}

Source::Source()
{
    location = Vector3f();
    intensity = 0;
}

Vector3f Source::get_location()
{
    return location;
}

float Source::get_intensity()
{
    return intensity;
}
#pragma endregion

#pragma region // Triangle
class Triangle : public Hittable
{
    Vector3f vertex_1, vertex_2, vertex_3;

public:
    Triangle(Vector3f, Vector3f, Vector3f, Vector3f, float, float, float);
    void set_values(Vector3f, Vector3f, Vector3f, Vector3f, float, float, float);
    float distance_from_p(Vector3f, Vector3f) override;
    Vector3f get_normal(Vector3f) override;
};
// constructor
Triangle::Triangle(Vector3f v1, Vector3f v2, Vector3f v3, Vector3f col, float dif_c, float spec_c, float p)
{
    vertex_1 = v1;
    vertex_2 = v2;
    vertex_3 = v3;
    object_color = col;
    diffusion_coeff = dif_c;
    specular_coeff = spec_c;
    phong_exp = p;
}
// setter
void Triangle::set_values(Vector3f v1, Vector3f v2, Vector3f v3, Vector3f col, float dif_c, float spec_c, float p)
{
    vertex_1 = v1;
    vertex_2 = v2;
    vertex_3 = v3;
    object_color = col;
    diffusion_coeff = dif_c;
    specular_coeff = spec_c;
    phong_exp = p;
}
// normal
Vector3f Triangle::get_normal(Vector3f point)
{
    Vector3f v1 = vertex_2 - vertex_1;
    Vector3f v2 = vertex_3 - vertex_1;
    Vector3f norm = v1.crossProduct(v2);
    return norm.normalized() * -1;
}
// distance from p - Moller Trumbore algorithm
float Triangle::distance_from_p(Vector3f point, Vector3f dir)
{
    const float epsilon = 0.00001;
    Vector3f edge1 = vertex_2 - vertex_1;
    Vector3f edge2 = vertex_3 - vertex_1;
    Vector3f h = dir.crossProduct(edge2);
    float a = edge1.innerProduct(h);
    // parralel ray
    // cout<<abs(a)<<endl;
    if (abs(a) < epsilon)
    {
        // cout << "caz1" << endl;
        return -1;
    }
    float f = 1.0 / a;
    Vector3f s = point - vertex_1;
    float u = f * s.innerProduct(h);
    if (u < 0.0 || u > 1.0)
    {
        // cout << "caz2" << endl;
        return -1;
    }
    Vector3f q = s.crossProduct(edge1);
    float v = f * dir.innerProduct(q);
    if (v < 0 || u + v > 1)
        return -1;
    float t = f * edge2.innerProduct(q);
    if (t > epsilon)
        return t;
    return -1;
}
#pragma endregion

int main()
{
    // antialiasing parameters
    int n_iter = 24;
    srand(time(0));

    // position of the origin
    Vector3f origin(0, 2, 0);

    // light sources
    float ambient_light = 0.1;
    int n_sources = 3;
    vector<Source> light_sources(n_sources);
    light_sources[0] = Source(Vector3f(3, 2, 0), 20);
    light_sources[1] = Source(Vector3f(-3, 6, 5), 8);
    light_sources[2] = Source(Vector3f(0, 3, -1), 10);

    // Some colors
    Vector3f black(0, 0, 0);
    Vector3f red(1, 0, 0);
    Vector3f green(0, 1, 0);
    Vector3f blue(0, 0, 1);
    Vector3f white(1, 1, 1);

    // coordinates of the screen
    float screen_dist = 0.5;
    float camera_z = screen_dist + origin.getZ();

    float min_x = -0.5;
    float max_x = 0.5;
    min_x += origin.getX();
    max_x += origin.getX();

    float min_y = -0.5;
    float max_y = 0.5;
    min_y += origin.getY();
    max_y += origin.getY();

    // resolution on each axis
    int nx = 800;
    int ny = 800;

    // color depth
    int depth = 4096;

    // list objects
    vector<Hittable *> objects;

    // Some example spheres
    objects.push_back(new Sphere(Vector3f(1, 2, 3), 0.5, red, 1, 0, 1));
    objects.push_back(new Sphere(Vector3f(-1, 2, 4), 0.5, green, 1, 0.25, 5));
    objects.push_back(new Sphere(Vector3f(-1, 1, 1.5), 0.5, blue, 0.03, 0.5, 10));
    objects.push_back(new Sphere(Vector3f(0, 3, 5), 1, white, 0.7, 1, 100));
    objects.push_back(new Plane(Vector3f(0, 1, 0), Vector3f(0, 0, 0), white, 0.3, 0.1, 10));
    objects.push_back(new Plane(Vector3f(0, 0, -1), Vector3f(0, 0, 6), white, 0.3, 0.5, 15));
    objects.push_back(new Triangle(Vector3f(1, 1, 2), Vector3f(0, 2, 2.5), Vector3f(-1, 1, 1.5), white, 0.3, 0.5, 15));
    int n_objects = objects.size();
    vector<vector<Vector3f>> image(nx, vector<Vector3f>(ny));

// antialiasing loop
#pragma omp parallel for schedule(dynamic)
    for (int iter = 0; iter < n_iter; iter++)
    { // iterate through pixels
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
            {
                // get normalised direction

                float random = (float)rand() / RAND_MAX - 0.5;
                float pixel_x = min_x + (max_x - min_x) / nx * (i + random);
                random = (float)rand() / RAND_MAX - 0.5;
                float pixel_y = min_y + (max_y - min_y) / ny * (j + random);

                float pixel_z = camera_z;

                Vector3f pixel(pixel_x, pixel_y, pixel_z);

                Vector3f direction = (pixel - origin);
                direction = direction.normalized();

                // get min dist along a given normalised direction
                float dist_min = 999;
                int obj_min = -1;

                for (int obj = 0; obj < n_objects; obj++)
                {
                    float current_distance = objects[obj]->distance_from_p(origin, direction);
                    if (current_distance > 0 && current_distance < dist_min)
                    {
                        dist_min = current_distance;
                        obj_min = obj;
                    }
                }

                int coordx = nx - 1 - j;
                int coordy = i;

                // point of contact with an object
                Vector3f point = origin + direction * dist_min;
                Vector3f normal;
                Vector3f color;
                float k;
                float sk;
                float pe;
                if (obj_min == -1)
                {
                    normal = Vector3f(0, 0, 0);
                    color = black;
                    k = 0;
                    sk = 0;
                    pe = 0;
                }
                else
                {
                    normal = objects[obj_min]->get_normal(point);
                    color = objects[obj_min]->get_color();
                    k = objects[obj_min]->get_diffc();
                    sk = objects[obj_min]->get_specc();
                    pe = objects[obj_min]->get_phongexp();
                }

                point = point + normal * 0.001;

                Vector3f total_color = color;
                float coeff = 0;
                coeff += ambient_light;

                for (int cur_source = 0; cur_source < n_sources; cur_source++)
                {
                    // Check if there is a direct path to the source
                    Vector3f light_v = (light_sources[cur_source].get_location() - point).normalized();
                    float source_dist = (light_sources[cur_source].get_location() - point).length();

                    bool shaded = false;

                    for (int obj = 0; obj < n_objects; obj++)
                    {
                        float current_distance = objects[obj]->distance_from_p(point, light_v);
                        if (current_distance > 0 && current_distance < source_dist)
                        {
                            shaded = true;
                            break;
                        }
                    }

                    if (!shaded)
                    {
                        // Cosine law
                        Vector3f light_v = (light_sources[cur_source].get_location() - point).normalized();
                        float cosine = normal.innerProduct(light_v);

                        if (cosine < 0)
                            cosine = 0;

                        // Blinn-Phong shading

                        Vector3f h = (light_v - direction) / 2;
                        h = h.normalized();

                        float bp_coeff = normal.innerProduct(h);
                        if (bp_coeff < 0)
                            bp_coeff = 0;
                        bp_coeff = pow(bp_coeff, pe);

                        float cur_coeff = k * cosine + sk * bp_coeff;

                        // inverse square law
                        float dist_to_source = (light_sources[cur_source].get_location() - point).length();
                        cur_coeff = cur_coeff * light_sources[cur_source].get_intensity() / (dist_to_source * dist_to_source);

                        coeff += cur_coeff;
                    }
                }
                image[coordx][coordy] = image[coordx][coordy] + total_color * coeff;
            }
    }

    ofstream o("image.ppm");
    o << "P3" << ' ' << ny << ' ' << nx << ' ' << depth << endl;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            image[i][j] = image[i][j] / float(n_iter);
            o << int(depth * image[i][j].getX()) << ' ' << int(depth * image[i][j].getY()) << ' ' << int(depth * image[i][j].getZ()) << endl;
        }
    }
}