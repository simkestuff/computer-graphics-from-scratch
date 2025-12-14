#ifndef VEC_H_
#define VEC_H_

#include <math.h>

typedef struct {
    float x, y, z;
} Vec;

float magnitude(Vec v);
Vec sub_points(Vec p, Vec q);
Vec add_point(Vec p, Vec v);
Vec add(Vec v1, Vec v2);
Vec scalar_product(Vec v, float a);
Vec scalar_divide(Vec v, float a);
float dot(Vec v1, Vec v2);
Vec cross(Vec v1, Vec v2);

#endif // VEC_H_


#ifdef VEC_IMPLEMENTATION

float magnitude(Vec v)
{
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vec sub_points(Vec p, Vec q)
{
    Vec v = { p.x - q.x, p.y - q.y, p.z - q.z };
    return v;
}

Vec sub(Vec p, Vec q)
{
    return (Vec){p.x - q.x, p.y - q.y, p.z - q.z};
}

Vec add_point(Vec p, Vec v)
{
    Vec q = { p.x + v.x, p.y + v.y, p.z + v.z };
    return q;
}

Vec add(Vec v1, Vec v2)
{
    Vec r = { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
    return r;
}

Vec scalar_product(Vec v, float a)
{
    Vec r = { v.x * a, v.y * a, v.z * a };
    return r;
}

Vec scalar_divide(Vec v, float a)
{
    return scalar_product(v, 1.0f/a);
}

float dot(Vec v1, Vec v2)
{
    return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

Vec cross(Vec v1, Vec v2)
{
    Vec v = { v1.y * v2.z - v1.z * v2.y,
	      v1.x * v2.z - v1.z * v2.x,
	      v1.x * v2.y - v1.y * v2.x };
    return v;
}

#endif // VEC_IMPLEMENTATION
