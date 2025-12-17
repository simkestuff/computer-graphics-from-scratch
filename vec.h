#ifndef VEC_H_
#define VEC_H_

#include <math.h>

typedef struct {
    float x, y, z;
} Vec;

typedef struct {
    float x, y, z, w;
} Vec4;

typedef struct {
    float m[4][4];
} Mat4;

float magnitude(Vec v);
Vec sub_points(Vec p, Vec q);
Vec add_point(Vec p, Vec v);
Vec add(Vec v1, Vec v2);
Vec scalar_product(Vec v, float a);
Vec scalar_divide(Vec v, float a);
float dot(Vec v1, Vec v2);
Vec cross(Vec v1, Vec v2);
Vec4 direction(float x, float y, float z);
Vec4 h_point(float x, float y, float z);
Vec4 mat4_mul_vec4(Mat4 m, Vec4 v);
Mat4 mat4_mul(Mat4 a, Mat4 b);
Mat4 mat4_identity(void);
Mat4 mat4_translation(float tx, float ty, float tz);
Mat4 mat4_rotate_y(float angle);
Mat4 mat4_scale(float sx, float sy, float sz);
Mat4 mat4_perspective(float d);
Vec4 perspective_divide(Vec4 v);

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

Vec4 direction(float x, float y, float z)
{
    return (Vec4){ .x = x, .y = y, .z = z, .w = 0.0f };
}

Vec4 h_point(float x, float y, float z)
{
    return (Vec4){ .x = x, .y = y, .z = z, .w = 1.0f };
}

Vec4 mat4_mul_vec4(Mat4 m, Vec4 v)
{
    return (Vec4){
        m.m[0][0]*v.x + m.m[0][1]*v.y + m.m[0][2]*v.z + m.m[0][3]*v.w,
        m.m[1][0]*v.x + m.m[1][1]*v.y + m.m[1][2]*v.z + m.m[1][3]*v.w,
        m.m[2][0]*v.x + m.m[2][1]*v.y + m.m[2][2]*v.z + m.m[2][3]*v.w,
        m.m[3][0]*v.x + m.m[3][1]*v.y + m.m[3][2]*v.z + m.m[3][3]*v.w
    };
}

Mat4 mat4_mul(Mat4 a, Mat4 b)
{
    Mat4 r = {0};

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                r.m[i][j] += a.m[i][k] * b.m[k][j];
            }
        }
    }

    return r;
}

Mat4 mat4_identity(void)
{
    Mat4 m = {0};
    m.m[0][0] = 1;
    m.m[1][1] = 1;
    m.m[2][2] = 1;
    m.m[3][3] = 1;
    return m;
}

Mat4 mat4_translation(float tx, float ty, float tz)
{
    Mat4 m = mat4_identity();
    m.m[0][3] = tx;
    m.m[1][3] = ty;
    m.m[2][3] = tz;
    return m;
}

Mat4 mat4_rotate_y(float angle)
{
    float c = cosf(angle);
    float s = sinf(angle);

    Mat4 m = mat4_identity();
    m.m[0][0] =  c;
    m.m[0][2] =  s;
    m.m[2][0] = -s;
    m.m[2][2] =  c;

    return m;
}

Mat4 mat4_scale(float sx, float sy, float sz)
{
    Mat4 m = {0};
    m.m[0][0] = sx;
    m.m[1][1] = sy;
    m.m[2][2] = sz;
    m.m[3][3] = 1;
    return m;
}

Mat4 mat4_perspective(float d)
{
    Mat4 m = {0};
    m.m[0][0] = 1;
    m.m[1][1] = 1;
    m.m[2][2] = 1;
    m.m[3][2] = 1.0f / d;
    return m;
}

Vec4 perspective_divide(Vec4 v)
{
    return (Vec4){
        v.x / v.w,
        v.y / v.w,
        v.z / v.w,
        1.0f
    };
}

#endif // VEC_IMPLEMENTATION
