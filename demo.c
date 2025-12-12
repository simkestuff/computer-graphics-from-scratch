#include <stdio.h>
#include <math.h>
#include <stdint.h>

#define VEC_IMPLEMENTATION
#include "vec.h"

//////////////////////////////////////////////////////////
// canvas - basically same thing as screen BUT instead of
// origin at top left corner canvas have one at middle
// horizontally and vertically.
//////////////////////////////////////////////////////////

// screen dimensions
#define WIDTH 600
#define HEIGHT 600
#define SPHERENUM 4
#define LIGHTNUM 3

// viewport dimension
const float viewport_width = 1.0f;
const float viewport_height = 1.0f;
const float distance = 1.0f; // projection distance

typedef struct {
    unsigned char r, g, b;
} Color;

typedef struct {
    Vec center;
    float radius;
    Color color;
    uint32_t specular; // how much shiny sufrace is
    float reflective;  // how much reflective surface is
} Sphere;

typedef enum {
    AMBIENT,
    POINT,
    DIRECTIONAL
} LightType;

typedef struct {
    LightType type;
    Vec *vec; // position or direction
    float intensity;
} Light;

Color FrameBuffer[HEIGHT][WIDTH];

Sphere spheres[SPHERENUM];

Light lights[LIGHTNUM];

void putPixel(int x, int y, Color color); // takes screen coord's!
Vec canvas_to_viewpoint(int cx, int cy);
Color TraceRay(Vec o, Vec vp, float tstart, float tstop, int rec_depth);
void intersect_ray_sphere(Vec o, Vec vp, Sphere *s, float *t1, float *t2);
Sphere *closest_intersection(Vec o,  Vec d, float tstart, float tstop, float *closest);
Vec reflect_ray(Vec ray, Vec normal);
float compute_lighting(Vec point, Vec normal, Vec camera, uint32_t specular);
Color color_lighting(Color c, float l);
float clamp_color(float value);
Color color_add(Color first, Color second);

void make_pic()
{
    FILE *fp = fopen("pic.ppm", "w");

    fprintf(fp, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
    fwrite(FrameBuffer, sizeof(FrameBuffer), 1, fp);
    
    fclose(fp);
}

int main(void)
{
    Sphere sr = { (Vec){0, -1, 3},    1.0f,    (Color){255, 0, 0},   500,  0.2f}; 
    Sphere sb = { (Vec){2, 0, 4},     1.0f,    (Color){0, 0, 255},   500,  0.3f}; 
    Sphere sg = { (Vec){-2, 0, 4},    1.0f,    (Color){0, 255, 0},   10,   0.4f}; 
    Sphere sy = { (Vec){0, -5001, 0}, 5000.0f, (Color){255, 255, 0}, 1000, 0.5f}; 
    
    spheres[0] = sr;
    spheres[1] = sb;
    spheres[2] = sg;
    spheres[3] = sy;

    Vec pointlight = {2.0f, 1.0f, 0.0f};
    Vec directionallight = {1.0f, 4.0f, 4.0f};
    Light light_a = { .type = AMBIENT, .intensity = 0.2f, .vec = NULL };
    Light light_p = { .type = POINT, .intensity = 0.6f, .vec = &pointlight };
    Light light_d = { .type = DIRECTIONAL, .intensity = 0.2f, .vec = &directionallight };
    lights[0] = light_a;
    lights[1] = light_p;
    lights[2] = light_d;

    Vec origin = { 0, 0, 0 };

    for (size_t i = 0; i < HEIGHT; ++i) {
	for (size_t j = 0; j < WIDTH; ++j) {

	    // convert pixel to canvas coord
	    int cx = j - (WIDTH / 2); 
	    int cy = (HEIGHT / 2) - i;

	    // map canvas coord onto viewport
	    Vec p = canvas_to_viewpoint(cx, cy);
	    p = scalar_divide(p, magnitude(p));

	    // get color out in the world
	    Color color = TraceRay(origin, p, 1.0, INFINITY, 4);
	    putPixel(j, i, color);
	}
    }
    make_pic();
   
    return 0;
}

Vec canvas_to_viewpoint(int cx, int cy)
{
    Vec p = { (viewport_width/WIDTH) * cx, (viewport_height/HEIGHT) * cy, distance };
    return p;
}

void putPixel(int x, int y, Color color)
{
    FrameBuffer[y][x] = color;
}

Sphere *closest_intersection(Vec o,  Vec d, float tstart, float tstop, float *res)
{
    float closest = INFINITY;
    Sphere *closest_sphere = NULL;
    float t1, t2;
    
    for (size_t i = 0; i < SPHERENUM; ++i) {
	
	intersect_ray_sphere(o, d, &spheres[i], &t1, &t2);
	
	if (t1 > tstart && t1 < tstop && t1 < closest) {
	    closest = t1;
	    closest_sphere = &spheres[i];
	}
	if (t2 > tstart && t2 < tstop && t2 < closest) {
	    closest = t2;
	    closest_sphere = &spheres[i];
	}
    }

    *res = closest;
    return closest_sphere;
}

Color TraceRay(Vec o, Vec vp, float tstart, float tstop, int rec_depth)
{
    float closest = 0.0f;
    Sphere *closest_sphere = closest_intersection(o, vp, tstart, tstop, &closest);
    if (closest_sphere == NULL)
	return (Color){255, 255, 255};

    // compute local color
    Vec p = add_point(o, scalar_product(vp, closest));
    Vec n = sub_points(p, closest_sphere->center);
    n = scalar_divide(n,  magnitude(n));
    
    Color local_color = color_lighting(closest_sphere->color, compute_lighting(p, n, scalar_product(vp, -1), closest_sphere->specular));

    // if we hit recursion limit or the object is not reflective we are done
    float reflective = closest_sphere->reflective;
    if (rec_depth <= 0 || reflective <= 0)
	return local_color;

    // compute reflected color
    Vec ref_ray = reflect_ray(scalar_product(vp, -1), n);
    ref_ray = scalar_divide(ref_ray, magnitude(ref_ray));
    Color reflected_color = TraceRay(p, ref_ray, 0.001f, INFINITY, rec_depth - 1);

    //return color_add(color_lighting(local_color, (1 - reflective)), color_lighting(reflected_color, reflective));
    return color_add(color_lighting(local_color, (1 - reflective)), color_lighting(reflected_color, reflective));
}

Vec reflect_ray(Vec ray, Vec normal)
{
    float dot_n_r = dot(normal, ray);
    Vec scaled = scalar_product(normal, dot_n_r * 2);
    return sub_points(scaled, ray);
}

void intersect_ray_sphere(Vec o, Vec vp, Sphere *s,  float *t1, float *t2)
{
    float radius = s->radius;
    Vec co = sub_points(o, s->center);
    //Vec ray = sub_points(vp, o);
    Vec ray = vp;

    float a = dot(ray, ray);
    float b = 2 * dot(co, ray);
    float c = dot(co, co) - radius*radius;

    float discriminant = b*b - 4*a*c;
    if (discriminant < 0) {
	*t1 = INFINITY;
	*t2 = INFINITY;
    } else {
	*t1 = (-b + sqrtf(discriminant)) / (2*a);
	*t2 = (-b - sqrtf(discriminant)) / (2*a);
    }
}

float compute_lighting(Vec point, Vec normal, Vec camera, uint32_t specular)
{
    camera = scalar_divide(camera, magnitude(camera));
    float li = 0.0f;
    float tmax = INFINITY;

    for (size_t i = 0; i < LIGHTNUM; ++i) {
	Vec direction;
	if (lights[i].type == AMBIENT)
	    li += lights[i].intensity;
	else {
	    if (lights[i].type == POINT) {
		direction = sub_points(*lights[i].vec, point);
		tmax = 1;
	    } else {
		direction = *lights[i].vec;
		tmax = INFINITY;
	    }

	    direction = scalar_divide(direction, magnitude(direction));
	    
	    float closest = 0.0f;
	    Sphere *shadow_sphere = closest_intersection(point, direction, 0.001f, tmax, &closest);
	    if (shadow_sphere != NULL) {
		continue;
	    }

	    // difuse  
	    float n_dot_l = dot(normal, direction);
	    if (n_dot_l > 0) // if <0 -> illuminate back side of surface
		li += lights[i].intensity * n_dot_l / (magnitude(normal) * magnitude(direction));

	    // specular
	    if (specular != UINT32_MAX) { // if not mate object
		Vec reflect = reflect_ray(direction, normal);
		reflect = scalar_divide(reflect, magnitude(reflect));
		float r_dot_v = dot(reflect, camera);
		if (r_dot_v > 0)
		    li += lights[i].intensity * pow(r_dot_v / (magnitude(reflect) * magnitude(camera)), specular);
	    }
	}
    }

    return li;
}

Color color_lighting(Color c, float l)
{
    float r = clamp_color(c.r * l);
    float g = clamp_color(c.g * l);
    float b = clamp_color(c.b * l);
    Color res = (Color) {r, g, b};
    return res;
}

float clamp_color(float value)
{
    if (value < 0.0f)
	return 0.0f;
    if (value > 255.0f)
	return 255.0f;
    return value;
}

Color color_add(Color first, Color second)
{
    float r = clamp_color(first.r + second.r);
    float g = clamp_color(first.g + second.g);
    float b = clamp_color(first.b + second.b);

    return (Color){r, g, b};
}
