#include <stdio.h>
#include <math.h>

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
Color TraceRay(Vec o, Vec vp, float tstart, float tstop);
void intersect_ray_sphere(Vec o, Vec vp, Sphere *s, float *t1, float *t2);

float compute_lighting(Vec point, Vec normal);
Color color_lighting(Color c, float l);

void make_pic()
{
    FILE *fp = fopen("pic.ppm", "w");

    fprintf(fp, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
    fwrite(FrameBuffer, sizeof(FrameBuffer), 1, fp);
    
    fclose(fp);
}

int main(void)
{
    Sphere sr = { (Vec){0, -1, 3}, 1.0f, (Color){255, 0, 0} };
    Sphere sb = { (Vec){2, 0, 4}, 1.0f, (Color){0, 0, 255} };
    Sphere sg = { (Vec){-2, 0, 4}, 1.0f, (Color){0, 255, 0} };
    Sphere sy = { (Vec){0, -5001, 0}, 5000.0f, (Color){255, 255, 0} };
    

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

	    // get color out in the world
	    Color color = TraceRay(origin, p, 1.0, INFINITY);
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

Color TraceRay(Vec o, Vec vp, float tstart, float tstop)
{
    float closest = INFINITY;
    Sphere *closest_sphere = NULL;
    float t1, t2;
    
    for (size_t i = 0; i < SPHERENUM; ++i) {
	
	intersect_ray_sphere(o, vp, &spheres[i], &t1, &t2);
	
	if (t1 > tstart && t1 < tstop && t1 < closest) {
	    closest = t1;
	    closest_sphere = &spheres[i];
	}
	if (t2 > tstart && t2 < tstop && t2 < closest) {
	    closest = t2;
	    closest_sphere = &spheres[i];
	}
    }
    
    if (closest_sphere == NULL)
	return (Color){255, 255, 255};

    Vec p = add_point(o, scalar_product(vp, closest));
    Vec n = sub_points(p, closest_sphere->center);
    n = scalar_divide(n,  magnitude(n));
    
    return color_lighting(closest_sphere->color, compute_lighting(p, n));
}

void intersect_ray_sphere(Vec o, Vec vp, Sphere *s,  float *t1, float *t2)
{
    float radius = s->radius;
    Vec co = sub_points(o, s->center);
    Vec ray = sub_points(vp, o);

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

float compute_lighting(Vec point, Vec normal)
{
    float li = 0.0f;

    for (size_t i = 0; i < LIGHTNUM; ++i) {
	Vec direction;
	if (lights[i].type == AMBIENT)
	    li += lights[i].intensity;
	else {
	    if (lights[i].type == POINT) {
		direction = sub_points(*lights[i].vec, point);
	    } else {
		direction = *lights[i].vec;
	    }

	    float n_dot_l = dot(normal, direction);
	    if (n_dot_l > 0) // if <0 -> illuminate back side of surface
		li += lights[i].intensity * n_dot_l / (magnitude(normal) * magnitude(direction));
	}
    }

    return li;
}

Color color_lighting(Color c, float l)
{
    Color res = (Color) {c.r * l, c.g * l, c.b * l};
    return res;
}
