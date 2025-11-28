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

// viewport dimension
const float viewport_width = 1.0f;
const float viewport_height = 1.0f;
const float distance = 1.0f; // projection distance

typedef struct {
    unsigned char r, g, b;
} Color;

typedef struct {
    Point center;
    float radius;
    Color color;
} Sphere;

Color FrameBuffer[HEIGHT][WIDTH];

Sphere spheres[3];

void putPixel(int x, int y, Color color); // takes screen coord's!
Point canvas_to_viewpoint(int cx, int cy);
Color TraceRay(Point o, Point vp, float tstart, float tstop);
void intersect_ray_sphere(Point o, Point vp, Sphere *s, float *t1, float *t2);


void make_pic()
{
    FILE *fp = fopen("pic.ppm", "w");

    fprintf(fp, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
    fwrite(FrameBuffer, sizeof(FrameBuffer), 1, fp);
    
    fclose(fp);
}

int main(void)
{
    Sphere sr = { (Point){0, -1, 3}, 1.0f, (Color){255, 0, 0} };
    Sphere sb = { (Point){2, 0, 4}, 1.0f, (Color){0, 0, 255} };
    Sphere sg = { (Point){-2, 0, 4}, 1.0f, (Color){0, 255, 0} };

    spheres[0] = sr;
    spheres[1] = sb;
    spheres[2] = sg;

    Point origin = { 0, 0, 0 };

    for (size_t i = 0; i < HEIGHT; ++i) {
	for (size_t j = 0; j < WIDTH; ++j) {

	    // convert pixel to canvas coord
	    int cx = j - (WIDTH / 2); 
	    int cy = (HEIGHT / 2) - i;

	    // map canvas coord onto viewport
	    Point p = canvas_to_viewpoint(cx, cy);

	    // get color out in the world
	    Color color = TraceRay(origin, p, 1.0, INFINITY);
	    putPixel(j, i, color);
	}
    }
    make_pic();
   
    return 0;
}

Point canvas_to_viewpoint(int cx, int cy)
{
    Point p = { (viewport_width/WIDTH) * cx, (viewport_height/HEIGHT) * cy, distance };
    return p;
}

void putPixel(int x, int y, Color color)
{
    FrameBuffer[y][x] = color;
}

Color TraceRay(Point o, Point vp, float tstart, float tstop)
{
    float closest = INFINITY;
    Sphere *closest_sphere = NULL;
    float t1, t2;
    
    for (size_t i = 0; i < 3; ++i) {
	
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
    
    return closest_sphere->color;
}

void intersect_ray_sphere(Point o, Point vp, Sphere *s,  float *t1, float *t2)
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
