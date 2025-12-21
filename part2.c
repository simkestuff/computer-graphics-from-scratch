#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#define VEC_IMPLEMENTATION
#include "vec.h"

#define DEG2RAD(d) ((d) * (M_PI / 180.0f))

#define swap(T,a,b)    \
   do {		       \
        T tmp;         \
        tmp = (a);     \
        (a) = (b);     \
        (b) = tmp;     \
    } while (0)


#define DA_INIT_CAP 256
typedef struct {
    float *items;
    size_t count;
    size_t capacity;
} FunValues;

#define da_reserve(da,capacity_expected)\
do {\
    if ((capacity_expected) >= (da)->capacity) {\
	if ((da)->capacity == 0) {\
	    (da)->capacity = DA_INIT_CAP > (capacity_expected) ? DA_INIT_CAP : (capacity_expected);\
	} else {\
	    while ((da)->capacity < (capacity_expected))\
		(da)->capacity *= 2;\
	}\
	(da)->items = realloc((da)->items, sizeof(*(da)->items) * (da)->capacity);\
	assert((da)->items != NULL && "memory allocation failed");\
    }\
} while (0)

#define da_append(da,item) \
do { \
    da_reserve((da), (da)->count + 1); \
    (da)->items[(da)->count++] = (item); \
} while (0)

#define da_append_many(da,append_items,append_items_count) \
do { \
    da_reserve((da), (da)->count + (append_items_count)); \
    memcpy((da)->items + (da)->count, (append_items), (append_items_count) * sizeof(*(da)->items)); \
    (da)->count += (append_items_count); \
} while (0)

typedef struct {
    unsigned char R, G, B;
} Color;

#define BLACK (Color){.R = 0, .G = 0, .B = 0}
#define WHITE (Color){.R = 255, .G = 255, .B = 255}
#define RED (Color){.R = 255, .G = 0, .B = 0}
#define GREEN (Color){.R = 0, .G = 255, .B = 0}
#define BLUE (Color){.R = 0, .G = 0, .B = 255}
#define PURPLE (Color){.R = 128, .G = 0, .B = 128}
#define YELLOW (Color){.R = 255, .G = 255, .B = 0}
#define CYAN (Color){.R = 0, .G = 255, .B = 255}
#define BACKGROUND_COLOR WHITE

#define CW 600  // canvas widht
#define CH 600  // canvas height
Color FrameBuffer[CH][CW];


typedef struct {
    int x, y;
    float h;  // intensity [0.0 - 1.0]
} Point;

// cx ∈ [-CW/2, CW/2)
// cy ∈ [-CH/2, CH/2)
void put_pixel(int x, int y, Color color)   // paint pixel to color
{
    int col = CW/2 + x;  
    int row = CH/2 -1 - y;

    if (col < 0 || col >= CW || row < 0 || row >= CH) return;
    
    FrameBuffer[row][col] = color;
}

void interpolate(int i0, float d0, int i1, float d1, FunValues *values)
{
    if (i0 == i1) {
	da_append(values, d0);
    } else {
	float a = (d1 - d0) / (i1 - i0);
	float value = d0;
	for (int i = i0; i <= i1; ++i) {
	    da_append(values, value);
	    value += a;
	}
    }
}

void draw_line(Point p0, Point p1, Color color)
{
    FunValues values = {0};
    
    if (abs(p1.x - p0.x) > abs(p1.y - p0.y)) {
	// line is horizontal-ish
	// make sure x0 < x1
	if (p0.x > p1.x) swap(Point, p0, p1);

	interpolate(p0.x, (float)p0.y, p1.x, (float)p1.y, &values);

	for (int x = p0.x; x <= p1.x; ++x) {
	    put_pixel(x, values.items[x - p0.x], color);
	}
    } else {
	// line is vertical-ish
	// make sure y0 < y1
	if (p0.y > p1.y) swap(Point, p0, p1);

	interpolate(p0.y, (float)p0.x, p1.y, (float)p1.x, &values);

	for (int y = p0.y; y <= p1.y; ++y) {
	    put_pixel(values.items[y-p0.y], y, color);
	}
    }
}

void draw_wireframe_triangle(Point p0, Point p1, Point p2, Color color)
{
    draw_line(p0, p1, color);
    draw_line(p1, p2, color);
    draw_line(p2, p0, color);
}

void draw_filled_triangle(Point p0, Point p1, Point p2, Color color)
{
    // sort the points so y0 < y1 < y2
    if (p1.y < p0.y) swap(Point, p0, p1);
    if (p2.y < p0.y) swap(Point, p0, p2);
    if (p2.y < p1.y) swap(Point, p1, p2);

    // compute x coordinates of triangle edges
    FunValues v01 = {0};
    FunValues v02 = {0};
    FunValues v12 = {0};
    interpolate(p0.y, (float)p0.x, p1.y, (float)p1.x, &v01);
    interpolate(p0.y, (float)p0.x, p2.y, (float)p2.x, &v02);
    interpolate(p1.y, (float)p1.x, p2.y, (float)p2.x, &v12);

    // concanate short sides
    FunValues v012 = {0};
    da_append_many(&v012, v01.items, v01.count - 1); // skip last one because overlapping
    da_append_many(&v012, v12.items, v12.count);
    
    // determine which is left and which is right
    int m = v012.count / 2;
    FunValues  *x_left;
    FunValues  *x_right;
    if (v02.items[m] < v012.items[m]) {
	x_left = &v02;
	x_right = &v012;
    } else {
	x_left = &v012;
	x_right = &v02;
    }

    // draw horizontal segments
    for (int y = p0.y; y < p2.y; ++y) {
	for (int x = x_left->items[y - p0.y]; x < x_right->items[y - p0.y]; ++x) {
	    put_pixel(x, y, color);
	}
    }
}

void draw_shaded_triangle(Point p0, Point p1, Point p2, Color color)
{
    // sort the points so y0 < y1 < y2
    if (p1.y < p0.y) swap(Point, p0, p1);
    if (p2.y < p0.y) swap(Point, p0, p2);
    if (p2.y < p1.y) swap(Point, p1, p2);

    // compute x coordinates of triangle edges
    FunValues v01 = {0};
    FunValues v02 = {0};
    FunValues v12 = {0};
    interpolate(p0.y, (float)p0.x, p1.y, (float)p1.x, &v01);
    interpolate(p0.y, (float)p0.x, p2.y, (float)p2.x, &v02);
    interpolate(p1.y, (float)p1.x, p2.y, (float)p2.x, &v12);

    // compute h(intensity) along triangle sides
    FunValues h01 = {0};
    FunValues h02 = {0};
    FunValues h12 = {0};
    interpolate(p0.y, p0.h, p1.y, p1.h, &h01);
    interpolate(p0.y, p0.h, p2.y, p2.h, &h02);
    interpolate(p1.y, p1.h, p2.y, p2.h, &h12);
    
    // concanate short sides
    FunValues v012 = {0};
    da_append_many(&v012, v01.items, v01.count - 1); // skip last one because overlapping
    da_append_many(&v012, v12.items, v12.count);

    FunValues h012 = {0};
    da_append_many(&h012, h01.items, h01.count - 1);
    da_append_many(&h012, h12.items, h12.count);
    
    // determine which is left and which is right
    int m = v012.count / 2;
    FunValues  *x_left;
    FunValues  *x_right;
    FunValues *h_left;
    FunValues *h_right;
    if (v02.items[m] < v012.items[m]) {
	x_left = &v02;
	h_left = &h02;
	x_right = &v012;
	h_right = &h012;
    } else {
	x_left = &v012;
	h_left = &h012;
	x_right = &v02;
	h_right = &h02;
    }

    // draw horizontal segments with color shading
    for (int y = p0.y; y < p2.y; ++y) {
	int x_l = x_left->items[y - p0.y];
	int x_r = x_right->items[y - p0.y];
	FunValues h_segment = {0};
	interpolate(x_l, h_left->items[y - p0.y], x_r, h_right->items[y - p0.y], &h_segment);
	for (int x = x_l; x < x_r; ++x) {
	    Color shaded_color = {.R = color.R * h_segment.items[x - x_l],
				  .G = color.G * h_segment.items[x - x_l],
				  .B = color.B * h_segment.items[x - x_l]};
	    put_pixel(x, y, shaded_color);
	}
    }
}

void paint_background(Color color)
{
    for (ssize_t cy = -CH/2; cy < CH/2; ++cy) {
	for (ssize_t cx = -CW/2; cx < CW/2; ++cx) {
	    put_pixel(cx, cy, color);
	}
    }    
}

void make_pic()
{
    FILE *fp = fopen("pic.ppm", "w");

    fprintf(fp, "P6\n%d %d\n255\n", CW, CH);
    fwrite(FrameBuffer, sizeof(FrameBuffer), 1, fp);
    
 fclose(fp);}

// Vertex je točka u prostoru, Vec je točka u matematičkom prostoru
typedef struct {
    Vec pos;
    float h; // intensity
} Vertex;

#define VIEWPORT_WIDTH 1.0f
#define VIEWPORT_HEIGHT 1.0f
#define DISTANCE 1.0f

Vertex make_vertex(float x, float y, float z, float intensity)
{
    Vec v = { .x = x, .y = y, .z = z };
    return (Vertex) { .pos = v, .h = intensity };
}

typedef struct {
    Vertex *items;
    size_t count;
    size_t capacity;
} Vertices;

typedef struct {
    size_t indices[3];
    Color color;
} Triangle;

Triangle mk_triangle(size_t a, size_t b, size_t c, Color color)
{
    Triangle t = {0};
    t.indices[0] = a;
    t.indices[1] = b;
    t.indices[2] = c;
    t.color = color;
    return t;
}

typedef struct {
    Triangle *items;
    size_t count;
    size_t capacity;
} Triangles;


typedef struct {
    Point *items;
    size_t count;
    size_t capacity;
} ProjectedPoints;


typedef enum {
    UNDEF,
    CUBE
} ModelType;

typedef struct {
    ModelType type;
    Vertices *vertices;
    Triangles *triangles;
} Model;

typedef struct {
    Vertex translation; // placement in space
    float scale;
    float rotation; // in degrees around os y
} Transform;

typedef struct {
    Model *model;
    Transform t;
} Instance;


typedef struct {
    Instance *items;
    size_t count;
    size_t capacity;
} Scene;

void render_triangle(Triangle t, ProjectedPoints *pp)
{
    draw_wireframe_triangle(pp->items[t.indices[0]], pp->items[t.indices[1]], pp->items[t.indices[2]], t.color);
}

// NDC means Normalized Device Coordinates.
// It’s the coordinate space after projection and perspective divide, and before mapping to screen pixels.
Point ndc_to_canvas(Vec4 v)
{
    return (Point){
        .x = (int)(v.x * (CW / 2)),
        .y = (int)(v.y * (CH / 2)),
        .h = 1.0f
    };
}


void render_model(Model *model, Mat4 transform)
{
    ProjectedPoints projected = {0};
    for (size_t i = 0; i < (model->vertices)->count; ++i) {
	Vertex vx = (model->vertices)->items[i];
	
	Vec4 vec4 = h_point(vx.pos.x, vx.pos.y, vx.pos.z);
	vec4 = mat4_mul_vec4(transform, vec4);
	vec4 = perspective_divide(vec4);
	
	vx.pos.x = vec4.x;
	vx.pos.y = vec4.y;
	vx.pos.z = vec4.z;
	
	da_append(&projected, ndc_to_canvas(vec4));
    }
    for (size_t i = 0; i < (model->triangles)->count; ++i) {
	Triangle t = (model->triangles)->items[i];
	render_triangle(t, &projected);
    }
}

void render_scene(Scene *scene)
{
    // camera part
    Vec camera = { 0, 0, 0 };
    Mat4 m_view = mat4_translation(-camera.x, -camera.y, -camera.z);
    
    for (size_t i = 0; i < scene->count; ++i) {
	Instance ins = scene->items[i];
	Mat4 m_scale = mat4_scale(ins.t.scale,
				  ins.t.scale,
				  ins.t.scale);
	Mat4 m_rotate = mat4_rotate_y(DEG2RAD(ins.t.rotation));
	Mat4 m_translate = mat4_translation(ins.t.translation.pos.x,
					    ins.t.translation.pos.y,
					    ins.t.translation.pos.z);
	Mat4 m_model = mat4_mul(m_translate, mat4_mul(m_rotate, m_scale));

	Mat4 m_perspective = mat4_perspective(DISTANCE);

	Mat4 MVP = mat4_mul(m_perspective, mat4_mul(m_view, m_model));

	render_model(ins.model, MVP);
    }
}

Model make_cube_model(void)
{
    Model cube = {0};
    cube.type = CUBE;

    // cube vertices  
    Vertex A = make_vertex( 1,  1,  1, 1);
    Vertex B = make_vertex(-1,  1,  1, 1);
    Vertex C = make_vertex(-1, -1,  1, 1);
    Vertex D = make_vertex( 1, -1,  1, 1);
    Vertex E = make_vertex( 1,  1, -1, 1);
    Vertex F = make_vertex(-1,  1, -1, 1);
    Vertex G = make_vertex(-1, -1, -1, 1);
    Vertex H = make_vertex( 1, -1, -1, 1);

    Vertices * vertices = malloc(sizeof(*vertices));
    *vertices = (Vertices){0};
    da_append(vertices, A);
    da_append(vertices, B);
    da_append(vertices, C);
    da_append(vertices, D);
    da_append(vertices, E);
    da_append(vertices, F);
    da_append(vertices, G);
    da_append(vertices, H);

    cube.vertices = vertices;

    Triangle t0  = mk_triangle(0, 1, 2, RED);
    Triangle t1  = mk_triangle(0, 2, 3, RED);
    Triangle t2  = mk_triangle(4, 0, 3, GREEN);
    Triangle t3  = mk_triangle(4, 3, 7, GREEN);
    Triangle t4  = mk_triangle(5, 4, 7, BLUE);
    Triangle t5  = mk_triangle(5, 7, 6, BLUE);
    Triangle t6  = mk_triangle(1, 5, 6, YELLOW);
    Triangle t7  = mk_triangle(1, 6, 2, YELLOW);
    Triangle t8  = mk_triangle(4, 5, 1, PURPLE);
    Triangle t9  = mk_triangle(4, 1, 0, PURPLE);
    Triangle t10 = mk_triangle(2, 6, 7, CYAN);
    Triangle t11 = mk_triangle(2, 7, 3, CYAN);
    
    Triangles *triangles = malloc(sizeof(*triangles));
    *triangles = (Triangles){0};
    da_append(triangles, t0);
    da_append(triangles, t1);
    da_append(triangles, t2);
    da_append(triangles, t3);
    da_append(triangles, t4);
    da_append(triangles, t5);
    da_append(triangles, t6);
    da_append(triangles, t7);
    da_append(triangles, t8);
    da_append(triangles, t9);
    da_append(triangles, t10);
    da_append(triangles, t11);

    cube.triangles = triangles;

    return cube;
}

int main(void)
{
    paint_background(WHITE);

    Scene scene = {0};

    Model cube_model = make_cube_model();

    Transform t1 = { .scale = 0.8f, .translation = make_vertex(1,1,7,1), .rotation = 45};
    Transform t2 = { .scale = 1.0f, .translation = make_vertex(-1.8,-1.5,5,1), .rotation = 0};

    Instance cubeA = (Instance) { .model = &cube_model, .t = t1 };
    Instance cubeB = (Instance) { .model = &cube_model, .t = t2 };

    da_append(&scene, cubeA);
    da_append(&scene, cubeB);
    
    render_scene(&scene);
    
    make_pic();
    
    return 0;
}
