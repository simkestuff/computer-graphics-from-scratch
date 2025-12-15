#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define swap(T,a,b)    \
   do {		       \
        T tmp;         \
        tmp = (a);     \
        (a) = (b);     \
        (b) = tmp;     \
    } while (0)

#define CW 600  // canvas widht
#define CH 600  // canvas height

#define FUN_VALUES_INITIAL_CAPACITY 256
typedef struct {
    float *items;
    size_t count;
    size_t capacity;
} FunValues;

void values_capacity_reserve(FunValues *fv, size_t capacity_expected)
{
    if (capacity_expected >= fv->capacity) {
	if (fv->capacity == 0) {
	    fv->capacity = FUN_VALUES_INITIAL_CAPACITY > capacity_expected ? FUN_VALUES_INITIAL_CAPACITY : capacity_expected;
	} else {
	    while (fv->capacity < capacity_expected)
		fv->capacity *= 2;
	}

	fv->items = realloc(fv->items, sizeof(*fv->items) * fv->capacity);

	assert(fv->items != NULL && "memory allocation failed");
    }

}

void values_append(FunValues *fv, float item)
{
    values_capacity_reserve(fv, fv->count + 1);

    fv->items[fv->count++] = item;
}

void values_append_many(FunValues *fv, float *items, size_t items_count)
{
    values_capacity_reserve(fv, fv->count + items_count);
    memcpy(fv->items + fv->count, items, items_count * sizeof(*fv->items));

    fv->count += items_count;
}

typedef struct {
    unsigned char R, G, B;
} Color;

#define BACKGROUND_COLOR (Color){.R = 255, .B = 255, .G = 255}

// this type of initialization needs gcc/clang extension
// use paint_background function as alternative
Color FrameBuffer[CH][CW] = {
    [0 ... CH-1][0 ... CW-1] = BACKGROUND_COLOR
};


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

    FrameBuffer[row][col] = color;
}

void interpolate(int i0, float d0, int i1, float d1, FunValues *values)
{
    if (i0 == i1) {
	values_append(values, d0);
    } else {
	float a = (d1 - d0) / (i1 - i0);
	float value = d0;
	for (int i = i0; i <= i1; ++i) {
	    values_append(values, value);
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
    values_append_many(&v012, v01.items, v01.count - 1); // skip last one because overlapping
    values_append_many(&v012, v12.items, v12.count);
    
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
    values_append_many(&v012, v01.items, v01.count - 1); // skip last one because overlapping
    values_append_many(&v012, v12.items, v12.count);

    FunValues h012 = {0};
    values_append_many(&h012, h01.items, h01.count - 1);
    values_append_many(&h012, h12.items, h12.count);
    
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

int main(void)
{
    Color color = {0};
    Color c = {0,255,0};
    draw_wireframe_triangle((Point){-200,-250,0.3}, (Point){200,50,0.2}, (Point){20,250,1.0}, color);
    draw_shaded_triangle((Point){-200,-250,0.4}, (Point){200,50,0.2}, (Point){20,250,1.0}, c);
    
    make_pic();
    
    return 0;
}
