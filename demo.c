#include <stdio.h>

typedef struct {
    ;
} Canvas;

typedef struct {
    unsigned short r, g, b;
} Color;

void putPixel(Canvas *canvas, int x, int y, Color color);

int main(void)
{
    printf("hello world!\n");
    return 0;
}
