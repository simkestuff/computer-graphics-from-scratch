CC = gcc
CPPFLAGS =
CFLAGS = -Wall -Wextra -g
LDFLAGS =

.PHONY = all clean

all: demo

demo: demo.c vec.h
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f *.o demo ./pic.ppm

