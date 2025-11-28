CC = gcc
CPPFLAGS =
CFLAGS = -Wall -Wextra -g
LDFLAGS =

.PHONY = all clean

all: demo

demo: demo.c
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f *.o demo
