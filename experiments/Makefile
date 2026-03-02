CC = gcc
CFLAGS = -O3 -march=native -Wall

all: A14_isochrone_mc

A14_isochrone_mc: A14_isochrone_mc.c
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f A14_isochrone_mc
