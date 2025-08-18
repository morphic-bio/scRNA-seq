#add library for openmp
#add library for math functions

CC=gcc
CFLAGS=-Ofast -Wall $(INCLUDES)
ifdef DEBUG
	CFLAGS=-ggdb -Wunused-function -fsanitize=thread $(INCLUDES)
	#CFLAGS=-ggdb -Wunused-function -fsanitize=address $(INCLUDES)
endif
LDFLAGS=-lz -lglib-2.0 -fopenmp -lm
#INCLUDES=-I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include -I/usr/lib64/glib-2.0/include -I-I/usr/lib/x86_64-linux-gnu/glib-2.0/include 
INCLUDES=-I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include -I/usr/lib64/glib-2.0/include -I/usr/lib/x86_64-linux-gnu/glib-2.0/include 

default: clean assignBarcodes
all: clean filterfastq assignBarcodes

assignBarcodes: assignBarcodes.c queue.c
	$(CC) $(CFLAGS) -o assignBarcodes assignBarcodes.c queue.c $(LDFLAGS) $(INCLUDES)
filterfastq: filterfastq.c
	$(CC) $(CFLAGS) -o filterfastq filterfastq.c $(LDFLAGS) $(INCLUDES)
clean:
	rm -f filterfastq assignBarcodes

