#include <stddef.h>

#include "triplets.h"


#ifndef HEAP_H
#define HEAP_H

typedef struct heap_t {
	size_t max_size;
	size_t length;
	triplet_t *data;
} heap_t;

heap_t *heap_init(size_t max_size);
void heap_free(heap_t *heap);

void heap_push(heap_t *heap, triplet_t element);
triplet_t heap_pop(heap_t *heap);
triplet_t heap_pushpop(heap_t *heap, triplet_t element);
void heap_push_over(heap_t *heap, triplet_t element);

#endif
