#include <stddef.h>

#include "triplets.h"


#ifndef HEAP_H
#define HEAP_H

typedef int (*element_comparator)(const void *, const void *);

typedef struct heap_t {
	size_t max_size;
	size_t length;
	size_t element_size;
	element_comparator comparator;
	void *data;

} heap_t;

heap_t *heap_init(size_t max_size, size_t element_size, element_comparator comparator);
void heap_free(heap_t *heap);

void heap_push(heap_t *heap, void *element);
void heap_pop(heap_t *heap, void *out_data);
void heap_pushpop(heap_t *heap, void *element);
void heap_push_over(heap_t *heap, void *element);

#endif
