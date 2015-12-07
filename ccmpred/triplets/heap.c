#include <stdlib.h>
#include <stdio.h>

#include "heap.h"

#define PARENT(i) ((i) - 1) / 2
#define LEFTCHILD(i) ((i) * 2 + 1)
#define RIGHTCHILD(i) ((i) * 2 + 2)

void swap(triplet_t *a, triplet_t *b) {
	triplet_t c = *a;
	*a = *b;
	*b = c;
}


void _siftdown(heap_t *heap, size_t startpos, size_t pos) {
	triplet_t newitem = heap->data[pos];
	while (pos > startpos) {
		size_t parentpos = PARENT(pos);
		triplet_t *parent = &heap->data[parentpos];
		
		if (!triplet_lt(&newitem, parent)) {
			break;
		}

		heap->data[pos] = *parent;
		pos = parentpos;
	}

	heap->data[pos] = newitem;
}

void _siftup(heap_t *heap, size_t pos) {
	size_t startpos = pos;
	triplet_t newitem = heap->data[pos];

	size_t childpos = LEFTCHILD(pos);
	size_t rightpos = RIGHTCHILD(pos);
	while (childpos < heap->length) {

		// set childpos to the index of the smaller child
		if (rightpos < heap->length && !triplet_lt(&heap->data[childpos], &heap->data[rightpos])) {
			childpos = rightpos;
		}

		// bubble up smaller child
		heap->data[pos] = heap->data[childpos];
		pos = childpos;

		childpos = LEFTCHILD(pos);
		rightpos = RIGHTCHILD(pos);
	}

	// we now have an empty leaf at pos -- add the item there
	// and bubble up to final pos by sifting down parents
	heap->data[pos] = newitem;
	_siftdown(heap, startpos, pos);
}


heap_t *heap_init(size_t max_size) {
	triplet_t *data = (triplet_t *)malloc(sizeof(triplet_t) * max_size);
	heap_t *heap = (heap_t *)malloc(sizeof(heap_t));

	heap->data = data;
	heap->max_size = max_size;
	heap->length = 0;

	return heap;
}

void heap_free(heap_t *heap) {
	free(heap->data);
	free(heap);
}


void heap_push(heap_t *heap, triplet_t element) {
	if(heap->length >= heap->max_size) {
		fprintf(stderr, "Cannot add element over heap size!\n");
		exit(1);
	}

	heap->data[heap->length] = element;
	heap->length++;

	_siftdown(heap, 0, heap->length - 1);
}

triplet_t heap_pop(heap_t *heap) {
	if(heap->length < 1) {
		fprintf(stderr, "Cannot remove element in empty heap!\n");
		exit(1);
	}
	
	triplet_t result = heap->data[0];
	// TODO rebuild heap

	return result;
}

triplet_t heap_pushpop(heap_t *heap, triplet_t element) {
	if(triplet_lt(&heap->data[0], &element)) {
		swap(&element, &heap->data[0]);
		_siftup(heap, 0);
	}

	return element;
}

void heap_push_over(heap_t *heap, triplet_t element) {
	if(heap-> length >= heap->max_size) {
		heap_pushpop(heap, element);
	} else {
		heap_push(heap, element);
	}

}
