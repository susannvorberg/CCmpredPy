#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "heap.h"

#define PARENT(i) ((i) - 1) / 2
#define LEFTCHILD(i) ((i) * 2 + 1)
#define RIGHTCHILD(i) ((i) * 2 + 2)

void swap(void *a, void *b, size_t size) {
	void *temp = malloc(size);
	memcpy(temp, a, size);
	memcpy(a, b, size);
	memcpy(b, temp, size);
	free(temp);
}


void _siftdown(heap_t *heap, size_t startpos, size_t pos) {
	void *newitem = malloc(heap->element_size);

	// newitem = heap->data[pos]
	memcpy(
		newitem,
		heap->data + (pos * heap->element_size),
		heap->element_size
	);

	while (pos > startpos) {
		size_t parentpos = PARENT(pos);

		// newitem >= parent
		if (heap->comparator(
				newitem,
				heap->data + (parentpos * heap->element_size)
			) <= 0) {
			break;
		}

		// heap->data[pos] = *parent;
		memcpy(
			heap->data + (pos * heap->element_size),
			heap->data + (parentpos * heap->element_size),
			heap->element_size
		);

		pos = parentpos;
	}

	// heap->data[pos] = newitem
	memcpy(
		heap->data + (pos * heap->element_size),
		newitem,
		heap->element_size
	);

	free(newitem);
}

void _siftup(heap_t *heap, size_t pos) {
	void *newitem = malloc(heap->element_size);

	// newitem = heap->data[pos];
	memcpy(
		newitem,
		heap->data + (pos * heap->element_size),
		heap->element_size
	);

	size_t startpos = pos;
	size_t childpos = LEFTCHILD(pos);
	size_t rightpos = RIGHTCHILD(pos);
	while (childpos < heap->length) {

		// set childpos to the index of the smaller child
		// rightpos < heap->length && heap->data[childpos] >= heap->data[rightpos]
		if (rightpos < heap->length &&
			heap->comparator(
				heap->data + (childpos * heap->element_size),
				heap->data + (rightpos * heap->element_size)
			) <= 0) {
			childpos = rightpos;
		}

		// bubble up smaller child
		// heap->data[pos] = heap->data[childpos];
		memcpy(
			heap->data + (pos * heap->element_size),
			heap->data + (childpos * heap->element_size),
			heap->element_size
		);

		pos = childpos;

		childpos = LEFTCHILD(pos);
		rightpos = RIGHTCHILD(pos);
	}

	// we now have an empty leaf at pos -- add the item there
	// and bubble up to final pos by sifting down parents
	// heap->data[pos] = newitem;
	memcpy(
		heap->data + (pos * heap->element_size),
		newitem,
		heap->element_size
	);

	_siftdown(heap, startpos, pos);

	free(newitem);
}


heap_t *heap_init(size_t max_size, size_t element_size, element_comparator comparator) {
	void *data = malloc(element_size * max_size);
	heap_t *heap = (heap_t *)malloc(sizeof(heap_t));

	heap->data = data;
	heap->max_size = max_size;
	heap->element_size = element_size;
	heap->length = 0;
	heap->comparator = comparator;

	return heap;
}

void heap_free(heap_t *heap) {
	free(heap->data);
	free(heap);
}


void heap_push(heap_t *heap, void *element) {
	if(heap->length >= heap->max_size) {
		fprintf(stderr, "Cannot add element over heap size!\n");
		exit(1);
	}

	// heap->data[heap->length] = element;
	memcpy(
		heap->data + (heap->length * heap->element_size),
		element,
		heap->element_size
	);

	heap->length++;

	_siftdown(heap, 0, heap->length - 1);
}

void heap_pop(heap_t *heap, void *out_data) {
	if(heap->length < 1) {
		fprintf(stderr, "Cannot remove element in empty heap!\n");
		exit(1);
	}

	if(out_data != NULL) {
		memcpy(out_data, heap->data, heap->element_size);
	}

	heap->length--;
}

void heap_pushpop(heap_t *heap, void *element) {
	// heap->data[0] < element
	if(heap->comparator(heap->data, element) > 0) {
		swap(element, heap->data, heap->element_size);
		_siftup(heap, 0);
	}
}

void heap_push_over(heap_t *heap, void *element) {
	if(heap->length >= heap->max_size) {
		heap_pushpop(heap, element);
	} else {
		heap_push(heap, element);
	}

}
