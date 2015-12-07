#include <stdint.h>
#include <stdlib.h>
#include "heap.h"

#define N_ALPHA 21
#define x2_index(i, j, a, b) ((((i) * ncol + j) * N_ALPHA + a) * N_ALPHA + b)
#define X2(i, j, a, b) x_pair[x2_index(i, j, a, b)]

bool triplet_lt(triplet_t *a, triplet_t *b) {
	return a->score < b->score;
}

int compare_triplets(const void *a, const void *b) {
	const triplet_t *ta = (triplet_t *) a;
	const triplet_t *tb = (triplet_t *) b;

	return (ta->score < tb->score) - (ta->score > tb->score);
}


uint64_t find_triplets(double *x_pair, uint32_t *triplets, uint32_t ncol, uint32_t ntriplets, uint32_t min_separation) {

	heap_t *heap = heap_init(ntriplets);

	for(uint32_t i = 0; i < ncol; i++) {
		for(uint32_t j = i + min_separation; j < ncol; j++) {
			for(uint32_t k = j + min_separation; k < ncol; k++) {

				for(uint8_t a = 0; a < N_ALPHA - 1; a++) {
					for(uint8_t b = 0; b < N_ALPHA - 1; b++) {
						for(uint8_t c = 0; c < N_ALPHA - 1; c++) {

							double score = X2(i, j, a, b) + X2(i, k, a, c) + X2(j, k, b, c);
							triplet_t trp = {i, j, k, a, b, c, score};
							heap_push_over(heap, trp);

						}
					}
				}

			}
		}
	}

	qsort(heap->data, heap->length, sizeof(triplet_t), compare_triplets);

	for(uint32_t t = 0; t < heap->length; t++) {
		triplets[t * 6 + 0] = heap->data[t].i;
		triplets[t * 6 + 1] = heap->data[t].j;
		triplets[t * 6 + 2] = heap->data[t].k;

		triplets[t * 6 + 3] = heap->data[t].a;
		triplets[t * 6 + 4] = heap->data[t].b;
		triplets[t * 6 + 5] = heap->data[t].c;
	}

	uint64_t len = heap->length;

	heap_free(heap);

	return len;

}
