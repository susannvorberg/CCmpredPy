#include <stdint.h>
#include <stdlib.h>
#include "heap.h"

#define N_ALPHA 21
#define x2_index(i, j, a, b) ((((i) * ncol + j) * N_ALPHA + a) * N_ALPHA + b)
#define X2(i, j, a, b) x_pair[x2_index(i, j, a, b)]

int compare_triplet6(const void *a, const void *b) {
	const triplet6_t *ta = (triplet6_t *) a;
	const triplet6_t *tb = (triplet6_t *) b;

	return (ta->score < tb->score) - (ta->score > tb->score);
}


int compare_triplet3(const void *a, const void *b) {
	const triplet3_t *ta = (triplet3_t *) a;
	const triplet3_t *tb = (triplet3_t *) b;

	return (ta->score < tb->score) - (ta->score > tb->score);
}


uint64_t find_triplets(double *x_pair, uint32_t *triplets, uint32_t ncol, uint32_t ntriplets, uint32_t min_separation) {

	heap_t *heap = heap_init(ntriplets, sizeof(triplet6_t), compare_triplet6);

	for(uint32_t i = 0; i < ncol; i++) {
		for(uint32_t j = i + min_separation; j < ncol; j++) {
			for(uint32_t k = j + min_separation; k < ncol; k++) {

				for(uint8_t a = 0; a < N_ALPHA - 1; a++) {
					for(uint8_t b = 0; b < N_ALPHA - 1; b++) {
						for(uint8_t c = 0; c < N_ALPHA - 1; c++) {

							double score = X2(i, j, a, b) + X2(i, k, a, c) + X2(j, k, b, c);
							triplet6_t trp = {i, j, k, a, b, c, score};

							heap_push_over(heap, &trp);

						}
					}
				}

			}
		}
	}

	qsort(heap->data, heap->length, sizeof(triplet6_t), compare_triplet6);

	for(uint32_t t = 0; t < heap->length; t++) {
		triplet6_t *trp = (triplet6_t *)(heap->data + t * sizeof(triplet6_t));

		triplets[t * 6 + 0] = trp->i;
		triplets[t * 6 + 1] = trp->j;
		triplets[t * 6 + 2] = trp->k;

		triplets[t * 6 + 3] = trp->a;
		triplets[t * 6 + 4] = trp->b;
		triplets[t * 6 + 5] = trp->c;
	}

	uint64_t len = heap->length;

	heap_free(heap);

	return len;

}
