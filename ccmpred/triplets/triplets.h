#include <stdint.h>
#include <stdbool.h>

#ifndef TRIPLET_H
#define TRIPLET_H

typedef struct triplet_t {
	uint32_t i;
	uint32_t j;
	uint32_t k;
	uint8_t a;
	uint8_t b;
	uint8_t c;
	double score;
} triplet_t;

bool triplet_lt(triplet_t *a, triplet_t *b);

#endif
