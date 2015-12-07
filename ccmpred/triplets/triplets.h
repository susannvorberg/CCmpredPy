#include <stdint.h>
#include <stdbool.h>

#ifndef TRIPLET_H
#define TRIPLET_H

typedef struct triplet6_t {
	uint32_t i;
	uint32_t j;
	uint32_t k;
	uint8_t a;
	uint8_t b;
	uint8_t c;
	double score;
} triplet6_t;

bool triplet6_lt(triplet6_t *a, triplet6_t *b);

typedef struct triplet3_t {
	uint32_t i;
	uint32_t j;
	uint32_t k;
	double score;
} triplet3_t;

bool triplet3_lt(triplet3_t *a, triplet3_t *b);

#endif
