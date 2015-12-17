#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include "triplet.h"

double evaluate_triplet_pll(
	const double *x,
	double *g,
	const double *weights,
	const uint8_t *msa,
	const uint32_t *triplets,
	const uint32_t nrow,
	const uint32_t ncol,
	const uint32_t ntriplets
) {
	// partition x and g into pointers x1, x2, x3 and g1, g2, g3
	const uint32_t x2pos = ncol * N_ALPHA;
	const uint32_t x3pos = x2pos + ncol * (ncol - 1) / 2 * N_ALPHA * N_ALPHA;
	const uint64_t nvar = x3pos + ntriplets;
	const double *x1 = x;
	const double *x2 = &x[x2pos];
	const double *x3 = &x[x3pos];
	double *g1 = g;
	double *g2 = &g[x2pos];
	double *g3 = &g[x3pos];

	// set fx and gradient to 0 initially
	double fx = 0.0;
	memset(g, 0, sizeof(double) * nvar);

	double *sum_pot = malloc(sizeof(double) * ncol * N_ALPHA);
	double *log_z   = malloc(sizeof(double) * ncol);
	double *p_cond  = malloc(sizeof(double) * ncol * N_ALPHA);

	for(uint32_t n = 0; n < nrow; n++) {
		double weight = weights[n];

		// compute sum_pot
		memset(sum_pot, 0, sizeof(double) * ncol * N_ALPHA);
		for(uint32_t i = 0; i < ncol; i++) {
			for(uint32_t a = 0; a < N_ALPHA - 1; a++) {
				SUMPOT(i, a) = X1(i, a);
			}
		}

		for(uint32_t i = 0, ij = 0; i < ncol; i++) {
			for(uint32_t j = i + 1; j < ncol; j++, ij++) {
				for(uint32_t a = 0; a < N_ALPHA - 1; a++) {
					SUMPOT(i, a) += X2(ij, a, X(n, j));
					SUMPOT(j, a) += X2(ij, X(n, i), a);
				}
			}
		}

		for(uint32_t t = 0; t < ntriplets; t++) {
			uint32_t i = triplets[t * 6];
			uint32_t j = triplets[t * 6 + 1];
			uint32_t k = triplets[t * 6 + 2];
			uint8_t a = triplets[t * 6 + 3];
			uint8_t b = triplets[t * 6 + 4];
			uint8_t c = triplets[t * 6 + 5];

			if(X(n, j) == b && X(n, k) == c) {
				SUMPOT(i, a) += X3(t);
			}

			if(X(n, i) == a && X(n, k) == c) {
				SUMPOT(j, b) += X3(t);
			}

			if(X(n, i) == a && X(n, j) == b) {
				SUMPOT(k, c) += X3(t);
			}
		}

		for(uint32_t i = 0; i < ncol; i++) {
			SUMPOT(i, N_ALPHA - 1) = 0;
		}

		// compute log_z
		memset(log_z, 0, sizeof(double) * ncol);
		for(uint32_t i = 0; i < ncol; i++) {
			for(uint32_t a = 0; a < N_ALPHA - 1; a++) {
				log_z[i] += exp(SUMPOT(i, a));
			}
			log_z[i] = log(log_z[i]);
		}

		// compute p_cond
		memset(p_cond, 0, sizeof(double) * ncol * N_ALPHA);
		for(uint32_t a = 0; a < N_ALPHA - 1; a++) {
			for(uint32_t i = 0; i < ncol; i++) {
				PCOND(i, a) = exp(SUMPOT(i, a) - log_z[i]);
			}
		}

		// compute fx
		for(uint32_t i = 0; i < ncol; i++) {
			fx -= weight * (SUMPOT(i, X(n, i)) - log_z[i]);

			if(X(n, i) == N_ALPHA - 1) {

				fx -= weight * log_z[i];

				for(int a = 0; a < N_ALPHA; a++) {
					PCOND(i, a) = 0;
				}

			} else {
				G1(i, X(n, i)) -= weight;
			}
		}

		
		// compute g1
		for(uint32_t i = 0; i < ncol; i++) {
			for(uint8_t a = 0; a < N_ALPHA - 1; a++) {
				G1(i, a) += weight * PCOND(i, a);
			}
		}

		// compute g2
		for(uint32_t i = 0, ij = 0; i < ncol; i++) {
			for(uint32_t j = i + 1; j < ncol; j++, ij++) {
				uint8_t xni = X(n, i);
				uint8_t xnj = X(n, j);

				G2(ij, xni, xnj) -= weight * 2;

				for(uint8_t a = 0; a < N_ALPHA - 1; a++) {
					G2(ij, a, xnj) += weight * PCOND(i, a);
					G2(ij, xni, a) += weight * PCOND(j, a);
				}
			}
		}

		// compute g3
		for(uint32_t t = 0; t < ntriplets; t++) {
			uint32_t i = triplets[t * 6];
			uint32_t j = triplets[t * 6 + 1];
			uint32_t k = triplets[t * 6 + 2];
			uint8_t a = triplets[t * 6 + 3];
			uint8_t b = triplets[t * 6 + 4];
			uint8_t c = triplets[t * 6 + 5];

			if(X(n, j) == b && X(n, k) == c) {
				G3(t) += weight * PCOND(i, a);
			}

			if(X(n, i) == a && X(n, k) == c) {
				G3(t) += weight * PCOND(j, b);
			}

			if(X(n, i) == a && X(n, j) == b) {
				G3(t) += weight * PCOND(k, c);
			}

			if(X(n, i) == a && X(n, j) == b && X(n, k) == c) {
				G3(t) -= weight * 3;
			}
		}

	}

	// zero out possible gap gradients in G1
	for(uint32_t i = 0; i < ncol; i++) {
		G1(i, N_ALPHA - 1) = 0;
	}

	// zero out possible gap gradients in G2
	uint32_t nij = ncol * (ncol - 1) / 2;
	for(uint32_t ij = 0; ij < nij; ij++) {
		for(uint8_t a = 0; a < N_ALPHA; a++) {
			G2(ij, a, N_ALPHA - 1) = 0;
			G2(ij, N_ALPHA - 1, a) = 0;
		}
	}

	free(sum_pot);
	free(log_z);
	free(p_cond);

	return fx;
}
