#ifndef TRIPLETS_PLL_H
#define TRIPLETS_PLL_H

#define N_ALPHA 21

#define x1_index(i,a) (i) * N_ALPHA + (a)
#define X1(i,a) x1[x1_index(i,a)]
#define G1(i,a) g1[x1_index(i,a)]

#define x2_index(ij,a,b) ((ij) * N_ALPHA + (a)) * N_ALPHA + (b)
#define X2(ij,a,b) x2[x2_index(ij,a,b)]
#define G2(ij,a,b) g2[x2_index(ij,a,b)]

#define x3_index(t,a,b,c) (((t) * N_ALPHA + (a)) * N_ALPHA + (b)) * N_ALPHA + (c)
#define X3(b,k,a,j) x3[x3_index(b,k,a,j)]
#define G3(b,k,a,j) g3[x3_index(b,k,a,j)]

#define msa_index(n,i) (n) * ncol + i
#define X(n,i) msa[msa_index(n,i)]

#define SUMPOT(i,a) sum_pot[(i) * N_ALPHA + (a)]
#define PCOND(i,a) p_cond[(i) * N_ALPHA + (a)]


#endif

