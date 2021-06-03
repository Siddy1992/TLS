#ifndef NTT_H
#define NTT_H

#include <stdint.h>

#define N 256

//for DILITHIUM 1
#define KYBER_DILITHIUM 0


#if KYBER_DILITHIUM==0 //KYBER
    #define LEN 2
    #define Q 3329
    #define QINV -3327
    #define F 1441  // mont^2/128
    #define INTT_K 127
	  #define INTT_LIMIT 128
    #define SHIFT 16
    typedef int16_t INT;

#elif KYBER_DILITHIUM==1  //DILITHIUM
    #define LEN 1
    #define Q 8380417
    #define QINV 58728449
    #define F 41978  // mont^2/256
    #define INTT_K 255
	  #define INTT_LIMIT 255
    #define SHIFT 32
    typedef int32_t INT;

#endif

extern const int16_t zetas_kyber[128];
extern const int32_t zetas_dilithium[256];
int16_t barrett_reduce(int16_t a);
int32_t montgomery_reduce_64(int64_t a);
int16_t montgomery_reduce(int32_t a);
void ntt(INT a[N]);
void invntt(INT a[N]);
void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta);

#endif
