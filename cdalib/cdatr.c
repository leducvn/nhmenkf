#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <netinet/in.h>
#include <stdint.h>

static void tr24(uint32_t *, uint16_t *, size_t);
static void tr21(uint16_t *, uint16_t *, size_t);
static void tr42(uint16_t *, uint32_t *, size_t);
static void tr12(uint16_t *, uint16_t *, size_t);

/* unpack */
void CDATR_I2I4(char *i4, char *i2, int *n) {
	if (*n > 0) tr24((uint32_t*)i4, (uint16_t*)i2, (size_t)*n);
}
void CDATR_I2R4(char *r4, char *i2, int *n) {
	if (*n > 0) tr24((uint32_t*)r4, (uint16_t*)i2, (size_t)*n);
}
void CDATR_I2C(char *c, char *i2, int *n) {
	if (*n > 0) tr21((uint16_t*)c, (uint16_t*)i2, (size_t)*n);
}
void MOVEC_I2I4(char *i4, int *offi4, char *i2, int *offi2, int *n) {
	if (*n > 0) tr24((uint32_t*)(i4+*offi4-1), (uint16_t*)(i2+*offi2-1), (size_t)(*n/sizeof(uint16_t)));
}
void MOVEC_I2R4(char *r4, int *offr4, char *i2, int *offi2, int *n) {
	if (*n > 0) tr24((uint32_t*)(r4+*offr4-1), (uint16_t*)(i2+*offi2-1), (size_t)(*n/sizeof(uint16_t)));
}
void MOVEC_I2C(char *c, int *offc, char *i2, int *offi2, int *n) {
	if (*n > 0) tr21((uint16_t*)(c+*offc-1), (uint16_t*)(i2+*offi2-1), (size_t)(*n/sizeof(uint16_t)));
}
void MOVEC_I2I1(int *i1, int *offi1, short *i2, int *offi2, int *n) {
	int i;
	if (*offi2 == 1) {
		*i1 = *i2 / 256;
	} else if (*offi2 == 2) {
		*i1 = *i2 % 256;
	}
	for (i = 0; i < 4-*offi1; i++) {
		*i1 *= 256;
	}
}
static void tr24(uint32_t *p, uint16_t *q, size_t n) {
	/*
	p : 4byte (out)
	q : 2byte (in)
	n : size of q
	
	ex.
		1234 5678 : original(big) in disk
		2143 6587 : byte swap for 2byte(little) in memory <- in
		1234 5678 : restore(big)
		4321 8765 : byte swap for 4byte(little) <- out
	
	*/
	int i;
	uint16_t sn[2];
	for (i = 0; i < n/2; i++) {
		sn[0] = htons(*(q + 2*i));
		sn[1] = htons(*(q + 2*i+1));
		*(p+i) = ntohl(*(uint32_t*)sn);
	}
}
static void tr21(uint16_t *p, uint16_t *q, size_t n) {
	/*
	p : 1byte (out)
	q : 2byte (in)
	n : size of q
	*/
	int i;
	for (i = 0; i < n; i++) {
		*(p+i) = htons(*(q+i));
	}
}
/* pack */
void CDATR_I4I2(char *i2, char *i4, int *n) {
	if (*n > 0) tr42((uint16_t*)i2, (uint32_t*)i4, (size_t)*n);
}
void CDATR_R4I2(char *i2, char *r4, int *n) {
	if (*n > 0) tr42((uint16_t*)i2, (uint32_t*)r4, (size_t)*n);
}
void CDATR_CI2(char *i2, char *c, int *n) {
	if (*n > 0) tr12((uint16_t*)i2, (uint16_t*)c, (size_t)*n);
}
void MOVEC_I4I2(char *i2, int *offi2, char *i4, int *offi4, int *n) {
	if (*n > 0) tr42((uint16_t*)(i2+*offi2-1), (uint32_t*)(i4+*offi4-1), (size_t)(*n/sizeof(uint16_t)));
}
void MOVEC_R4I2(char *i2, int *offi2, char *r4, int *offr4, int *n) {
	if (*n > 0) tr42((uint16_t*)(i2+*offi2-1), (uint32_t*)(r4+*offr4-1), (size_t)(*n/sizeof(uint16_t)));
}
void MOVEC_CI2(char *i2, int *offi2, char *c, int *offc, int *n) {
	if (*n > 0) tr12((uint16_t*)(i2+*offi2-1), (uint16_t*)(c+*offc-1), (size_t)(*n/sizeof(uint16_t)));
}
static void tr42(uint16_t *p, uint32_t *q, size_t n) {
	/*
	p : 2byte (out)
	q : 4byte (in)
	n : size of p
	
	ex.
		4321 8765 : 4byte(little) in memory <- in
		1234 5678 : byte swap for 4byte(big)
		2143 6587 : byte swap for 2byte(little) in memory <- out
		1234 5678 : original(big) in disk
	*/
	int i;
	uint16_t sn[2];
	for (i = 0; i < n/2; i++) {
		*(uint32_t*)sn = htonl(*(q + i));
		*(p + 2*i) = ntohs(*sn);
		*(p + 2*i+1) = ntohs(*(sn + 1));
	}
}
static void tr12(uint16_t *p, uint16_t *q, size_t n) {
	/*
	p : 2byte (out)
	q : 1byte (in)
	n : size of p
	*/
	int i;
	for (i = 0; i < n; i++) {
		*(p+i) = ntohs(*(q+i));
	}
}


/* fortran interface */

void cdatr_i2i4_(char *i4, char *i2, int *n) {
  CDATR_I2I4(i4, i2, n);
}
void cdatr_i2i4__(char *i4, char *i2, int *n) {
  CDATR_I2I4(i4, i2, n);
}

void cdatr_i2r4_(char *r4, char *i2, int *n) {
  CDATR_I2R4(r4, i2, n);
}
void cdatr_i2r4__(char *r4, char *i2, int *n) {
  CDATR_I2R4(r4, i2, n);
}

void cdatr_i2c_(char *c, char *i2, int *n) {
  CDATR_I2C(c, i2, n);
}
void cdatr_i2c__(char *c, char *i2, int *n) {
  CDATR_I2C(c, i2, n);
}


void movec_i2i4_(char *i4, int *offi4, char *i2, int *offi2, int *n) {
  MOVEC_I2I4(i4, offi4, i2, offi2, n);
}
void movec_i2i4__(char *i4, int *offi4, char *i2, int *offi2, int *n) {
  MOVEC_I2I4(i4, offi4, i2, offi2, n);
}


void movec_i2r4_(char *r4, int *offr4, char *i2, int *offi2, int *n) {
  MOVEC_I2R4(r4, offr4, i2, offi2, n);
}
void movec_i2r4__(char *r4, int *offr4, char *i2, int *offi2, int *n) {
  MOVEC_I2R4(r4, offr4, i2, offi2, n);
}


void movec_i2c_(char *c, int *offc, char *i2, int *offi2, int *n) {
  MOVEC_I2C(c, offc, i2, offi2, n);
}
void movec_i2c__(char *c, int *offc, char *i2, int *offi2, int *n) {
  MOVEC_I2C(c, offc, i2, offi2, n);
}

void movec_i2i1_(int *i1, int *offi1, short *i2, int *offi2, int *n) {
  MOVEC_I2I1(i1, offi1, i2, offi2, n);
}
void movec_i2i1__(int *i1, int *offi1, short *i2, int *offi2, int *n) {
  MOVEC_I2I1(i1, offi1, i2, offi2, n);
}



void cdatr_i4i2_(char *i2, char *i4, int *n) {
  CDATR_I4I2(i2, i4, n);
}
void cdatr_i4i2__(char *i2, char *i4, int *n) {
  CDATR_I4I2(i2, i4, n);
}

void cdatr_r4i2_(char *i2, char *r4, int *n) {
  CDATR_R4I2(i2, r4, n);
}
void cdatr_r4i2__(char *i2, char *r4, int *n) {
  CDATR_R4I2(i2, r4, n);
}


void cdatr_ci2_(char *i2, char *c, int *n) {
  CDATR_CI2(i2, c, n);
}
void cdatr_ci2__(char *i2, char *c, int *n) {
  CDATR_CI2(i2, c, n);
}

void movec_i4i2_(char *i2, int *offi2, char *i4, int *offi4, int *n) {
  MOVEC_I4I2(i2, offi2, i4, offi4, n);
}
void movec_i4i2__(char *i2, int *offi2, char *i4, int *offi4, int *n) {
  MOVEC_I4I2(i2, offi2, i4, offi4, n);
}

void movec_r4i2_(char *i2, int *offi2, char *r4, int *offr4, int *n) {
  MOVEC_R4I2(i2, offi2, r4, offr4, n);
}
void movec_r4i2__(char *i2, int *offi2, char *r4, int *offr4, int *n) {
  MOVEC_R4I2(i2, offi2, r4, offr4, n);
}

void movec_ci2_(char *i2, int *offi2, char *c, int *offc, int *n) {
  MOVEC_CI2(i2, offi2, c, offc, n);
}
void movec_ci2__(char *i2, int *offi2, char *c, int *offc, int *n) {
  MOVEC_CI2(i2, offi2, c, offc, n);
}

