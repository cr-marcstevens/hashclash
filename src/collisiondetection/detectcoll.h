/**************************************************************************\
|
|    Copyright (C) 2012 Marc Stevens
|    FOR RESEARCH COLLABORATION ONLY, NOT FOR REDISTRIBUTION
|
\**************************************************************************/

#include <stdint.h>

typedef struct {
	uint64_t total;
	uint32_t ihv[4];
	unsigned char buffer[64];
	int bigendian;
	int found_collision;
	int safe_hash;

	uint32_t ihv1[4];
	uint32_t ihv2[4];
	uint32_t m2[16];
	uint32_t states[260];
	uint32_t statesmsb[260];
	uint32_t tmpihv[4];
	uint32_t tmpblock[16];
	uint32_t previhv[4];
	uint32_t prevblock[16];
} MD5_CTX;

void MD5_Init(MD5_CTX*); // outputs MD5 hash if no collision was found and a modified-MD5 hash otherwise
void MD5_Init_unsafe(MD5_CTX*); // always outputs MD5 hash
void MD5_Update(MD5_CTX*, const char* buf, unsigned len);
int  MD5_Finish(unsigned char hash[16], MD5_CTX*); // returns: 0 = no collision, otherwise = collision found => warn user for active attack

void md5compress_states(uint32_t ihv[4], const uint32_t block[16], uint32_t states[260]);
int md5recompress_fast(unsigned t, uint32_t ihv[4], const uint32_t block[16], const uint32_t state[4], const uint32_t rihv[4]);
int detect_coll(const uint32_t block1[16], const uint32_t states[260], const uint32_t statesmsb[260], const uint32_t tihv[4], uint32_t ihv2[4], uint32_t block2[16]);

typedef struct {
	uint32_t msgdiff[16];
	unsigned t;
	int negate;
	int zero;
	int msb;
} msgdiff_tuples_t;
extern msgdiff_tuples_t msgdiff_tuples[];



typedef struct {
	uint64_t total;
	uint32_t ihv[5];
	unsigned char buffer[64];
	int bigendian;
	int found_collision;
	int safe_hash;

	uint32_t ihv1[5];
	uint32_t ihv2[5];
	uint32_t m1[80];
	uint32_t m2[80];
	uint32_t states[81*5];
} SHA1_CTX;	

void SHA1_Init(SHA1_CTX*); // outputs SHA-1 hash if no collision was found and a modified-SHA-1 hash otherwise
void SHA1_Init_unsafe(SHA1_CTX*); // always outputs SHA-1 hash
void SHA1_Update(SHA1_CTX*, const char* buf, unsigned len);
int  SHA1_Finish(unsigned char hash[20], SHA1_CTX*); // returns: 0 = no collision, otherwise = collision found => warn user for active attack

void sha1compress_me(const uint32_t block[16], uint32_t me[80]);
void sha1compress_states(uint32_t ihv[5], const uint32_t me[80], uint32_t states[81*5]);
int sha1recompress_fast(unsigned t, uint32_t ihv[5], const uint32_t me[80], const uint32_t state[5], const uint32_t rihv[5]);

