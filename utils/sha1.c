/*
SHA-1 in C
By Steve Reid <steve@edmweb.com>
100% Public Domain

   small edits for cado-nfs by E. Thomé. Still PD.

Test Vectors (from FIPS PUB 180-1)
"abc"
  A9993E36 4706816A BA3E2571 7850C26C 9CD0D89D
"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"
  84983E44 1C3BD26E BAAE4AA1 F95129E5 E54670F1
A million repetitions of "a"
  34AA973C D4C4DAA4 F61EEB2B DBAD2731 6534016F
*/

/* #define LITTLE_ENDIAN * This should be #define'd already, if true. */

#include "cado.h" // IWYU pragma: keep
#include <string.h>
#include <stdint.h> /* for uint32_t */

#include "sha1.h"

#define rol(value, bits) (((value) << (bits)) | ((value) >> (32U - (bits))))

/* blk0() and blk() perform the initial expand. */
/* I got the idea of expanding during the round function from SSLeay */
#if BYTE_ORDER == LITTLE_ENDIAN
#define blk0(i) (block->l[i] = (rol(block->l[i],24U)&0xFF00FF00U) \
    |(rol(block->l[i],8U)&0x00FF00FFU))
#elif BYTE_ORDER == BIG_ENDIAN
#define blk0(i) block->l[i]
#else
#error "Endianness not defined!"
#endif
#define blk(i) (block->l[(i)&15U] = rol(block->l[((i)+13U)&15U] ^ \
                                        block->l[((i)+8U)&15U]  ^ \
                                        block->l[((i)+2U)&15U]   ^ \
                                        block->l[(i)&15U],1U))

/* (R0+R1), R2, R3, R4 are the different operations used in SHA1 */
#define R0(v,w,x,y,z,i) z+=(((w)&((x)^(y)))^(y))+blk0(i)+0x5A827999U+rol(v,5U);(w)=rol(w,30U);
#define R1(v,w,x,y,z,i) z+=(((w)&((x)^(y)))^(y))+blk(i)+0x5A827999U+rol(v,5U);(w)=rol(w,30U);
#define R2(v,w,x,y,z,i) z+=((w)^(x)^(y))+blk(i)+0x6ED9EBA1U+rol(v,5U);(w)=rol(w,30U);
#define R3(v,w,x,y,z,i) z+=((((w)|(x))&(y))|((w)&(x)))+blk(i)+0x8F1BBCDCU+rol(v,5U);(w)=rol(w,30U);
#define R4(v,w,x,y,z,i) z+=((w)^(x)^(y))+blk(i)+0xCA62C1D6U+rol(v,5U);(w)=rol(w,30U);


/* Hash a single 512-bit block. This is the core of the algorithm. */

void SHA1Transform(uint32_t * state, const unsigned char * buffer)
{
    uint32_t a, b, c, d, e;

    union {
	unsigned char c[64];
	uint32_t l[16];
    } block[1];	/* use array to appear as a pointer */

    // NOLINTNEXTLINE(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling)
    memcpy(block, buffer, 64);

    /* Copy context->state[] to working vars */
    a = state[0];
    b = state[1];
    c = state[2];
    d = state[3];
    e = state[4];
    /* 4 rounds of 20 operations each. Loop unrolled. */
    R0(a,b,c,d,e, 0U); R0(e,a,b,c,d, 1U); R0(d,e,a,b,c, 2U); R0(c,d,e,a,b, 3U);
    R0(b,c,d,e,a, 4U); R0(a,b,c,d,e, 5U); R0(e,a,b,c,d, 6U); R0(d,e,a,b,c, 7U);
    R0(c,d,e,a,b, 8U); R0(b,c,d,e,a, 9U); R0(a,b,c,d,e,10U); R0(e,a,b,c,d,11U);
    R0(d,e,a,b,c,12U); R0(c,d,e,a,b,13U); R0(b,c,d,e,a,14U); R0(a,b,c,d,e,15U);
    R1(e,a,b,c,d,16U); R1(d,e,a,b,c,17U); R1(c,d,e,a,b,18U); R1(b,c,d,e,a,19U);
    R2(a,b,c,d,e,20U); R2(e,a,b,c,d,21U); R2(d,e,a,b,c,22U); R2(c,d,e,a,b,23U);
    R2(b,c,d,e,a,24U); R2(a,b,c,d,e,25U); R2(e,a,b,c,d,26U); R2(d,e,a,b,c,27U);
    R2(c,d,e,a,b,28U); R2(b,c,d,e,a,29U); R2(a,b,c,d,e,30U); R2(e,a,b,c,d,31U);
    R2(d,e,a,b,c,32U); R2(c,d,e,a,b,33U); R2(b,c,d,e,a,34U); R2(a,b,c,d,e,35U);
    R2(e,a,b,c,d,36U); R2(d,e,a,b,c,37U); R2(c,d,e,a,b,38U); R2(b,c,d,e,a,39U);
    R3(a,b,c,d,e,40U); R3(e,a,b,c,d,41U); R3(d,e,a,b,c,42U); R3(c,d,e,a,b,43U);
    R3(b,c,d,e,a,44U); R3(a,b,c,d,e,45U); R3(e,a,b,c,d,46U); R3(d,e,a,b,c,47U);
    R3(c,d,e,a,b,48U); R3(b,c,d,e,a,49U); R3(a,b,c,d,e,50U); R3(e,a,b,c,d,51U);
    R3(d,e,a,b,c,52U); R3(c,d,e,a,b,53U); R3(b,c,d,e,a,54U); R3(a,b,c,d,e,55U);
    R3(e,a,b,c,d,56U); R3(d,e,a,b,c,57U); R3(c,d,e,a,b,58U); R3(b,c,d,e,a,59U);
    R4(a,b,c,d,e,60U); R4(e,a,b,c,d,61U); R4(d,e,a,b,c,62U); R4(c,d,e,a,b,63U);
    R4(b,c,d,e,a,64U); R4(a,b,c,d,e,65U); R4(e,a,b,c,d,66U); R4(d,e,a,b,c,67U);
    R4(c,d,e,a,b,68U); R4(b,c,d,e,a,69U); R4(a,b,c,d,e,70U); R4(e,a,b,c,d,71U);
    R4(d,e,a,b,c,72U); R4(c,d,e,a,b,73U); R4(b,c,d,e,a,74U); R4(a,b,c,d,e,75U);
    R4(e,a,b,c,d,76U); R4(d,e,a,b,c,77U); R4(c,d,e,a,b,78U); R4(b,c,d,e,a,79U);
    /* Add the working vars back into context.state[] */
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    /* Wipe variables */
    // NOLINTNEXTLINE(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling)
    memset(block, 0, sizeof(block));
    // This is an ill-advised attempt at clearing the registers. The C
    // model actually cannot forbid the compiler to notice that these
    // statements are not needed, and elide them. Let's just NOT do this.
    // Of course, we don't careabout security in our use case.
    // a=b=c=d=e=0;
}


/* SHA1Init - Initialize new context */

void SHA1Init(SHA1_CTX * context)
{
    /* SHA1 initialization constants */
    context->state[0] = 0x67452301;
    context->state[1] = 0xEFCDAB89;
    context->state[2] = 0x98BADCFE;
    context->state[3] = 0x10325476;
    context->state[4] = 0xC3D2E1F0;
    context->count[0] = context->count[1] = 0;
}

/* Run your data through this. */

void SHA1Update(SHA1_CTX * context, const unsigned char *data, uint32_t len)
{
    uint32_t i;
    uint32_t j;

    // NOLINTBEGIN(bugprone-*)
    j = context->count[0];
    if ((context->count[0] += len << 3U) < j)
    // NOLINTEND(bugprone-*)
	context->count[1]++;
    context->count[1] += (len >> 29U);
    j = (j >> 3U) & 63U;
    if ((j + len) > 63) {
        // NOLINTNEXTLINE(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling)
	memcpy(&context->buffer[j], data, (i = 64 - j));
	SHA1Transform(context->state, context->buffer);
	for (; i + 63 < len; i += 64) {
	    SHA1Transform(context->state, &data[i]);
	}
	j = 0;
    } else {
	i = 0;
    }
    // NOLINTNEXTLINE(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling)
    memcpy(&context->buffer[j], &data[i], len - i);
}


/* Add padding and return the message digest. */

void SHA1Final(unsigned char digest[20], SHA1_CTX * context)
{
    unsigned i;
    unsigned char finalcount[8];
    unsigned char c;

    for (i = 0; i < 8; i++) {
	finalcount[i] = (unsigned char) ((context->count[(i >= 4 ? 0 : 1)] >> ((3 - (i & 3U)) * 8)) & 255U);	/* Endian independent */
    }
    c = 0200;
    SHA1Update(context, &c, 1);
    while ((context->count[0] & 504U) != 448U) {
	c = 0000;
	SHA1Update(context, &c, 1);
    }
    SHA1Update(context, finalcount, 8);	/* Should cause a SHA1Transform() */
    for (i = 0; i < 20; i++) {
	digest[i] = (unsigned char)
	    ((context->state[i >> 2U] >> ((3 - (i & 3U)) * 8)) & 255U);
    }
    /* Wipe variables */
    // NOLINTNEXTLINE(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling)
    memset(context, '\0', sizeof(*context));
    // NOLINTNEXTLINE(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling)
    memset(&finalcount, '\0', sizeof(finalcount));
}

void SHA1(char *hash_out, const char *str, int len)
{
    SHA1_CTX ctx;
    SHA1Init(&ctx);
    SHA1Update(&ctx, (const unsigned char *) str, len);
    SHA1Final((unsigned char *) hash_out, &ctx);
    hash_out[20] = '\0';
}
