#ifndef SHA1_H_
#define SHA1_H_

/*
   SHA-1 in C
   By Steve Reid <steve@edmweb.com>
   100% Public Domain

   small edits for cado-nfs by E. Thomé. Still PD.
 */

#include <stdint.h>

typedef struct
{
    uint32_t state[5];
    uint32_t count[2];
    unsigned char buffer[64];
} SHA1_CTX;

#ifdef __cplusplus
extern "C" {
#endif

void SHA1Transform(uint32_t * state, const unsigned char * buffer);
void SHA1Init(SHA1_CTX * context); 
void SHA1Update(SHA1_CTX * context, const unsigned char *data, uint32_t len);
void SHA1Final(unsigned char digest[20], SHA1_CTX * context);
void SHA1(char *hash_out, const char *str, int len);

#ifdef __cplusplus
}
#endif

/* end "SHA-1 in C" PD interface */

#ifdef __cplusplus
#include <streambuf>
#include <ostream>
#include <string>
/* mock streambuf implementation that just computes the sha-1 sum of what
 * it gets fed.
 *
 * idea of hacking around streambufs picked from
 * http://wordaligned.org/articles/cpp-streambufs (PD)
 *
 */
template <typename char_type,
          typename traits = std::char_traits<char_type> >
class sha1_checksumming_streambuf: public std::basic_streambuf<char_type, traits>
{
    SHA1_CTX ctx;
public:
    typedef typename traits::int_type int_type;
    sha1_checksumming_streambuf() { SHA1Init(&ctx); }
    void checksum(char out[41]) const
    {
        SHA1_CTX ctx1;
        memcpy(&ctx1, &ctx, sizeof(SHA1_CTX));
        uint8_t dest[20];
        SHA1Final(dest, &ctx1);
        for(int i = 0 ; i < 20 ; i++)
            snprintf(out + 2*i, 3, "%02x", (unsigned int) dest[i]);
    }
    std::string digest() const
    {
        char x[41] = { '\0' };
        checksum(x);
        return std::string(x);
    }
private:
    virtual int_type overflow(int_type c) {
        if (traits::eq_int_type(c, traits::eof()))
            return traits::not_eof(c);
        char_type const ch = traits::to_char_type(c);
        SHA1Update(&ctx, (unsigned char*) &ch, sizeof(char_type));
        return c;
    }
};

class sha1_checksumming_stream : public std::ostream
{
    sha1_checksumming_streambuf<char> cbuf;
public:
    sha1_checksumming_stream() : std::ostream(&cbuf) {}
    void checksum(char out[41]) const { cbuf.checksum(out); }
    std::string digest() const { return cbuf.digest(); }
};
#endif


#endif	/* SHA1_H_ */

