#ifndef FIX_ENDIANNESS_H_
#define FIX_ENDIANNESS_H_

#ifdef __cplusplus
extern "C" {
#endif

size_t fread32_little(uint32_t * ptr, size_t nmemb, FILE * stream);
size_t fwrite32_little(const uint32_t * ptr, size_t nmemb, FILE * stream);
size_t fread64_little(uint64_t * ptr, size_t nmemb, FILE * stream);
size_t fwrite64_little(const uint64_t * ptr, size_t nmemb, FILE * stream);

#ifdef __cplusplus
}
#endif

#endif	/* FIX_ENDIANNESS_H_ */
