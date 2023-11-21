#include <sys/endian.h>
#include <stdint.h>
#include <stdlib.h>

int main()
{
    uint32_t a = 0x12345678;
    a = bswap32(a);
    return a == 0x78563412 ? EXIT_SUCCESS : EXIT_FAILURE;
}


