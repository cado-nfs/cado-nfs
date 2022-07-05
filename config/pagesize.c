#include <stdio.h>
#include <unistd.h>
#include <stddef.h>

int main()
{
    size_t s = sysconf(_SC_PAGESIZE);
    printf("%zu\n", s);
    return 0;
}


