#include <stdlib.h>
int main()
{
    void *res = NULL;
    for (size_t i = 1; i <= 2340; i++)
        {
            int rc = posix_memalign (&res, 64, i);
            if (rc != 0 || (((unsigned long) res) % 64) != 0)
                __builtin_abort ();
        }
    return 0;
}
