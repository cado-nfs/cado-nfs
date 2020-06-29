#include <stdlib.h>
int main()
{
    return aligned_alloc(64, 1024) != NULL;
}
