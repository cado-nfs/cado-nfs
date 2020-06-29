    #include <stdlib.h>
    int main()
    {
        srand(1);
        unsigned long a = rand();
        srand(1);
        unsigned long b = rand();
        return a == b ? EXIT_SUCCESS : EXIT_FAILURE;
    }
