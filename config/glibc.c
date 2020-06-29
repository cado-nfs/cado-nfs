#include <stdio.h>
#ifndef __GLIBC__
#error "not glibc"
#endif
int main() { printf("Hi!\n"); return 0; }
