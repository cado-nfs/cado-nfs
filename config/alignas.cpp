#include <cstdint>
int main()
{
    alignas(8) uint32_t x = 0;
    return x;
}
