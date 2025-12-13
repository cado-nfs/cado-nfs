#include <map>
#include <stdint.h>

struct foo {
	uint64_t z[18];
};

using stupidmap = std::multimap<int, foo>;

void do_bug(stupidmap & s, int a, foo const & r)
{
    s.emplace(a, r);
}

int main()
{
    return 0;
}

/*

g++ -O2 -mavx512f -mavx512vl  -c bug30027.cpp

/tmp/cc1AweNa.s: Assembler messages:
/tmp/cc1AweNa.s:37: Error: unsupported instruction `vmovdqu'

gcc version 8.3.0 (Debian 8.3.0-6)

 */
