#include "cado.h" // IWYU pragma: keep
#include <cstdlib>
#include <iostream>
#include <string>
#include <unistd.h>
#include "gzip.h"

// coverity[root_function]
int main(int argc, char * argv[])
{
    const char * filename = "test.gz";
    const char * t = getenv("wdir");

    if (t)
        chdir(t);

    if (argc > 1)
        filename = argv[1];

    ofstream_maybe_compressed os(filename);
    os << "Hello, world\n";
    os.close();

    ifstream_maybe_compressed is(filename);
    std::string s;
    getline(is, s);
    std::cout << s << "\n";
}
