#include "cado.h" // IWYU pragma: keep
#include <iostream>
#include <string>
#include <unistd.h>
#include "gzip.h"

// coverity[root_function]
int main(int argc, char * argv[])
{
    const char * filename = "test.gz";

    if (argc > 2 && std::string(argv[1]) == "--wdir") {
        chdir(argv[2]);
        argc--,argv++;
        argc--,argv++;
    }

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
