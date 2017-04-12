#include "cado.h"
#include <iostream>
#include <iterator>
#include "compressible_heap.hpp"

int main()
{
    int blah[20];
    compressible_heap<int, int, 32> H;

    for(int i = 0 ; i < 2000 ; i++) {
        int n = rand() % 20;
        for(int j = 0 ; j < n ; j++)
            blah[j] = rand() % 97;
        H.push_back(blah, n);
    }

    /* This is not really leaving 2000-1900 = 100 rows in the heap,
     * because erasing by absolute indices makes it possible to erase a
     * relation which is already no longer present.
     *
     * Instead, the expected number of remaining rows is
     * 2000*exp(-1900/2000), which is about 773.
     */
    for(int i = 0 ; i < 1900 ; i++)
        H.kill(rand() % H.size());

    /* add a few more */
    for(int i = 0 ; i < 2000 ; i++) {
        int n = rand() % 20;
        for(int j = 0 ; j < n ; j++)
            blah[j] = rand() % 97;
        H.push_back(blah, n);
    }

    /* and kill some again.
     *
     * same idea about expected population at the end, except that we've
     * done some erasures already. Half of our deletes will stack on top
     * of those we've done already on the first batch.
     *
     * So we expect to have 2000*(exp(-9900/2000)+exp(-8000/2000))
     *
     * which is about 50 relations.
     */
    for(int i = 0 ; i < 16000 ; i++)
        H.kill(rand() % H.size());

    size_t n = 0;
    for(auto it = H.begin() ; it != H.end() ; ++it, n++) {
        std::cout << "row " << it.index() << ": "
                << "length " << it->second << ": ";
        std::copy(it->first, it->first + it->second, std::ostream_iterator<int>(std::cout, " "));
        std::cout << "\n";
    }
    std::cout << "Found "<<n<<" rows in total\n";
}
