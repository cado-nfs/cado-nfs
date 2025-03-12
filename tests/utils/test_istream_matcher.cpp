#include "cado.h" // IWYU pragma: keep

#include <sstream>
#include <iostream>
#include <string>

#include "istream_matcher.hpp"
#include "macros.h"


int main()
{
    {
        const std::string data = "Catch 22";

        {
            std::istringstream is(data);
            int x;
            istream_matcher m(is);
            ASSERT_ALWAYS(m >> "Catch" >> x && x == 22);
        }

        {
            std::istringstream is(data);
            int x;
            istream_matcher m(is);
            ASSERT_ALWAYS(!(m >> "Catch" >> std::noskipws >> x));
        }

        {
            std::istringstream is(data);
            const std::string prefix("Catch");
            int x;
            istream_matcher m(is);
            ASSERT_ALWAYS(!(m >> prefix >> std::noskipws >> x));
        }
    }

    {
        const std::string data = "Catch22";

        {
            std::istringstream is(data);
            int x;
            istream_matcher m(is);
            ASSERT_ALWAYS(m >> "Catch" >> x && x == 22);
        }

        {
            std::istringstream is(data);
            int x;
            istream_matcher m(is);
            ASSERT_ALWAYS(m >> "Catch" >> std::noskipws >> x && x == 22);
        }
    }
}
