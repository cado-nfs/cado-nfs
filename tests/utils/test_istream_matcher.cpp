#include "cado.h"

#include <sstream>
#include <iostream>
#include <stdexcept>
#include "istream_matcher.hpp"
#include "macros.h"


int main()
{
    {
        std::string data = "Catch 22";

        {
            std::istringstream is(data);
            int x;
            istream_matcher<char> m(is);
            bool b = (m >> "Catch" >> x && x == 22);
            ASSERT_ALWAYS(b);
        }

        {
            std::istringstream is(data);
            int x;
            istream_matcher<char> m(is);
            bool b = (m >> "Catch" >> std::noskipws >> x && x == 22);
            ASSERT_ALWAYS(!b);
        }

        {
            std::istringstream is(data);
            std::string prefix("Catch");
            int x;
            istream_matcher<char> m(is);
            bool b = (m >> (std::string const &) prefix >> std::noskipws >> x && x == 22);
            ASSERT_ALWAYS(!b);
        }
    }

    {
        std::string data = "Catch22";

        {
            std::istringstream is(data);
            int x;
            istream_matcher<char> m(is);
            bool b = (m >> "Catch" >> x && x == 22);
            ASSERT_ALWAYS(b);
        }

        {
            std::istringstream is(data);
            int x;
            istream_matcher<char> m(is);
            bool b = (m >> "Catch" >> std::noskipws >> x && x == 22);
            ASSERT_ALWAYS(b);
        }
    }
}
