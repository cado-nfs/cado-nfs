#include "cado.h" // IWYU pragma: keep

#include <cstdlib>

#include <iostream>
#include <sstream>
#include <vector>
#include <utility>
#include <string>

#include "json.hpp"

static const std::vector<std::pair<std::string, std::string>> examples {
    { R"({ "Image": { "Width":  800, "Height": 600, "Title":  "View from 15th Floor", "Thumbnail": { "Url":    "http://www.example.com/image/481989943", "Height": 125, "Width":  100 }, "Animated" : false, "IDs": [116, 943, 234, 38793] }})", "" },
    { R"([ { "precision": "zip", "Latitude":  37.7668, "Longitude": -122.3959, "Address":   "", "City":      "SAN FRANCISCO", "State":     "CA", "Zip":       "94107", "Country":   "US" }, { "precision": "zip", "Latitude":  37.371991, "Longitude": -122.026020, "Address":   "", "City":      "SUNNYVALE", "State":     "CA", "Zip":       "94085", "Country":   "US" } ])", "" },
};

static const std::vector<std::string> expected_failures {
    "[,]", "[],{}", "{a:2}",
};
int main()
{
    for(auto const & example : examples) {
        json f;
        if (!(std::istringstream(example.first) >> f)) {
            std::cerr << "cannot parse json: " << example.first << "\n";
            exit(EXIT_FAILURE);
        }
        std::cout << f << "\n";
    }
    for(auto const & example : expected_failures) {
        json f;
        if ((std::istringstream(example) >> f)) {
            std::cerr << "unexpected success while parsing bad json\n";
            exit(EXIT_FAILURE);
        }
    }
}
