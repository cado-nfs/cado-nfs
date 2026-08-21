#include <valarray>
#include <ranges>
import llvm_bug_190333;
#include <vector>

int main() {
    std::vector<int> v {{1, 2, 3}};
    auto p = [](auto const & e) { return e + 1; };
    auto c = v | std::views::transform(p);
    return 0;
}
