#include <compare>
#include <functional>
#include <iostream>

struct S
{
    int val;

    S(int v) : val(v) {}

    operator const void *() const { std::cout << "cast\n"; return &val; }

    friend bool operator==(S lhs, S rhs) noexcept 
    { std::cout << "op==\n"; return lhs.val == rhs.val; }
    friend std::strong_ordering operator<=>(S lhs, S rhs) noexcept 
    { std::cout << "op<=>\n"; return lhs.val <=> rhs.val; }  
};

int main()
{
    const S arr[] = {S{2}, S{1}};
    // In C++20 mode it compares pointers, and so considers that arr[1] > arr[0],
    // which is wrong!
    return std::greater_equal<>{}(arr[0], arr[1]) ? 0 : 1;
}
