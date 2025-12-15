#include "cado.h"       // IWYU pragma: keeep

#include "work_queue.hpp"

#include "fmt/base.h"

static void say_hi()
{
    fmt::print("Hello world\n");
}

static void say_hi2(int a)
{
    fmt::print("Hello world ({})\n", a);
}

struct bob {
    int x;
    void doit() const {
        fmt::print("bob({})\n", x);
    }
};


int main()
{
    bob B { 42 };
    cado::work_queue Q(4);

    for(int i = 0 ; i < 100 ; i++) {
        Q.push_task(say_hi);
        Q.push_task(say_hi2, 12 + i);
        Q.push_task([&B](){B.doit();});
    }
}


