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
    void doit(cado::work_queue & Q) const {
        fmt::print("bob({})\n", x);
        /* This one is funny since we might actually be calling
         * ourselves. */
        if (Q.do_one_task())
            fmt::print("stolen one task! :-)\n");
    }
};


int main()
{
    bob B { 42 };
    cado::work_queue Q(4);

    for(int i = 0 ; i < 100 ; i++) {
        Q.push_task(say_hi);
        Q.push_task(say_hi2, 12 + i);
        Q.push_task([&B, &Q](){B.doit(Q);});
    }
}


