#ifndef TIME_BBLAS_COMMON_HPP_
#define TIME_BBLAS_COMMON_HPP_

#include <ctime>
#include <cstdio>
// #include <type_traits>
#include <utility>
#include <map>
#include <string>
#include <functional>

/* We used to do this with C macros. It's not harder to do the same with
 * C++ templates.
 */
struct bblas_timer {
    clock_t measuring_time;
    clock_t t0=0, t1=0;
    int j=0;
    double maxtime;
    std::string const & name;
    double t;
    const char * unit="s";
    bblas_timer(double maxtime, std::string const & name)
        : maxtime(maxtime)
          , name(name)
    {
        measuring_time = maxtime * CLOCKS_PER_SEC;
        if (test_bblas_base::test_accel) measuring_time /= 100;
    }
    void get_unit()
    {
        /* run this after testit() */
        t = t1;
        t /= CLOCKS_PER_SEC;
        if (j) t /= j; else t = 0;
        unit = "s";
        if (t < 1.0e-7) { unit = "ns"; t *= 1.0e9;			
        } else if (t < 1.0e-4) { unit = "micros"; t *= 1.0e6;
        } else if (t < 1.0e-1) { unit = "ms"; t *= 1.0e3; }
    }
    template<typename T, typename... Args>
        void time1_common(T const & f, Args&&... args)
        {
            auto F = std::bind(f, args...);

            t0 = clock();
            for (j = 0; ; j++) {
                F();
                t1 = clock() - t0;
                if (j && t1 > measuring_time)
                    break;
            }
            get_unit();
        }
    template<typename T, typename... Args>
        void time1(T const & f, Args&&... args)
        {
            time1_common(f, std::forward<Args>(args)...);
            printf("%s \t%d times in %.4f %s each\n", name.c_str(), j, t, unit);
        }
    template<typename T, typename... Args>
        void time1n(int n, T const & f, Args&&... args)
        {
            time1_common(f, std::forward<Args>(args)...);
            printf("%s(n=%d) \t%d times in %.4f %s each\n", name.c_str(), n, j, t, unit);
        }

    template<typename R, typename T, typename... Args>
        void time1n_classify(int n, R const & rr, T const & f, Args&&... args) {
            // typedef decltype(std::declval<T>()(f(std::declval<Args>()...))) U;
            std::map<int, std::pair<int, clock_t>> ts;
            t0 = clock();
            clock_t tx = t0;
            clock_t fence = maxtime * CLOCKS_PER_SEC;
            if (test_bblas_base::test_accel) fence /= 100;
            fence += t0;
            for (j = 0; ; j++) {
                rr();
                int ret = f(std::forward<Args>(args)...);
                ts[ret].first++;
                ts[ret].second+= (t1 = clock()) - tx;
                tx = t1;
                if (j && t1 > fence)					
                    break;
            }
            int nch=0;
            for(auto const & s : ts) {
                int key = s.first;
                if (nch == 0) nch=printf("%s(n=%d)", name.c_str(), n);
                else for(int k = nch ; k-- ; putchar(' '));		
                j  = s.second.first;
                t1 = s.second.second;
                get_unit();
                get_unit();
                printf(" \t[%d] %d times in %.4f %s each\n",key,j,t,unit);
            }
        }
};

#define TIME1(maxtime, what, ...) bblas_timer(maxtime, #what).time1(what, __VA_ARGS__)
#define TIME1N(maxtime, what, ...) bblas_timer(maxtime, #what).time1n(n, what, __VA_ARGS__)
#define TIME1N_SPINS(randomizer_code, maxtime, what, ...) bblas_timer(maxtime, #what).time1n_classify(n, [&](){randomizer_code;}, what, __VA_ARGS__)

#endif	/* TIME_BBLAS_COMMON_HPP_ */
