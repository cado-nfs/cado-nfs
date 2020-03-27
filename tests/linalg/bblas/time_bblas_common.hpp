#ifndef TIME_BBLAS_COMMON_HPP_
#define TIME_BBLAS_COMMON_HPP_

#include <ctime>
#include <cstdio>
// #include <type_traits>
#include <utility>
#include <map>

extern int test_accel;  /* declaration is in test_bblas */

#if 0
#define t_and_unit_from_clock__bare(t, unit, t1, j)                     \
    double t = t1;							\
    t /= CLOCKS_PER_SEC;						\
    if (j) t /= j; else t = 0;						\
    const char * unit = "s";						\
    if (t < 1.0e-7) { unit = "ns"; t *= 1.0e9;			        \
    } else if (t < 1.0e-4) { unit = "micros"; t *= 1.0e6;		\
    } else if (t < 1.0e-1) { unit = "ms"; t *= 1.0e3; }                 \
    do { } while (0)

#define TIME1__bare(maxtime, what, args) 		        	\
    clock_t measuring_time = maxtime * CLOCKS_PER_SEC / test_accel;	\
    clock_t t0, t1;							\
    int j;								\
    t0 = clock();							\
    for (j = 0; ; j++) {						\
        what args;							\
        t1 = clock() - t0;						\
        if (j && t1 > measuring_time)					\
            break;							\
    }									\
    t_and_unit_from_clock__bare(t, unit, t1, j);

#define TIME1(maxtime, what, args) do {			        	\
    TIME1__bare(maxtime, what, args)                                    \
    printf(#what " \t%d times in %.4f %s each\n",       		\
            j, t, unit);		                        	\
} while (0)

#define TIME1N(maxtime, what, args) do {		        	\
    TIME1__bare(maxtime, what, args)                                    \
    printf(#what "(n=%d) \t%d times in %.4f %s each\n", n,     	        \
            j, t, unit);		                        	\
} while (0)

#define TIME1N_spins(rexpr, maxtime, what, args, spinexpr, spinmax) do {        \
    clock_t ts[spinmax];						\
    int ns[spinmax];                                                    \
    for(int s = 0 ; s < spinmax ; s++) ts[s] = ns[s] = 0;		\
    clock_t t0, t1;							\
    int j;								\
    t0 = clock();							\
    clock_t t = t0;                                                     \
    clock_t fence = t0 + maxtime * CLOCKS_PER_SEC / test_accel;		\
    for (j = 0; ; j++) {						\
        rexpr;                                                          \
        int ret = what args;						\
        int s = spinexpr;                                               \
        ts[s] += (t1 = clock()) - t;                                    \
        ns[s] ++;                                                       \
        t = t1;                                                         \
        if (j && t1 > fence)					        \
            break;							\
    }									\
    int nch=0;                                                          \
    for(int s = 0 ; s < spinmax ; s++) {				\
        if (s == 0) nch=printf(#what"(n=%d)", n);			\
        else for(int k = nch ; k-- ; putchar(' '));		        \
        t_and_unit_from_clock__bare(t, unit, ts[s], ns[s]);		\
        printf(" \t[%d] %d times in %.4f %s each\n",s,ns[s],t,unit);	\
    }									\
} while (0)
#else

struct bblas_timer {
    clock_t measuring_time;
    clock_t t0, t1;
    int j;
    double maxtime;
    const char * name;
    double t;
    const char * unit;
    bblas_timer(double maxtime, const char * name)
        : maxtime(maxtime)
          , name(name)
    {
        measuring_time = maxtime * CLOCKS_PER_SEC;
        if (test_accel) measuring_time /= 100;
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
            t0 = clock();
            for (j = 0; ; j++) {
                f(std::forward<Args>(args)...);
                t1 = clock() - t0;
                if (j && t1 > measuring_time)
                    break;
            }
            get_unit();
        }
    template<typename T, typename... Args>
        void time1(T const & f, Args... args)
        {
            time1_common(f, std::forward<Args>(args)...);
            printf("%s \t%d times in %.4f %s each\n", name, j, t, unit);
        }
    template<typename T, typename... Args>
        void time1n(int n, T const & f, Args... args)
        {
            time1_common(f, std::forward<Args>(args)...);
            printf("%s(n=%d) \t%d times in %.4f %s each\n", name, n, j, t, unit);
        }

    template<typename R, typename T, typename... Args>
        void time1n_classify(int n, R const & rr, T const & f, Args... args) {
            // typedef decltype(std::declval<T>()(f(std::declval<Args>()...))) U;
            std::map<int, std::pair<int, clock_t>> ts;
            t0 = clock();
            clock_t tx = t0;
            clock_t fence = t0 + maxtime * CLOCKS_PER_SEC / test_accel;
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
                if (nch == 0) nch=printf("%s(n=%d)", name, n);
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

#endif

#endif	/* TIME_BBLAS_COMMON_HPP_ */
