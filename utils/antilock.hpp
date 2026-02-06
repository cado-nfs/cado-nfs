#ifndef UTILS_ANTILOCK_HPP_
#define UTILS_ANTILOCK_HPP_

#include <mutex>
#include <system_error>

namespace cado {
    template<class Mutex>
        struct antilock {
            /* This can be used to temporarily release a unique lock,
             * e.g. like what condition_variable::wait does.
             */
            std::unique_lock<Mutex> & u;
            antilock(antilock const &) = delete;
            antilock& operator=(antilock const &) = delete;
            antilock(antilock &&) = default;
            antilock& operator=(antilock &&) = default;
            explicit antilock(std::unique_lock<Mutex> & u)
                : u(u)
            {
                if (!u.owns_lock())
                    throw std::system_error(std::make_error_code(std::errc::operation_not_permitted));
                u.unlock();
            }
            ~antilock() { u.lock(); }
        };
} /* namespace cado */

#endif	/* UTILS_ANTILOCK_HPP_ */
