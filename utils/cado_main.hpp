#ifndef CADO_MAIN_HPP
#define CADO_MAIN_HPP

/* Programs that report a fatal error by throwing -- with cado::error, or
 * with cxx_param_list::fail() -- need something to catch it, or the
 * process dies on std::terminate with no useful message and no chance to
 * flush anything. Wrap the body of main() like this:
 *
 *   static int main_(int argc, char const * argv[]);
 *
 *   int main(int argc, char const * argv[])
 *   {
 *       return cado::main_wrapper(main_, argc, argv);
 *   }
 *
 *   static int main_(int argc, char const * argv[])
 *   {
 *       ... unchanged body ...
 *   }
 */

namespace cado {
int main_wrapper(int (*f)(int, char const **), int argc, char const * argv[]);
}

#endif /* CADO_MAIN_HPP */
