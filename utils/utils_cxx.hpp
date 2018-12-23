#ifndef UTILS_CXX_HPP_
#define UTILS_CXX_HPP_

/* Base class with private copy-constructor and assignment operator.
   Classes which are not copy-constructible can inherit this with:
   private NonCopyable.
   
   See also:

   https://www.boost.org/doc/libs/1_67_0/libs/core/doc/html/core/noncopyable.html
   */
struct NonCopyable {
   NonCopyable() {}
   ~NonCopyable() {}
   NonCopyable(NonCopyable&&) = delete;
   NonCopyable& operator=(NonCopyable&&) = delete;
   NonCopyable(NonCopyable const &) = delete;
   NonCopyable& operator=(NonCopyable const &) = delete;
};

/* This is an example, but not a class one can inherit from: inheriting
 * does not have the same effect as having this kind of ctors declared...
 *
class MoveOnly
{
public:
   MoveOnly() {}
   ~MoveOnly() {}
   MoveOnly(const MoveOnly&) = delete;
   MoveOnly& operator=(const MoveOnly&) = delete;
   MoveOnly(MoveOnly&&) = default;
   MoveOnly& operator=(MoveOnly&&) = default;
};
 */

/* Generic tool to cause an arbitrary closure to be called at destructor
 * time. See detached_cofac for a use case.
 */
template<typename T> struct call_dtor_s {
    T x;
    call_dtor_s(T x): x(x) {}
    ~call_dtor_s() { x(); }
};
template<typename T> call_dtor_s<T> call_dtor(T x) { return call_dtor_s<T>(x); }


#endif	/* UTILS_CXX_HPP_ */
