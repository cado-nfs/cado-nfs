include(CheckCXXSourceCompiles)

set(check_shared_ptr_prefix "#include <memory>")

set(check_shared_ptr_suffix "
    struct foo {
        std::shared_ptr<int> a;
    };
    int main()
    {
        foo f;
        f.a = std::make_shared<int>();
        return 0;
    }
")

set(check_shared_ptr_copy_ctor_suffix "
class nocopy {
        protected:
                nocopy() {}
                ~nocopy() {}
        private:
                nocopy(const nocopy&);
                nocopy& operator=(const nocopy&);
};
class foo: private nocopy {
        int val;
        public:
        foo(int val) : val(val) {}
        int value() const { return val; }
};
std::shared_ptr<foo> blah()
{
        const int val = 65536;
        return std::make_shared<foo>(val);
}
int main()
{
    return blah()->value();
}
")

CHECK_CXX_SOURCE_COMPILES("${check_shared_ptr_prefix}
${check_shared_ptr_suffix}" HAVE_STD_SHARED_PTR_BASIC)

if(HAVE_STD_SHARED_PTR_BASIC)
    CHECK_CXX_SOURCE_COMPILES("${check_shared_ptr_prefix}
    ${check_shared_ptr_copy_ctor_suffix}" HAS_NOT_BUG_21397_STD_VERSION)
    if(HAS_NOT_BUG_21397_STD_VERSION)
        set(HAVE_STD_SHARED_PTR 1)
    else()
        message(STATUS "See bug #21397 -- std::shared_ptr can't be used, trying boost:: now")
    endif()
endif()

set(check_shared_ptr_boost_infix "
    #include <boost/shared_ptr.hpp>
    #include <boost/make_shared.hpp>
    namespace std {
        using boost::shared_ptr;
        using boost::make_shared;
    }
")

if(NOT HAS_NOT_BUG_21397_STD_VERSION)
    # In a sense, we're ready to believe that the boost version is almost
    # certainly fine, but we haven't checked whether it's always so. The
    # std:: version, however, is definitely not fine under certain
    # circumstances. See bug #21397
    CHECK_CXX_SOURCE_COMPILES("${check_shared_ptr_prefix}
    ${check_shared_ptr_boost_infix}
    ${check_shared_ptr_suffix}" HAVE_BOOST_SHARED_PTR_BASIC)

    if(HAVE_BOOST_SHARED_PTR_BASIC)
        CHECK_CXX_SOURCE_COMPILES("${check_shared_ptr_prefix}
        ${check_shared_ptr_boost_infix}
        ${check_shared_ptr_copy_ctor_suffix}"
        HAS_NOT_BUG_21397_BOOST_VERSION)
        if(HAS_NOT_BUG_21397_BOOST_VERSION)
            set(HAVE_BOOST_SHARED_PTR 1)
        else()
            message(STATUS "See bug #21397 -- boost::shared_ptr can't be used either")
        endif()
    endif()
endif()

if (NOT HAVE_STD_SHARED_PTR AND NOT HAVE_BOOST_SHARED_PTR)
    message(FATAL_ERROR "Need either a C++11 compiler, or boost libraries installed. You may simply unpack the boost headers at the root of the cado tree, that will probably do. Also, most OS distributions include packages for boost, which can be easily installed (names may vary, e.g.  libboost-dev or boost-devel")
endif()
