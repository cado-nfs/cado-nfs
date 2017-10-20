set (SQ_SRC special-q_sieve.cpp mat_int64.c int64_vector.c uint64_array.c
    mat_Z.c ideal.c makefb.c utils_mpz.c int64_poly.c sieving_bound.c
    array.c utils_int64.c factor_base.c utils_lattice.c double_vector.c
    mat_double.c list_int64_vector.c list_double_vector.c
    list_int64_vector_index.c utils_norm.c utils_cofactorisation.cpp)

add_executable(special-q_sieve ${SQ_SRC})

set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS SLLL_SAFE)

include_directories (${CMAKE_SOURCE_DIR}/sieve)
target_link_libraries(special-q_sieve facul utils)

# Makefb

set(MAKEFB_SRC makefb.c factor_base.c ideal.c utils_int64.c)

add_executable(makefb-hd ${MAKEFB_SRC})

set_property(TARGET makefb-hd PROPERTY COMPILE_DEFINITIONS MAIN)

target_link_libraries(makefb-hd utils)

# Polynomial selection

#set (POL_SRC polyselect-hd.c polynomials.c sieving_bound.c utils_mpz.c)

#add_executable(polyselect-hd ${POL_SRC})

#target_link_libraries(polyselect-hd utils polyselect_common)

#set_property(TARGET polyselect-hd PROPERTY COMPILE_DEFINITIONS EEA_BOUND)

# Parameters selection

#set (PAR_SRC paraselect-hd.c parameters.c polynomials.c sieving_bound.c utils_mpz.c)

#add_executable(paraselect-hd ${PAR_SRC})

#target_link_libraries(paraselect-hd utils polyselect_common)

# Alpha 3d
set (ALP_SRC alpha3d.c)
add_executable(alpha3d ${ALP_SRC})
set_property(TARGET alpha3d PROPERTY COMPILE_DEFINITIONS MAIN_ALPHA3D)
target_link_libraries(alpha3d utils)

# Murphy-E
set (MUR_SRC murphyE3d.c alpha3d.c)
add_executable(murphyE3d ${MUR_SRC})
set_property(TARGET murphyE3d PROPERTY COMPILE_DEFINITIONS MAIN_MURPHYE3D)
target_link_libraries(murphyE3d utils polyselect_common)

# Test LLL int64
set (LLL_INT64 utils_lattice.c utils_int64.c mat_int64.c mat_Z.c int64_vector.c int64_poly.c list_int64_vector.c list_int64_vector_index.c double_vector.c list_double_vector.c array.c)
add_executable(lll_int64 EXCLUDE_FROM_ALL ${LLL_INT64})
set_property(TARGET lll_int64 PROPERTY COMPILE_DEFINITIONS MAIN_LLL_INT64)
target_link_libraries(lll_int64 utils)
