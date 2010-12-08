
set (sources 
        cado_poly.c relation.c timing.c plain_poly.c getprime.c mod_ul.c
        modredc_ul.c modredc_15ul.c modredc_2ul2.c modul_poly.c gmp_aux.c
        mpz_poly.c fpoly.c rootfinder.c params.c gcd_int64.c gcd_uint64.c
        discriminant.c random_generation.c mpz_array.c gzip.c misc.c
        hashpair.c poly.c crc.c purgedfile.c bit_vector.c
        barrier.c usp.c fix-endianness.c  
)

add_library (utils ${sources})
target_link_libraries(utils gmp)


add_executable(check_rels check_rels.c)
target_link_libraries(check_rels utils)

add_executable(complete_rels ${EXCLUDE_FROM_RELEASE_BUILD} complete_rels.c)
target_link_libraries(complete_rels utils)

# install(TARGETS complete_rels check_rels RUNTIME DESTINATION bin/utils)

add_executable(test_mod_ul ${EXCLUDE_FROM_RELEASE_BUILD} test_mod_ul.c)
add_executable(test_modredc_ul ${EXCLUDE_FROM_RELEASE_BUILD} test_modredc_ul.c)
add_executable(test_modredc_15ul ${EXCLUDE_FROM_RELEASE_BUILD} test_modredc_15ul.c)
add_executable(test_modredc_2ul2 ${EXCLUDE_FROM_RELEASE_BUILD} test_modredc_2ul2.c)
add_executable(gmp_prob_prime gmp_prob_prime.c)

target_link_libraries(test_mod_ul utils)
target_link_libraries(test_modredc_ul utils)
target_link_libraries(test_modredc_15ul utils)
target_link_libraries(test_modredc_2ul2 utils)
target_link_libraries(gmp_prob_prime gmp)

include(CMakeLists-nodist.txt OPTIONAL)
