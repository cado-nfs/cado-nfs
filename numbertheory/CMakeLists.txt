
add_library (ant ant.cpp)
target_link_libraries (ant utils)

add_executable(badideals badideals.cpp)
target_link_libraries (badideals ant)

install(TARGETS badideals RUNTIME DESTINATION ${LIBSUFFIX}/numbertheory)
if(ENABLE_SHARED)
install(TARGETS ant DESTINATION ${LIBSUFFIX}/numbertheory)
endif()
