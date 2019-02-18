# check for mallopt in "malloc.h"

include(${PROJECT_SOURCE_DIR}/config/search_for_function.cmake)
search_for_function(mallopt HAVE_MALLOPT)

