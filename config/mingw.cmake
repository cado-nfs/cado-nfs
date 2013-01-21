# check if we are under MinGW

message(STATUS "Checking for MinGW")
if(MINGW)
   message(STATUS "Checking for MinGW -- yes")
   set(HAVE_MINGW 1)
else(MINGW)
   message(STATUS "Checking for MinGW -- no")
   set(HAVE_MINGW 0)
endif(MINGW)
