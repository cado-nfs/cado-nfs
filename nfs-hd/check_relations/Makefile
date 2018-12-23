CXX=g++
HEADERS=/usr/local/include
build_tree=/tmp
PATH_LIBS=/usr/local/lib
LIBS=-lgmp -lpthread
FLAGS=-Wall -Wextra -Werror -pedantic -Wpedantic -Wconversion -O3 -DNDEBUG

all:
	$(CXX) -I$(HEADERS) -I$(HEADERS_CADO_H) -I$(HEADERS_CADO_CONFIG_H) check_relations.cpp -o $(build_tree)/check_relations -L$(PATH_LIBS) -lntl $(LIBS) $(FLAGS)
