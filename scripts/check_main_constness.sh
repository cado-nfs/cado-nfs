#!/bin/bash

if git grep '\bmain.*int.*argc.*argv' | grep -v '\(config\|out-of-core\|gf2x\)/' | grep -v 'char const' ; then
    echo "ERRROR: the functions above should have the prototype"
    echo "      int main (int argc, char const * argv[])"
    echo "Rationale: true, it's not mandated by C, but more often"
    echo "than not, doing so fosters good const-coorectness in the"
    echo "functions that are called"
    exit 1
fi
