/*
 * error.cc
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <iostream>
#include <cstdlib>
#include "error.hpp"

void error(string s)
{
    cerr << s << '\n';
    std::abort();
}

