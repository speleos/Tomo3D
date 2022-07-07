/*
 * util3d.h based on util.h
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#ifndef _JK_UTIL3D_H_
#define _JK_UTIL3D_H_

#include "geom3d.h"
#include "array.h"

const int MaxStr=512;
const char SepChars[]=" ";

void readBC(char* ifn, Array1d<RegularBC3d*>& velBC);
int countLines(const char *fn);
double intp1d(Array1d<double>& x, Array1d<double>& f, double xval);

#endif /* _JK_UTIL3D_H_ */
