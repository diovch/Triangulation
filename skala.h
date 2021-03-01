#ifndef SKALA_H
#define SKALA_H

#include "r2geom.h"
#include "R3Graph.h"
#include "Triangulation.h"

void skalaMethod(
    double (*f)(const R3Graph::R3Point&),
    const R3Graph::R3Box& box,
    int numX, int numY, int numZ,
    Triangulation& triangulation
);

R3Graph::R3Vector gradientR3(
    double (*f)(const R3Graph::R3Point&),
    const R3Graph::R3Point& p,
    double h
);

#endif
