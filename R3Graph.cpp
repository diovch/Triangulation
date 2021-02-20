#include "R3Graph.h"

namespace R3Graph {

bool intersectPlanes(
    const R3Point& p0, const R3Vector& n0,
    const R3Point& p1, const R3Vector& n1,
    R3Point& p, R3Vector& v
) {
    if (n0.area(n1) <= R3_EPSILON)
        return false;
    v = n0.vectorProduct(n1).normalize();
    R3Vector v1 = n0.vectorProduct(v);

    // p = p0 + v1*t
    // ((p0 + v1*t) - p1, n1) = 0
    // (p0 - p1, n1) + (v1, n1)*t = 0
    // t = (p1 - p0, n1) / (v1, n1)
    double s = v1.scalarProduct(n1);
    if (fabs(s) <= R3_EPSILON)
        return false;
    double t = (p1 - p0).scalarProduct(n1) / s;
    p = p0 + v1 * t;
    return true;
}

bool intersectPlaneAndLine(
    const R3Point& p0, const R3Vector& n,
    const R3Point& p1, const R3Vector& v,
    R3Point& p
) {
    double s = n.scalarProduct(v);
    if (fabs(s) <= R3_EPSILON)
        return false;
    // p = p1 + v*t;
    // (p - p0, n) = 0
    // (p1 + v*t - p0, n) = 0
    // (p0 - p1, n) = (v, n)*t
    // t = (p0 - p1, n) / (v, n)
    double t = (p0 - p1).scalarProduct(n) / s;
    p = p1 + v*t;
    return true;
}

} // end of namespace R3Graph
