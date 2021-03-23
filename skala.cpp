#include <vector>
#include <map>
#include <cstring>
#include <cassert>

#include "skala.h"

using namespace R3Graph;

static double dh = 0.01;
static const double EPS2 = R3_EPSILON * 10.;

class ScalaMatrix {
private:
    int nx;
    int ny;
    double* v;

private:
    ScalaMatrix() :
        nx(0),
        ny(0),
        v(0)
    {}

public:
    ScalaMatrix(int sx, int sy) :
        nx(sx),
        ny(sy),
        v(new double[sx * sy])
    {}

    ScalaMatrix(const ScalaMatrix& m) :
        nx(m.nx),
        ny(m.ny),
        v(new double[nx * ny])
    {
        memmove(v, m.v, nx * ny * sizeof(double));
    }

    ~ScalaMatrix() { delete[] v; }

    ScalaMatrix& operator=(const ScalaMatrix& m) {
        if (nx * ny > m.nx * m.ny) {
            delete[] v;
            v = new double[m.nx * m.ny];
        }
        nx = m.nx;
        ny = m.ny;
        memmove(v, m.v, nx * ny * sizeof(double));
        return *this;
    }

    int sizeX() const { return nx; }
    int sizeY() const { return ny; }

    double* operator[](int iy) {
        return v + iy * nx;
    }
    const double* operator[](int iy) const {
        return v + iy * nx;
    }

    double& at(int x, int y) {
        return v[y * nx + x];
    }
    const double& at(int x, int y) const {
        return v[y * nx + x];
    }
};

class TPoint {
public:
    int ix;
    int iy;
    int iz;

    TPoint() :
        ix(0),
        iy(0),
        iz(0)
    {}

    TPoint(int x, int y, int z) :
        ix(x),
        iy(y),
        iz(z)
    {}

    bool operator==(const TPoint& tp) const {
        return (ix == tp.ix && iy == tp.iy && iz == tp.iz);
    }

    bool operator!=(const TPoint& tp) const {
        return (ix != tp.ix || iy != tp.iy || iz != tp.iz);
    }

    bool operator<(const TPoint& tp) const {
        return (
            ix < tp.ix || (
                ix == tp.ix && (
                    iy < tp.iy || (
                        iy == tp.iy &&
                        iz < tp.iz
                        )
                    )
                )
            );
    }

    bool operator<=(const TPoint& tp) const {
        return (
            ix < tp.ix || (
                ix == tp.ix && (
                    iy < tp.iy || (
                        iy == tp.iy &&
                        iz <= tp.iz
                        )
                    )
                )
            );
    }

    bool operator>(const TPoint& tp) const {
        return !operator<=(tp);
    }

    bool operator>=(const TPoint& tp) const {
        return !operator<(tp);
    }

    R3Point coord(
        const R3Point& origin,
        double dx, double dy, double dz
    ) const {
        return R3Point(
            origin.x + dx * (double)ix,
            origin.y + dy * (double)iy,
            origin.z + dz * (double)iz
        );
    }
};

class TEdge {
public:
    TPoint point[2];
    int vertexIdx;

    TEdge() :
        vertexIdx(-1)
    {}

    TEdge(const TPoint& p0, const TPoint& p1) :
        vertexIdx(-1)
    {
        assert(p0 != p1);
        if (p0 <= p1) {
            point[0] = p0;
            point[1] = p1;
        }
        else {
            point[0] = p1;
            point[1] = p0;
        }
        assert(point[0] < point[1]);
    }

    bool operator==(const TEdge& te) const {
        return (
            point[0] == te.point[0] &&
            point[1] == te.point[1]
            );
    }

    bool operator!=(const TEdge& te) const {
        return (
            point[0] != te.point[0] ||
            point[1] != te.point[1]
            );
    }

    bool operator<(const TEdge& te) const {
        return (
            point[0] < te.point[0] || (
                point[0] == te.point[0] &&
                point[1] < te.point[1]
                )
            );
    }

    bool operator<=(const TEdge& te) const {
        return (
            point[0] < te.point[0] || (
                point[0] == te.point[0] &&
                point[1] <= te.point[1]
                )
            );
    }

    bool operator>(const TEdge& te) const {
        return !operator<=(te);
    }

    bool operator>=(const TEdge& te) const {
        return !operator<(te);
    }

    bool computeVertex(
        double v0, double v1,
        const R3Point& origin,
        double dx, double dy, double dz,
        R3Point& vertex
    ) {
        if (
            (v0 > 0. && v1 > 0.) ||
            (v0 < 0. && v1 < 0.)
            ) {
            return false;
        }
        double w0 = fabs(v0);
        if (w0 <= R3_EPSILON) {
            vertex = point[0].coord(
                origin, dx, dy, dz
            );
            return true;
        }
        double w1 = fabs(v1);
        if (w1 <= R3_EPSILON) {
            vertex = point[1].coord(
                origin, dx, dy, dz
            );
            return true;
        }
        R3Point p0 = point[0].coord(
            origin, dx, dy, dz
        );
        R3Point p1 = point[1].coord(
            origin, dx, dy, dz
        );
        // vertex = p0 + (p1 - p0)*(w0/(w0 + w1));
        vertex = p0 + (p1 - p0) * (v0 / (v0 - v1));
        return true;
    }
};

class TBox : public R3Box {
public:
    int nx;
    int ny;
    int nz;

    double dx;
    double dy;
    double dz;

private:
    TBox();

public:
    TBox(
        const R3Point& o, const R3Vector& s,
        int numX, int numY, int numZ
    ) :
        R3Box(o, s),
        nx(numX),
        ny(numY),
        nz(numZ)
    {
        assert(nx > 0);
        dx = size.x / (double)nx;

        assert(ny > 0);
        dy = size.y / (double)ny;

        assert(nz > 0);
        dz = size.z / (double)nz;
    }

    void setNx(int n) {
        assert(n > 0);
        nx = n;
        dx = size.x / (double)n;
    }

    void setNy(int n) {
        assert(n > 0);
        ny = n;
        dy = size.y / (double)n;
    }

    void setNz(int n) {
        assert(n > 0);
        nz = n;
        dz = size.z / (double)n;
    }

    R3Point gridPoint(int ix, int iy, int iz) const {
        return origin + R3Vector(dx * ix, dy * iy, dz * iz);
    }
};

void processEdgeX(
    TEdge& edge,
    const ScalaMatrix& redMatrix0,
    const ScalaMatrix& blueMatrix0,
    const ScalaMatrix& blueMatrix1,
    const TBox& box,
    std::map<TEdge, int>& edges,
    Triangulation& triangulation
);

void processEdgeY(
    TEdge& edge,
    const ScalaMatrix& redMatrix0,
    const ScalaMatrix& blueMatrix0,
    const ScalaMatrix& blueMatrix1,
    const TBox& box,
    std::map<TEdge, int>& edges,
    Triangulation& triangulation
);

void processEdgeZ(
    TEdge& edge,
    const ScalaMatrix& redMatrix0,
    const ScalaMatrix& redMatrix1,
    const ScalaMatrix& blueMatrix1,
    const TBox& box,
    std::map<TEdge, int>& edges,
    Triangulation& triangulation
);

void processTetrahedron(
    const TEdge& redEdge,
    double redValue0, double redValue1,
    const TEdge& blueEdge,
    double blueValue0, double blueValue1,
    const TBox& box,
    std::map<TEdge, int>& edges,
    Triangulation& triangulation
);

void skalaMethod(
    double (*f)(const R3Point&),
    const R3Box& box,
    int numX, int numY, int numZ,
    Triangulation& triangulation
) {
    triangulation.clear();

    double dx = box.size.x / (numX + 1.);
    double dy = box.size.y / (numY + 1.);
    double dz = box.size.z / (numZ + 1.);

    R3Vector vdx(dx, 0., 0.);
    R3Vector vdy(0., dy, 0.);
    R3Vector vdz(0., 0., dz);

    double dx2 = dx / 2.;
    double dy2 = dy / 2.;
    double dz2 = dz / 2.;

    R3Vector vdx2(dx2, 0., 0.);
    R3Vector vdy2(0., dy2, 0.);
    R3Vector vdz2(0., 0., dz2);

    int numX2 = numX * 2;
    int numY2 = numY * 2;
    int numZ2 = numZ * 2;

    R3Point blueOrigin = box.origin;
    R3Point redOrigin = blueOrigin + R3Vector(dx2, dy2, dz2);

    // For gradient calculation
    if (dh > dx / 10.)
        dh = dx / 10.;

    TBox tbox(
        redOrigin, box.size,
        numX2 + 2, numY2 + 2, numZ2 + 2
    );

    // Values at the bottom
    ScalaMatrix redMatrix0(numX + 1, numY + 1);
    ScalaMatrix redMatrix1(numX + 1, numY + 1);
    for (int iy = 0; iy <= numY; ++iy) {
        R3Point row0 = redOrigin + vdy * iy;
        R3Point row1 = row0 + vdz;
        for (int ix = 0; ix <= numX; ++ix) {
            //... redMatrix0.at(ix, iy) = (*f)(row0);
            double vr0 = (*f)(row0);
            if (fabs(vr0) <= R3_EPSILON)
                vr0 = EPS2;
            redMatrix0.at(ix, iy) = vr0;

            //... redMatrix1.at(ix, iy) = (*f)(row1);
            double vr1 = (*f)(row1);
            if (fabs(vr1) <= R3_EPSILON)
                vr1 = EPS2;
            redMatrix1.at(ix, iy) = vr1;

            row0 += vdx;
            row1 += vdx;
        }
    }

    ScalaMatrix blueMatrix0(numX + 2, numY + 2);
    ScalaMatrix blueMatrix1(numX + 2, numY + 2);
    for (int iy = 0; iy <= numY + 1; ++iy) {
        R3Point row0 = blueOrigin + vdy * iy;
        R3Point row1 = row0 + vdz;
        for (int ix = 0; ix <= numX + 1; ++ix) {
            //... blueMatrix0.at(ix, iy) = (*f)(row0);
            double vr0 = (*f)(row0);
            if (fabs(vr0) <= R3_EPSILON)
                vr0 = EPS2;
            blueMatrix0.at(ix, iy) = vr0;

            //... blueMatrix1.at(ix, iy) = (*f)(row1);
            double vr1 = (*f)(row1);
            if (fabs(vr1) <= R3_EPSILON)
                vr1 = EPS2;
            blueMatrix1.at(ix, iy) = vr1;

            row0 += vdx;
            row1 += vdx;
        }
    }

    std::map<TEdge, int> edges;

    R3Point redOriginZ = redOrigin + vdz;
    R3Point blueOriginZ = blueOrigin + vdz;

    for (int iz = 1; iz < numZ2; iz += 2) {
        for (int iy = 1; iy < numY2; iy += 2) {
            // Horizontal edges X: *---*---*--- ... *---*
            for (int ix = 1; ix < numX2 - 2; ix += 2) {
                TEdge edgeX(
                    TPoint(ix, iy, iz),
                    TPoint(ix + 2, iy, iz)
                );
                processEdgeX(
                    edgeX,
                    redMatrix0,
                    blueMatrix0,
                    blueMatrix1,
                    tbox,
                    edges,
                    triangulation
                );
            } // end for (ix...

            // Horizontal edges Y: |  |  | ... |
            if (iy != numY2 - 1) {
                for (int ix = 1; ix < numX2; ix += 2) {
                    TEdge edgeY(
                        TPoint(ix, iy, iz),
                        TPoint(ix, iy + 2, iz)
                    );
                    processEdgeY(
                        edgeY,
                        redMatrix0,
                        blueMatrix0,
                        blueMatrix1,
                        tbox,
                        edges,
                        triangulation
                    );
                } // end for (ix...
            } // end if

            // Vertical edges Z
            if (iz != numZ2 - 1) {
                for (int ix = 1; ix < numX2; ix += 2) {
                    TEdge edgeZ(
                        TPoint(ix, iy, iz),
                        TPoint(ix, iy, iz + 2)
                    );
                    processEdgeZ(
                        edgeZ,
                        redMatrix0,
                        redMatrix1,
                        blueMatrix1,
                        tbox,
                        edges,
                        triangulation
                    );
                } // end for (ix...
            } // end if
        } // end for (iy...

        redMatrix0 = redMatrix1;
        blueMatrix0 = blueMatrix1;

        redOriginZ += vdz;
        for (int iy = 0; iy <= numY; ++iy) {
            R3Point row1 = redOriginZ + vdy * iy;
            for (int ix = 0; ix <= numX; ++ix) {
                //... redMatrix1.at(ix, iy) = (*f)(row1);
                double vr1 = (*f)(row1);
                if (fabs(vr1) <= R3_EPSILON)
                    vr1 = EPS2;
                redMatrix1.at(ix, iy) = vr1;

                row1 += vdx;
            }
        }

        blueOriginZ += vdz;
        for (int iy = 0; iy <= numY; ++iy) {
            R3Point row1 = blueOriginZ + vdy * iy;
            for (int ix = 0; ix <= numX; ++ix) {
                //... blueMatrix1.at(ix, iy) = (*f)(row1);
                double vr1 = (*f)(row1);
                if (fabs(vr1) <= R3_EPSILON)
                    vr1 = EPS2;
                blueMatrix1.at(ix, iy) = vr1;

                row1 += vdx;
            }
        }

    } // end for (iz...

    // Define gradients
    size_t num_vertices = triangulation.vertices.size();
    for (size_t i = 0; i < num_vertices; ++i) {
        Triangulation::Vertex& v =
            triangulation.vertices.at(i);
        v.normal = (-gradientR3(
            f, v.point, dh
        ));
        v.normal.normalize();
    }

    triangulation.orientate();
}

void processEdgeX(
    TEdge& edgeX,
    const ScalaMatrix& redMatrix0,
    const ScalaMatrix& blueMatrix0,
    const ScalaMatrix& blueMatrix1,
    const TBox& tbox,
    std::map<TEdge, int>& edges,
    Triangulation& triangulation
) {
    assert(
        (edgeX.point[0].ix & 1) != 0 && // Odd number
        (edgeX.point[0].iy & 1) != 0 &&
        (edgeX.point[1].iy & 1) != 0 &&
        (edgeX.point[1].iy & 1) != 0
    );

    assert(
        edgeX.point[0].iy == edgeX.point[1].iy &&
        edgeX.point[0].iz == edgeX.point[1].iz
    );

    int ixRed0 = (edgeX.point[0].ix - 1) / 2;
    int iyRed0 = (edgeX.point[0].iy - 1) / 2;

    int ixRed1 = (edgeX.point[1].ix - 1) / 2;
    int iyRed1 = (edgeX.point[1].iy - 1) / 2;

    double v0 = redMatrix0.at(ixRed0, iyRed0);
    double v1 = redMatrix0.at(ixRed1, iyRed1);
    bool vertexOnEdgeExists = false;
    R3Point vertexOnEdge;

    if (
        (v0 <= 0. && v1 >= 0.) ||
        (v0 >= 0. && v1 <= 0.)
        ) {
        if (edges.count(edgeX) == 0) {
            vertexOnEdgeExists = edgeX.computeVertex(
                v0, v1,
                tbox.origin,
                tbox.dx, tbox.dy, tbox.dz,
                vertexOnEdge
            );
            assert(vertexOnEdgeExists);
            assert(tbox.contains(vertexOnEdge));

            // Add a vertex to triangulation
            triangulation.vertices.push_back(vertexOnEdge);
            int vertexIdx = (int)(triangulation.vertices.size() - 1);
            edgeX.vertexIdx = vertexIdx;        // for safety
            edges[edgeX] = vertexIdx;
        }
        else {
            if (edgeX.vertexIdx < 0) {
                edgeX.vertexIdx = edges[edgeX];
            }
            assert(edgeX.vertexIdx == edges[edgeX]);
        }
    }

    // Surrounding edges
    TEdge edgesYZ[4];
    double vertexValue[4];
    double edgeValue[4][2];
    int blue_ix = edgeX.point[0].ix + 1;
    int blue_iy = edgeX.point[0].iy - 1;
    int blue_iz = edgeX.point[0].iz - 1;
    assert(
        (blue_ix & 1) == 0 &&   // even number
        (blue_iy & 1) == 0 &&
        (blue_iz & 1) == 0
    );

    vertexValue[0] = blueMatrix0.at(
        blue_ix / 2, blue_iy / 2
    );
    vertexValue[1] = blueMatrix0.at(
        blue_ix / 2, blue_iy / 2 + 1
    );
    vertexValue[2] = blueMatrix1.at(
        blue_ix / 2, blue_iy / 2
    );
    vertexValue[3] = blueMatrix1.at(
        blue_ix / 2, blue_iy / 2 + 1
    );

    edgesYZ[0] = TEdge(
        TPoint(blue_ix, blue_iy, blue_iz),
        TPoint(blue_ix, blue_iy + 2, blue_iz)
    );
    edgeValue[0][0] = vertexValue[0];
    edgeValue[0][1] = vertexValue[1];

    edgesYZ[1] = TEdge(
        TPoint(blue_ix, blue_iy, blue_iz + 2),
        TPoint(blue_ix, blue_iy + 2, blue_iz + 2)
    );
    edgeValue[1][0] = vertexValue[2];
    edgeValue[1][1] = vertexValue[3];

    edgesYZ[2] = TEdge(
        TPoint(blue_ix, blue_iy, blue_iz),
        TPoint(blue_ix, blue_iy, blue_iz + 2)
    );
    edgeValue[2][0] = vertexValue[0];
    edgeValue[2][1] = vertexValue[2];

    edgesYZ[3] = TEdge(
        TPoint(blue_ix, blue_iy + 2, blue_iz),
        TPoint(blue_ix, blue_iy + 2, blue_iz + 2)
    );
    edgeValue[3][0] = vertexValue[1];
    edgeValue[3][1] = vertexValue[3];

    R3Point vertexYZ;
    for (int i = 0; i < 4; ++i) {
        double w0 = edgeValue[i][0];
        double w1 = edgeValue[i][1];

        if (
            (w0 <= 0. && w1 >= 0.) ||
            (w0 >= 0. && w1 <= 0.)
            ) {
            if (edges.count(edgesYZ[i]) == 0) {
                bool vertexOK = edgesYZ[i].computeVertex(
                    w0, w1,
                    tbox.origin,
                    tbox.dx, tbox.dy, tbox.dz,
                    vertexYZ
                );
                assert(vertexOK);
                assert(tbox.contains(vertexYZ));

                // Add a vertex to triangulation
                triangulation.vertices.push_back(vertexYZ);
                int vertexIdx = (int)(triangulation.vertices.size() - 1);
                edgesYZ[i].vertexIdx = vertexIdx;       // for safety
                edges[edgesYZ[i]] = vertexIdx;
            }
            else {
                if (edgesYZ[i].vertexIdx < 0) {
                    edgesYZ[i].vertexIdx = edges[edgesYZ[i]];
                }
                assert(edgesYZ[i].vertexIdx == edges[edgesYZ[i]]);
            }
        }

        processTetrahedron(
            edgeX, v0, v1,
            edgesYZ[i], w0, w1,
            tbox,
            edges,
            triangulation
        );
    } // end for
}

void processEdgeY(
    TEdge& edgeY,
    const ScalaMatrix& redMatrix0,
    const ScalaMatrix& blueMatrix0,
    const ScalaMatrix& blueMatrix1,
    const TBox& tbox,
    std::map<TEdge, int>& edges,
    Triangulation& triangulation
) {
    assert(
        (edgeY.point[0].ix & 1) != 0 && // Odd number
        (edgeY.point[0].iy & 1) != 0 &&
        (edgeY.point[1].iy & 1) != 0 &&
        (edgeY.point[1].iy & 1) != 0
    );

    assert(
        edgeY.point[0].ix == edgeY.point[1].ix &&
        edgeY.point[0].iz == edgeY.point[1].iz
    );

    int ixRed0 = (edgeY.point[0].ix - 1) / 2;
    int iyRed0 = (edgeY.point[0].iy - 1) / 2;

    int ixRed1 = (edgeY.point[1].ix - 1) / 2;
    int iyRed1 = (edgeY.point[1].iy - 1) / 2;

    double v0 = redMatrix0.at(ixRed0, iyRed0);
    double v1 = redMatrix0.at(ixRed1, iyRed1);
    bool vertexOnEdgeExists = false;
    R3Point vertexOnEdge;

    if (
        (v0 <= 0. && v1 >= 0.) ||
        (v0 >= 0. && v1 <= 0.)
        ) {
        if (edges.count(edgeY) == 0) {
            vertexOnEdgeExists = edgeY.computeVertex(
                v0, v1,
                tbox.origin,
                tbox.dx, tbox.dy, tbox.dz,
                vertexOnEdge
            );
            assert(vertexOnEdgeExists);
            assert(tbox.contains(vertexOnEdge));

            // Add a vertex to triangulation
            triangulation.vertices.push_back(vertexOnEdge);
            int vertexIdx = (int)(triangulation.vertices.size() - 1);
            edgeY.vertexIdx = vertexIdx;        // for safety
            edges[edgeY] = vertexIdx;
        }
        else {
            if (edgeY.vertexIdx < 0) {
                edgeY.vertexIdx = edges[edgeY];
            }
            assert(edgeY.vertexIdx == edges[edgeY]);
        }
    }

    // Surrounding edges
    TEdge edgesXZ[4];
    double vertexValue[4];
    double edgeValue[4][2];
    int blue_ix = edgeY.point[0].ix - 1;
    int blue_iy = edgeY.point[0].iy + 1;
    int blue_iz = edgeY.point[0].iz - 1;
    assert(
        (blue_ix & 1) == 0 &&   // even number
        (blue_iy & 1) == 0 &&
        (blue_iz & 1) == 0
    );

    vertexValue[0] = blueMatrix0.at(
        blue_ix / 2, blue_iy / 2
    );
    vertexValue[1] = blueMatrix0.at(
        blue_ix / 2 + 1, blue_iy / 2
    );
    vertexValue[2] = blueMatrix1.at(
        blue_ix / 2, blue_iy / 2
    );
    vertexValue[3] = blueMatrix1.at(
        blue_ix / 2 + 1, blue_iy / 2
    );

    edgesXZ[0] = TEdge(
        TPoint(blue_ix, blue_iy, blue_iz),
        TPoint(blue_ix + 2, blue_iy, blue_iz)
    );
    edgeValue[0][0] = vertexValue[0];
    edgeValue[0][1] = vertexValue[1];

    edgesXZ[1] = TEdge(
        TPoint(blue_ix, blue_iy, blue_iz + 2),
        TPoint(blue_ix + 2, blue_iy, blue_iz + 2)
    );
    edgeValue[1][0] = vertexValue[2];
    edgeValue[1][1] = vertexValue[3];

    edgesXZ[2] = TEdge(
        TPoint(blue_ix, blue_iy, blue_iz),
        TPoint(blue_ix, blue_iy, blue_iz + 2)
    );
    edgeValue[2][0] = vertexValue[0];
    edgeValue[2][1] = vertexValue[2];

    edgesXZ[3] = TEdge(
        TPoint(blue_ix + 2, blue_iy, blue_iz),
        TPoint(blue_ix + 2, blue_iy, blue_iz + 2)
    );
    edgeValue[3][0] = vertexValue[1];
    edgeValue[3][1] = vertexValue[3];

    R3Point vertexXZ;
    for (int i = 0; i < 4; ++i) {
        double w0 = edgeValue[i][0];
        double w1 = edgeValue[i][1];

        if (
            (w0 <= 0. && w1 >= 0.) ||
            (w0 >= 0. && w1 <= 0.)
            ) {
            if (edges.count(edgesXZ[i]) == 0) {
                bool vertexOK = edgesXZ[i].computeVertex(
                    w0, w1,
                    tbox.origin,
                    tbox.dx, tbox.dy, tbox.dz,
                    vertexXZ
                );
                assert(vertexOK);
                assert(tbox.contains(vertexXZ));

                // Add a vertex to triangulation
                triangulation.vertices.push_back(vertexXZ);
                int vertexIdx = (int)(triangulation.vertices.size() - 1);
                edgesXZ[i].vertexIdx = vertexIdx;       // for safety
                edges[edgesXZ[i]] = vertexIdx;
            }
            else {
                if (edgesXZ[i].vertexIdx < 0) {
                    edgesXZ[i].vertexIdx = edges[edgesXZ[i]];
                }
                assert(edgesXZ[i].vertexIdx == edges[edgesXZ[i]]);
            }
        }

        processTetrahedron(
            edgeY, v0, v1,
            edgesXZ[i], w0, w1,
            tbox,
            edges,
            triangulation
        );
    } // end for
}

void processEdgeZ(
    TEdge& edgeZ,
    const ScalaMatrix& redMatrix0,
    const ScalaMatrix& redMatrix1,
    const ScalaMatrix& blueMatrix1,
    const TBox& tbox,
    std::map<TEdge, int>& edges,
    Triangulation& triangulation
) {
    assert(
        (edgeZ.point[0].ix & 1) != 0 && // Odd number
        (edgeZ.point[0].iy & 1) != 0 &&
        (edgeZ.point[1].iy & 1) != 0 &&
        (edgeZ.point[1].iy & 1) != 0
    );

    assert(
        edgeZ.point[0].ix == edgeZ.point[1].ix &&
        edgeZ.point[0].iy == edgeZ.point[1].iy
    );

    int ixRed0 = (edgeZ.point[0].ix - 1) / 2;
    int iyRed0 = (edgeZ.point[0].iy - 1) / 2;

    int ixRed1 = (edgeZ.point[1].ix - 1) / 2;
    int iyRed1 = (edgeZ.point[1].iy - 1) / 2;

    assert(ixRed0 == ixRed1 && iyRed0 == iyRed1);

    double v0 = redMatrix0.at(ixRed0, iyRed0);
    double v1 = redMatrix1.at(ixRed1, iyRed1);
    bool vertexOnEdgeExists = false;
    R3Point vertexOnEdge;

    if (
        (v0 <= 0. && v1 >= 0.) ||
        (v0 >= 0. && v1 <= 0.)
        ) {
        if (edges.count(edgeZ) == 0) {
            vertexOnEdgeExists = edgeZ.computeVertex(
                v0, v1,
                tbox.origin,
                tbox.dx, tbox.dy, tbox.dz,
                vertexOnEdge
            );
            assert(vertexOnEdgeExists);
            assert(tbox.contains(vertexOnEdge));

            // Add a vertex to triangulation
            triangulation.vertices.push_back(vertexOnEdge);
            int vertexIdx = (int)(triangulation.vertices.size() - 1);
            edgeZ.vertexIdx = vertexIdx;        // for safety
            edges[edgeZ] = vertexIdx;
        }
        else {
            if (edgeZ.vertexIdx < 0) {
                edgeZ.vertexIdx = edges[edgeZ];
            }
            assert(edgeZ.vertexIdx == edges[edgeZ]);
        }
    }

    // Surrounding edges
    TEdge edgesXY[4];
    double vertexValue[4];
    double edgeValue[4][2];
    int blue_ix = edgeZ.point[0].ix - 1;
    int blue_iy = edgeZ.point[0].iy - 1;
    int blue_iz = edgeZ.point[0].iz + 1;
    assert(
        (blue_ix & 1) == 0 &&   // even numbers
        (blue_iy & 1) == 0 &&
        (blue_iz & 1) == 0
    );

    vertexValue[0] = blueMatrix1.at(
        blue_ix / 2, blue_iy / 2
    );
    vertexValue[1] = blueMatrix1.at(
        blue_ix / 2 + 1, blue_iy / 2
    );
    vertexValue[2] = blueMatrix1.at(
        blue_ix / 2, blue_iy / 2 + 1
    );
    vertexValue[3] = blueMatrix1.at(
        blue_ix / 2 + 1, blue_iy / 2 + 1
    );

    edgesXY[0] = TEdge(
        TPoint(blue_ix, blue_iy, blue_iz),
        TPoint(blue_ix + 2, blue_iy, blue_iz)
    );
    edgeValue[0][0] = vertexValue[0];
    edgeValue[0][1] = vertexValue[1];

    edgesXY[1] = TEdge(
        TPoint(blue_ix, blue_iy + 2, blue_iz),
        TPoint(blue_ix + 2, blue_iy + 2, blue_iz)
    );
    edgeValue[1][0] = vertexValue[2];
    edgeValue[1][1] = vertexValue[3];

    edgesXY[2] = TEdge(
        TPoint(blue_ix, blue_iy, blue_iz),
        TPoint(blue_ix, blue_iy + 2, blue_iz)
    );
    edgeValue[2][0] = vertexValue[0];
    edgeValue[2][1] = vertexValue[2];

    edgesXY[3] = TEdge(
        TPoint(blue_ix + 2, blue_iy, blue_iz),
        TPoint(blue_ix + 2, blue_iy + 2, blue_iz)
    );
    edgeValue[3][0] = vertexValue[1];
    edgeValue[3][1] = vertexValue[3];

    R3Point vertexXY;
    for (int i = 0; i < 4; ++i) {
        double w0 = edgeValue[i][0];
        double w1 = edgeValue[i][1];

        if (
            (w0 <= 0. && w1 >= 0.) ||
            (w0 >= 0. && w1 <= 0.)
            ) {
            if (edges.count(edgesXY[i]) == 0) {
                bool vertexOK = edgesXY[i].computeVertex(
                    w0, w1,
                    tbox.origin,
                    tbox.dx, tbox.dy, tbox.dz,
                    vertexXY
                );
                assert(vertexOK);
                assert(tbox.contains(vertexXY));

                // Add a vertex to triangulation
                triangulation.vertices.push_back(vertexXY);
                int vertexIdx = (int)(triangulation.vertices.size() - 1);
                edgesXY[i].vertexIdx = vertexIdx;       // for safety
                edges[edgesXY[i]] = vertexIdx;
            }
            else {
                if (edgesXY[i].vertexIdx < 0) {
                    edgesXY[i].vertexIdx = edges[edgesXY[i]];
                }
                assert(edgesXY[i].vertexIdx == edges[edgesXY[i]]);
            }
        }

        processTetrahedron(
            edgeZ, v0, v1,
            edgesXY[i], w0, w1,
            tbox,
            edges,
            triangulation
        );
    } // end for
}

void processTetrahedron(
    const TEdge& redEdge,
    double redValue0, double redValue1,
    const TEdge& blueEdge,
    double blueValue0, double blueValue1,
    const TBox& tbox,
    std::map<TEdge, int>& edges,
    Triangulation& triangulation
) {
    assert(
        (redEdge.point[0].ix & 1) != 0 &&
        (redEdge.point[0].iy & 1) != 0 &&
        (redEdge.point[1].ix & 1) != 0 &&
        (redEdge.point[1].iy & 1) != 0
    );
    assert(
        (blueEdge.point[0].ix & 1) == 0 &&
        (blueEdge.point[0].iy & 1) == 0 &&
        (blueEdge.point[1].ix & 1) == 0 &&
        (blueEdge.point[1].iy & 1) == 0
    );

    assert(redEdge.point[0] != redEdge.point[1]);
    assert(blueEdge.point[0] != blueEdge.point[1]);

    assert(
        redEdge.point[0].ix == redEdge.point[1].ix ||
        redEdge.point[0].iy == redEdge.point[1].iy ||
        redEdge.point[0].iz == redEdge.point[1].iz
    );
    assert(
        blueEdge.point[0].ix == blueEdge.point[1].ix ||
        blueEdge.point[0].iy == blueEdge.point[1].iy ||
        blueEdge.point[0].iz == blueEdge.point[1].iz
    );

    // Consider all 16 cases
    if (
        (
            redValue0 < 0. &&
            redValue1 < 0. &&
            blueValue0 < 0. &&
            blueValue1 < 0.
            )
        ||
        (
            redValue0 > 0. &&
            redValue1 > 0. &&
            blueValue0 > 0. &&
            blueValue1 > 0.
            )
        )
        return; // Nothing to do

    TPoint p[4];
    double w[4];
    int numPoints = 0;
    int vertexIdx[4];

    if (
        (
            redValue0 > 0. &&   // <--
            redValue1 < 0. &&
            blueValue0 < 0. &&
            blueValue1 < 0.
            )
        ||
        (
            redValue0 < 0. &&   // <--
            redValue1 > 0. &&
            blueValue0 > 0. &&
            blueValue1 > 0.
            )
        ) {
        numPoints = 1;
        p[0] = redEdge.point[0]; w[0] = redValue0;
        p[1] = redEdge.point[1]; w[1] = redValue1;
        p[2] = blueEdge.point[0]; w[2] = blueValue0;
        p[3] = blueEdge.point[1]; w[3] = blueValue1;
    }
    else if (
        (
            redValue0 < 0. &&
            redValue1 > 0. &&   // <--
            blueValue0 < 0. &&
            blueValue1 < 0.
            )
        ||
        (
            redValue0 > 0. &&
            redValue1 < 0. &&   // <--
            blueValue0 > 0. &&
            blueValue1 > 0.
            )
        ) {
        numPoints = 1;
        p[0] = redEdge.point[1]; w[0] = redValue1;
        p[1] = redEdge.point[0]; w[1] = redValue0;
        p[2] = blueEdge.point[0]; w[2] = blueValue0;
        p[3] = blueEdge.point[1]; w[3] = blueValue1;
    }
    else if (
        (
            redValue0 < 0. &&
            redValue1 < 0. &&
            blueValue0 > 0. &&  // <--
            blueValue1 < 0.
            )
        ||
        (
            redValue0 > 0. &&
            redValue1 > 0. &&
            blueValue0 < 0. &&  // <--
            blueValue1 > 0.
            )
        ) {
        numPoints = 1;
        p[0] = blueEdge.point[0]; w[0] = blueValue0;
        p[1] = blueEdge.point[1]; w[1] = blueValue1;
        p[2] = redEdge.point[0]; w[2] = redValue0;
        p[3] = redEdge.point[1]; w[3] = redValue1;
    }
    else if (
        (
            redValue0 < 0. &&
            redValue1 < 0. &&
            blueValue0 < 0. &&
            blueValue1 > 0.     // <--
            )
        ||
        (
            redValue0 > 0. &&
            redValue1 > 0. &&
            blueValue0 > 0. &&
            blueValue1 < 0.     // <--
            )
        ) {
        numPoints = 1;
        p[0] = blueEdge.point[1]; w[0] = blueValue1;
        p[1] = blueEdge.point[0]; w[1] = blueValue0;
        p[2] = redEdge.point[0]; w[2] = redValue0;
        p[3] = redEdge.point[1]; w[3] = redValue1;
    }
    else if (
        (
            redValue0 > 0. &&   // <--
            redValue1 > 0. &&   // <--
            blueValue0 < 0. &&
            blueValue1 < 0.
            )
        ||
        (
            redValue0 < 0. &&   // <--
            redValue1 < 0. &&   // <--
            blueValue0 > 0. &&
            blueValue1 > 0.
            )
        ) {
        numPoints = 2;
        p[0] = redEdge.point[0]; w[0] = redValue0;
        p[1] = redEdge.point[1]; w[1] = redValue1;
        p[2] = blueEdge.point[0]; w[2] = blueValue0;
        p[3] = blueEdge.point[1]; w[3] = blueValue1;
    }
    else if (
        (
            redValue0 > 0. &&   // <--
            redValue1 < 0. &&
            blueValue0 > 0. &&  // <--
            blueValue1 < 0.
            )
        ||
        (
            redValue0 < 0. &&   // <--
            redValue1 > 0. &&
            blueValue0 < 0. &&  // <--
            blueValue1 > 0.
            )
        ) {
        numPoints = 2;
        p[0] = redEdge.point[0]; w[0] = redValue0;
        p[1] = blueEdge.point[0]; w[1] = blueValue0;
        p[2] = redEdge.point[1]; w[2] = redValue1;
        p[3] = blueEdge.point[1]; w[3] = blueValue1;
    }
    else if (
        (
            redValue0 > 0. &&   // <--
            redValue1 < 0. &&
            blueValue0 < 0. &&
            blueValue1 > 0.     // <--
            )
        ||
        (
            redValue0 < 0. &&   // <--
            redValue1 > 0. &&
            blueValue0 > 0. &&
            blueValue1 < 0.     // <--
            )
        ) {
        numPoints = 2;
        p[0] = redEdge.point[0]; w[0] = redValue0;
        p[1] = blueEdge.point[1]; w[1] = blueValue1;
        p[2] = redEdge.point[1]; w[2] = redValue1;
        p[3] = blueEdge.point[0]; w[3] = blueValue0;
    }

    if (numPoints == 1) {
        TEdge star[3];
        double v[3][2];
        int starLen = 0;
        for (int i = 1; i < 4; ++i) {
            if (p[0] < p[i]) {
                star[starLen] = TEdge(p[0], p[i]);
                v[starLen][0] = w[0];
                v[starLen][1] = w[i];
            }
            else {
                star[starLen] = TEdge(p[i], p[0]);
                v[starLen][0] = w[i];
                v[starLen][1] = w[0];
            }
            ++starLen;
        }
        assert(starLen == 3);

        for (int i = 0; i < 3; ++i) {
            if (edges.find(star[i]) == edges.end()) {
                // Must be a diagonal edge
                assert(
                    star[i].point[0].ix != star[i].point[1].ix &&
                    star[i].point[0].iy != star[i].point[1].iy &&
                    star[i].point[0].iz != star[i].point[1].iz
                );

                R3Point vertexOnEdge;
                bool edgeHasVertex = star[i].computeVertex(
                    v[i][0], v[i][1],
                    tbox.origin,
                    tbox.dx, tbox.dy, tbox.dz,
                    vertexOnEdge
                );
                assert(edgeHasVertex);
                assert(tbox.contains(vertexOnEdge));

                // Add a vertex to triangulation
                triangulation.vertices.push_back(vertexOnEdge);
                vertexIdx[i] = (int)(triangulation.vertices.size() - 1);
                edges[star[i]] = vertexIdx[i];
            }
            else {
                vertexIdx[i] = edges[star[i]];
            }
        }

        /*
        triangulation.triangles.push_back(
            Triangulation::Triangle(
                vertexIdx[0], vertexIdx[1], vertexIdx[2]
            )
        );
        */
        R3Point t0 = p[0].coord(
            tbox.origin,
            tbox.dx, tbox.dy, tbox.dz
        );
        R3Point t1 = triangulation.vertices[vertexIdx[0]].point;
        R3Point t2 = triangulation.vertices[vertexIdx[1]].point;
        R3Point t3 = triangulation.vertices[vertexIdx[2]].point;
        R3Vector n = (t2 - t1).vectorProduct(t3 - t1);
        double dotprod = n * (t0 - t1);
        //... if ((dotprod >= 0. && w[0] <= 0.) || (dotprod <= 0. && w[0] >= 0.)) {
        if (
            (dotprod >= 0. && w[0] >= 0.) ||
            (dotprod <= 0. && w[0] <= 0.)
            ) {
            triangulation.triangles.push_back(
                Triangulation::Triangle(
                    vertexIdx[0], vertexIdx[1], vertexIdx[2]
                )
            );
        }
        else {
            triangulation.triangles.push_back(
                Triangulation::Triangle(
                    vertexIdx[0], vertexIdx[2], vertexIdx[1]
                )
            );
        }
    }
    else if (numPoints == 2) {
        TEdge star[4];
        double v[4][2];
        int starLen = 0;
        for (int i = 2; i <= 3; ++i) {
            if (p[0] < p[i]) {
                star[starLen] = TEdge(p[0], p[i]);
                v[starLen][0] = w[0];
                v[starLen][1] = w[i];
            }
            else {
                star[starLen] = TEdge(p[i], p[0]);
                v[starLen][0] = w[i];
                v[starLen][1] = w[0];
            }
            ++starLen;
        }
        for (int i = 2; i <= 3; ++i) {
            if (p[1] < p[i]) {
                star[starLen] = TEdge(p[1], p[i]);
                v[starLen][0] = w[1];
                v[starLen][1] = w[i];
            }
            else {
                star[starLen] = TEdge(p[i], p[1]);
                v[starLen][0] = w[i];
                v[starLen][1] = w[1];
            }
            ++starLen;
        }
        assert(starLen == 4);

        for (int i = 0; i < 4; ++i) {
            if (edges.find(star[i]) == edges.end()) {
                // Must be a diagonal edge
                assert(
                    star[i].point[0].ix != star[i].point[1].ix &&
                    star[i].point[0].iy != star[i].point[1].iy &&
                    star[i].point[0].iz != star[i].point[1].iz
                );

                R3Point vertexOnEdge;
                bool edgeHasVertex = star[i].computeVertex(
                    v[i][0], v[i][1],
                    tbox.origin,
                    tbox.dx, tbox.dy, tbox.dz,
                    vertexOnEdge
                );
                assert(edgeHasVertex);
                assert(tbox.contains(vertexOnEdge));

                // Add a vertex to triangulation
                triangulation.vertices.push_back(vertexOnEdge);
                vertexIdx[i] = (int)(triangulation.vertices.size() - 1);
                edges[star[i]] = vertexIdx[i];
            }
            else {
                vertexIdx[i] = edges[star[i]];
            }
        }

        /*...
        triangulation.triangles.push_back(
            Triangulation::Triangle(
                vertexIdx[0], vertexIdx[1], vertexIdx[3]
            )
        );
        triangulation.triangles.push_back(
            Triangulation::Triangle(
                vertexIdx[0], vertexIdx[3], vertexIdx[2]
            )
        );
        ...*/

        R3Point t0 = p[0].coord(
            tbox.origin,
            tbox.dx, tbox.dy, tbox.dz
        );
        R3Point t1 = triangulation.vertices[vertexIdx[0]].point;
        R3Point t2 = triangulation.vertices[vertexIdx[1]].point;
        R3Point t3 = triangulation.vertices[vertexIdx[2]].point;
        R3Point t4 = triangulation.vertices[vertexIdx[3]].point;

        // Select the shortest diagonal
        double diag03 = (t4 - t1).length2();
        double diag12 = (t3 - t2).length2();
        if (diag03 <= diag12) {
            // Select the 0--3 diagonal
            R3Vector n = (t2 - t1).vectorProduct(t4 - t1);
            double dotprod = n * (t0 - t1);

            if (
                (dotprod >= 0. && w[0] >= 0.) ||
                (dotprod <= 0. && w[0] <= 0.)
                ) {
                triangulation.triangles.push_back(
                    Triangulation::Triangle(
                        vertexIdx[0], vertexIdx[1], vertexIdx[3]
                    )
                );
                triangulation.triangles.push_back(
                    Triangulation::Triangle(
                        vertexIdx[0], vertexIdx[3], vertexIdx[2]
                    )
                );
            }
            else {
                triangulation.triangles.push_back(
                    Triangulation::Triangle(
                        vertexIdx[0], vertexIdx[3], vertexIdx[1]
                    )
                );
                triangulation.triangles.push_back(
                    Triangulation::Triangle(
                        vertexIdx[0], vertexIdx[2], vertexIdx[3]
                    )
                );
            }
        }
        else {
            // Select the 1--2 diagonal
            R3Vector n = (t2 - t1).vectorProduct(t3 - t1);
            double dotprod = n * (t0 - t1);

            if (
                (dotprod >= 0. && w[0] >= 0.) ||
                (dotprod <= 0. && w[0] <= 0.)
                ) {
                triangulation.triangles.push_back(
                    Triangulation::Triangle(
                        vertexIdx[0], vertexIdx[1], vertexIdx[2]
                    )
                );
                triangulation.triangles.push_back(
                    Triangulation::Triangle(
                        vertexIdx[1], vertexIdx[3], vertexIdx[2]
                    )
                );
            }
            else {
                triangulation.triangles.push_back(
                    Triangulation::Triangle(
                        vertexIdx[0], vertexIdx[2], vertexIdx[1]
                    )
                );
                triangulation.triangles.push_back(
                    Triangulation::Triangle(
                        vertexIdx[1], vertexIdx[2], vertexIdx[3]
                    )
                );
            }
        } // end if (diag...

        /*??????????????????????????????????????????????
        // Define the correct order of points
        R3Point quad[4];
        for (int i = 0; i < 4; ++i) {
            quad[i] =
                triangulation.vertices[i].point;
        }

        // Look for a diagonal
        R3Vector diag = quad[1] - quad[0];
        R3Vector e1 = diag.vectorProduct(quad[2] - quad[0]);
        R3Vector n = e1.vectorProduct(diag);
        double s0 = (quad[2] - quad[0])*n;
        double s1 = (quad[3] - quad[0])*n;
        if (
            (s0 > 0. && s1 < 0.) ||
            (s0 < 0. && s1 > 0.)
        ) {
            // Diagonal is found!
            triangulation.triangles.push_back(
                Triangulation::Triangle(
                    vertexIdx[0], vertexIdx[1], vertexIdx[2]
                )
            );
            triangulation.triangles.push_back(
                Triangulation::Triangle(
                    vertexIdx[0], vertexIdx[1], vertexIdx[3]
                )
            );
            return;
        }

        diag = quad[2] - quad[0];
        e1 = diag.vectorProduct(quad[1] - quad[0]);
        n = e1.vectorProduct(diag);
        s0 = (quad[1] - quad[0])*n;
        s1 = (quad[3] - quad[0])*n;
        if (
            (s0 > 0. && s1 < 0.) ||
            (s0 < 0. && s1 > 0.)
        ) {
            // Diagonal is found!
            triangulation.triangles.push_back(
                Triangulation::Triangle(
                    vertexIdx[0], vertexIdx[2], vertexIdx[1]
                )
            );
            triangulation.triangles.push_back(
                Triangulation::Triangle(
                    vertexIdx[0], vertexIdx[2], vertexIdx[3]
                )
            );
            return;
        }

        triangulation.triangles.push_back(
            Triangulation::Triangle(
                vertexIdx[0], vertexIdx[3], vertexIdx[1]
            )
        );
        triangulation.triangles.push_back(
            Triangulation::Triangle(
                vertexIdx[0], vertexIdx[3], vertexIdx[2]
            )
        );
        ??????????????????????????????????????????????*/

    }   // end if
}

R3Vector gradientR3(
    double (*f)(const R3Point&),
    const R3Point& p,
    double h
) {
    double h2 = h * 2.;
    double dx = (*f)(p + R3Vector(h, 0., 0.)) -
        (*f)(p - R3Vector(h, 0., 0.));
    dx /= h2;
    double dy = (*f)(p + R3Vector(0, h, 0.)) -
        (*f)(p - R3Vector(0., h, 0.));
    dy /= h2;
    double dz = (*f)(p + R3Vector(0., 0., h)) -
        (*f)(p - R3Vector(0., 0., h));
    dz /= h2;
    R3Vector n(dx, dy, dz);
    n.normalize();
    return n;
}
