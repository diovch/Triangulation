//#include <QFile>
//#include <QDataStream>
#include "voxelset.h"
#include "Triangulation.h"

#include <map>

using namespace R3Graph;

const int NEIGHBOURS6[6][3] = {         // slice, x, y
    { 0, -1, 0 },
    { 0, 1, 0 },
    { 0, 0, -1 },
    { 0, 0, 1 },    // First 4 element MUST BE on the same slice!!!
    //-------------------------------------------------------------
    { -1, 0, 0 },
    { 1, 0, 0 }
};

const int NEIGHBOURS4_3D[4][3] = {      // slice, x, y
    { 0, -1, 0 },
    { 0, 1, 0 },
    { 0, 0, -1 },
    { 0, 0, 1 }
};

const int Voxel::FACE_DIRECTIONS[6][3] = {
    {  0, -1,  0 }, // FACE_LEFT    slice, x, y
    {  0,  1,  0 }, // FACE_RIGHT
    {  0,  0, -1 }, // FACE_FRONT
    {  0,  0,  1 }, // FACE_BACK
    { -1,  0,  0 }, // FACE_BOTTOM
    {  1,  0,  0 }  // FACE_TOP
};

void VoxelSet::initialize(
   int x_max, int y_max, int max_slices
) {
    xMax = x_max;;
    yMax = y_max;
    maxSlices = max_slices;
    bitmasks.resize(maxSlices);
    numVoxels = 0.;
    for (int i = 0; i < maxSlices; ++i) {
        bitmasks.at(i).resizeBitmask(x_max, y_max);
        bitmasks.at(i).clear();
    }
    detected3D = false;
}

void VoxelSet::clear() {
    for (int i = 0; i < maxSlices; ++i) {
        bitmasks.at(i).clear();
    }
    numVoxels = 0.;
}

void VoxelSet::computeVoxelBox() {
    numVoxels = 0.;
    for (int slice = 0; slice < maxSlices; ++slice) {
        const Bitmask& bitmask = bitmaskAt(slice);
        for (int y = 0; y < yMax; ++y) {
            for (int x = 0; x < xMax; ++x) {
                if (bitmask.pixelAt(x, y) != 0) {
                    if (numVoxels == 0) {
                        // Initialize a voxel box
                        voxelBox.origin = Voxel(slice, x, y);
                        voxelBox.width = 0;
                        voxelBox.depth = 0;
                        voxelBox.height = 0;
                    } else {
                        if (x < voxelBox.xMin()) {
                            voxelBox.width += voxelBox.xMin() - x;
                            voxelBox.setXMin(x);
                        } else if (x > voxelBox.xMax()) {
                            voxelBox.width += x - voxelBox.xMax();
                        }

                        if (y < voxelBox.yMin()) {
                            voxelBox.depth += voxelBox.yMin() - y;
                            voxelBox.setYMin(y);
                        } else if (y > voxelBox.yMax()) {
                            voxelBox.depth += y - voxelBox.yMax();
                        }

                        if (slice < voxelBox.sliceMin()) {
                            voxelBox.height += voxelBox.sliceMin() - slice;
                            voxelBox.setSliceMin(slice);
                        } else if (slice > voxelBox.sliceMax()) {
                            voxelBox.height += slice - voxelBox.sliceMax();
                        }
                    }
                    numVoxels += 1.;
                } // end if
            } // end for (x...
        } // end for (y...
    } // end for (slice...
}

//bool VoxelSet::save(QString path) const {
//    QFile f(path);
//    if (!f.open(QIODevice::WriteOnly))
//        return false;
//    QDataStream out(&f);
//
//    out << threshold;
//    out << xMax;
//    out << yMax;
//    out << maxSlices;
//    if (out.status() == QDataStream::WriteFailed)
//        return false;
//    for (int s = 0; s < maxSlices; ++s) {
//        const Bitmask& bitmask = bitmaskAt(s);
//        if (!bitmask.write(out))
//            return false;
//    }
//    return true;
//}
//
//bool VoxelSet::load(QString path) {
//    QFile f(path);
//    if (!f.open(QIODevice::ReadOnly))
//        return false;
//    QDataStream in(&f);
//    in >> threshold;
//    in >> xMax;
//    in >> yMax;
//    in >> maxSlices;
//    if (in.status() != QDataStream::Ok)
//        return false;
//    for (int s = 0; s < maxSlices; ++s) {
//        Bitmask& bitmask = bitmaskAt(s);
//        if (!bitmask.read(in))
//            return false;
//    }
//    return true;
//}

PackedVoxelSet::PackedVoxelSet(const VoxelSet& vs) {
    operator=(vs);
}

PackedVoxelSet& PackedVoxelSet::operator=(const VoxelSet& vs) {
    xMax = vs.xMax;
    yMax = vs.yMax;
    maxSlices = vs.maxSlices;
    numVoxels = vs.numVoxels;
    voxelBox = vs.voxelBox;
    threshold = vs.threshold;
    detected3D = vs.detected3D;

    packedBitmasks.resize(vs.bitmasks.size());
    assert(maxSlices == (int) vs.bitmasks.size());
    for (int slice = 0; slice < maxSlices; ++slice) {
        packedBitmasks.at(slice) = vs.bitmasks.at(slice);
    }

    return *this;
}

VoxelSet& PackedVoxelSet::unpack(VoxelSet& vs) const {
    vs.xMax = xMax;
    vs.yMax = yMax;
    vs.maxSlices = maxSlices;
    vs.numVoxels = numVoxels;
    vs.voxelBox = voxelBox;
    vs.threshold = threshold;
    vs.detected3D = detected3D;

    vs.bitmasks.resize(packedBitmasks.size());
    assert(maxSlices == (int) packedBitmasks.size());
    for (int slice = 0; slice < maxSlices; ++slice) {
        packedBitmasks.at(slice).unpack(vs.bitmasks.at(slice));
    }

    return vs;
}

bool VoxelSet::faceOpen(const Voxel& v, int face) const {
    assert(voxelAt(v) != 0);
    Voxel neighbour(
        v.slice + Voxel::FACE_DIRECTIONS[face][0],
        v.point.x + Voxel::FACE_DIRECTIONS[face][1],
        v.point.y + Voxel::FACE_DIRECTIONS[face][2]
    );
    return (
        !voxelInVolume(neighbour) ||
        voxelAt(neighbour) == 0
    );
}

void computeTriangulationOfVoxelSet(
    Triangulation& triangulation,
    const VoxelSet& voxelSet,
    const R3Point& origin,
    double dx, double dy, double dz
) {
    const VoxelBox& voxelBox = voxelSet.voxelBox;
    int sliceStart = voxelBox.origin.slice - 2;
    if (sliceStart < 0)
        sliceStart = 0;

    int sliceFinish = voxelBox.origin.slice + voxelBox.height + 1;
    if (sliceFinish >= voxelSet.maxSlices)
        sliceFinish = voxelSet.maxSlices - 1;

    int ixmin = voxelBox.origin.point.x - 2;
    if (ixmin < 0)
        ixmin = 0;

    int ixmax = voxelBox.origin.point.x + voxelBox.width + 1;
    if (ixmax >= voxelSet.xMax)
        ixmax = voxelSet.xMax - 1;

    int iymin = voxelBox.origin.point.y - 2;
    if (iymin < 0)
        iymin = 0;

    int iymax = voxelBox.origin.point.y + voxelBox.depth + 1;
    if (iymax >= voxelSet.yMax)
        iymax = voxelSet.yMax - 1;

    triangulation.clear();

    // Indices of extended voxels vertices in array
    // Each voxel produce 8 vertices == extended voxels
    // voxel (slice, x, y) -> 8 extended voxels:
    //       (2*slice + s0, 2*x + s1, 2*y + s2), where
    //                                           si = +-1
    std::map<Voxel, int> vertexIndices;

    for (int slice = sliceStart; slice <= sliceFinish; ++slice) {
        for (int iy = iymin; iy <= iymax; ++iy) {
            for (int ix = ixmin; ix <= ixmax; ++ix) {
                if (voxelSet.voxelAt(slice, ix, iy) == 0)
                    continue;

                // Enumeration of cube vertices and faces:
                //        7         6
                //       +---------+          z
                //      /   top   /|        ^
                //   4 / |       / |        |
                //    +---------+5 |        |
                //Left|         |  | Right  |    ^ y
                //    |  |  back|  |        |   /
                //    |  + - - -|- +        |  /
                //    |  3      | /2        | /
                //    |/  Front |/          |/
                //    +---------+           ------> x
                //   0  bottom   1

                R3Point cubeVertices[8];
                //... R3Vector cubeGradients[8];

                Voxel cube(slice, ix, iy);
                R3Point cubeCenter = voxel3DCoord(
                    cube, origin, dx, dy, dz
                );

                // Bottom
                Voxel neighborVoxel =
                    Voxel(cube.slice - 1, cube.point.x, cube.point.y);
                R3Point bottomCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                // Top
                neighborVoxel.slice += 2;
                R3Point topCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                // Left
                neighborVoxel =
                    Voxel(cube.slice, cube.point.x - 1, cube.point.y);
                R3Point leftCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                // Right
                neighborVoxel.point.x += 2;
                R3Point rightCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                // Front
                neighborVoxel =
                    Voxel(cube.slice, cube.point.x, cube.point.y - 1);
                R3Point frontCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                // Back
                neighborVoxel.point.y += 2;
                R3Point backCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                cubeVertices[0] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;
                cubeVertices[1] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;
                cubeVertices[2] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;
                cubeVertices[3] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;

                cubeVertices[4] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;
                cubeVertices[5] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;
                cubeVertices[6] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;
                cubeVertices[7] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;

                Voxel extendedVoxels[8];
                for (int iv = 0; iv < 8; ++iv) {
                    extendedVoxels[iv] = Voxel(
                        cube.slice*2, cube.point.x*2, cube.point.y*2
                    );
                }
                --(extendedVoxels[0].slice);
                --(extendedVoxels[1].slice);
                --(extendedVoxels[2].slice);
                --(extendedVoxels[3].slice);
                ++(extendedVoxels[4].slice);
                ++(extendedVoxels[5].slice);
                ++(extendedVoxels[6].slice);
                ++(extendedVoxels[7].slice);

                --(extendedVoxels[0].point.x);
                --(extendedVoxels[3].point.x);
                --(extendedVoxels[4].point.x);
                --(extendedVoxels[7].point.x);
                ++(extendedVoxels[1].point.x);
                ++(extendedVoxels[2].point.x);
                ++(extendedVoxels[5].point.x);
                ++(extendedVoxels[6].point.x);

                --(extendedVoxels[0].point.y);
                --(extendedVoxels[1].point.y);
                --(extendedVoxels[4].point.y);
                --(extendedVoxels[5].point.y);
                ++(extendedVoxels[2].point.y);
                ++(extendedVoxels[3].point.y);
                ++(extendedVoxels[6].point.y);
                ++(extendedVoxels[7].point.y);

                int indices[8]; // Indices of vertices in array
                for (int iv = 0; iv < 8; ++iv)
                    indices[iv] = (-1);

                // Front face
                if (voxelSet.faceOpen(cube, Voxel::FACE_FRONT)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[0]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(
                            cubeVertices[0]
                        );
                        indices[0] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[0]] = indices[0];
                    } else {
                        // Point is already in the array
                        indices[0] = vertexIndices[extendedVoxels[0]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[1]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[1]);
                        indices[1] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[1]] = indices[1];
                    } else {
                        // Point is already in the array
                        indices[1] = vertexIndices[extendedVoxels[1]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[4]
                        ) == 0
                    )
                    {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[4]);
                        indices[4] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[4]] = indices[4];
                    } else {
                        // Point is already in the array
                        indices[4] = vertexIndices[extendedVoxels[4]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[5]
                        ) == 0
                    )
                    {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[5]);
                        indices[5] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[5]] = indices[5];
                    } else {
                        // Point is already in the array
                        indices[5] = vertexIndices[extendedVoxels[5]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[0] >= 0 &&
                        indices[1] >= 0 &&
                        indices[4] >= 0 &&
                        indices[5] >= 0
                    );

                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[0], indices[1], indices[5]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[0], indices[5], indices[4]
                        )
                    );
                } // end if (... FRONT_FACE

                // Back face
                if (voxelSet.faceOpen(cube, Voxel::FACE_BACK)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[2]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(
                            cubeVertices[2]
                        );
                        indices[2] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[2]] = indices[2];
                    } else {
                        // Point is already in the array
                        indices[2] = vertexIndices[extendedVoxels[2]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[3]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(
                            cubeVertices[3]
                        );
                        indices[3] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[3]] = indices[3];
                    } else {
                        // Point is already in the array
                        indices[3] = vertexIndices[extendedVoxels[3]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[6]
                        ) == 0
                    )
                    {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(
                            cubeVertices[6]
                        );
                        indices[6] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[6]] = indices[6];
                    } else {
                        // Point is already in the array
                        indices[6] = vertexIndices[extendedVoxels[6]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[7]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(
                            cubeVertices[7]
                        );
                        indices[7] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[7]] = indices[7];
                    } else {
                        // Point is already in the array
                        indices[7] = vertexIndices[extendedVoxels[7]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[2] >= 0 &&
                        indices[3] >= 0 &&
                        indices[6] >= 0 &&
                        indices[7] >= 0
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[2], indices[3], indices[6]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[6], indices[3], indices[7]
                        )
                    );
                } // end if (... BACK_FACE

                // Left face
                if (voxelSet.faceOpen(cube, Voxel::FACE_LEFT)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[0]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[0]);
                        indices[0] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[0]] = indices[0];
                    } else {
                        // Point is already in the array
                        indices[0] = vertexIndices[extendedVoxels[0]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[3]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[3]);
                        indices[3] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[3]] = indices[3];
                    } else {
                        // Point is already in the array
                        indices[3] = vertexIndices[extendedVoxels[3]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[4]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[4]);
                        indices[4] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[4]] = indices[4];
                    } else {
                        // Point is already in the array
                        indices[4] = vertexIndices[extendedVoxels[4]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[7]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[7]);
                        indices[7] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[7]] = indices[7];
                    } else {
                        // Point is already in the array
                        indices[7] = vertexIndices[extendedVoxels[7]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[0] >= 0 &&
                        indices[3] >= 0 &&
                        indices[4] >= 0 &&
                        indices[7] >= 0
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[0], indices[7], indices[3]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[7], indices[0], indices[4]
                        )
                    );
                } // end if (... LEFT_FACE

                // Right face
                if (voxelSet.faceOpen(cube, Voxel::FACE_RIGHT)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[1]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[1]);
                        indices[1] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[1]] = indices[1];
                    } else {
                        // Point is already in the array
                        indices[1] = vertexIndices[extendedVoxels[1]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[2]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[2]);
                        indices[2] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[2]] = indices[2];
                    } else {
                        // Point is already in the array
                        indices[2] = vertexIndices[extendedVoxels[2]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[5]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[5]);
                        indices[5] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[5]] = indices[5];
                    } else {
                        // Point is already in the array
                        indices[5] = vertexIndices[extendedVoxels[5]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[6]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[6]);
                        indices[6] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[6]] = indices[6];
                    } else {
                        // Point is already in the array
                        indices[6] = vertexIndices[extendedVoxels[6]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[1] >= 0 &&
                        indices[2] >= 0 &&
                        indices[5] >= 0 &&
                        indices[6] >= 0
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[1], indices[2], indices[6]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[6], indices[5], indices[1]
                        )
                    );
                } // end if (... RIGHT_FACE

                // Bottom face
                if (voxelSet.faceOpen(cube, Voxel::FACE_BOTTOM)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[0]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[0]);
                        indices[0] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[0]] = indices[0];
                    } else {
                        // Point is already in the array
                        indices[0] = vertexIndices[extendedVoxels[0]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[1]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[1]);
                        indices[1] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[1]] = indices[1];
                    } else {
                        // Point is already in the array
                        indices[1] = vertexIndices[extendedVoxels[1]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[2]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[2]);
                        indices[2] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[2]] = indices[2];
                    } else {
                        // Point is already in the array
                        indices[2] = vertexIndices[extendedVoxels[2]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[3]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[3]);
                        indices[3] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[3]] = indices[3];
                    } else {
                        // Point is already in the array
                        indices[3] = vertexIndices[extendedVoxels[3]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[0] >= 0 &&
                        indices[1] >= 0 &&
                        indices[2] >= 0 &&
                        indices[3] >= 0
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[0], indices[2], indices[1]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[0], indices[3], indices[2]
                        )
                    );
                } // end if (... BOTTOM_FACE

                // Top face
                if (voxelSet.faceOpen(cube, Voxel::FACE_TOP)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[4]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(
                            cubeVertices[4]
                        );
                        indices[4] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[4]] = indices[4];
                    } else {
                        // Point is already in the array
                        indices[4] = vertexIndices[extendedVoxels[4]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[5]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(
                            cubeVertices[5]
                        );
                        indices[5] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[5]] = indices[5];
                    } else {
                        // Point is already in the array
                        indices[5] = vertexIndices[extendedVoxels[5]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[6]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(
                            cubeVertices[6]
                        );
                        indices[6] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[6]] = indices[6];
                    } else {
                        // Point is already in the array
                        indices[6] = vertexIndices[extendedVoxels[6]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[7]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(
                            cubeVertices[7]
                        );
                        indices[7] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[7]] = indices[7];
                    } else {
                        // Point is already in the array
                        indices[7] = vertexIndices[extendedVoxels[7]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[4] >= 0 &&
                        indices[5] >= 0 &&
                        indices[6] >= 0 &&
                        indices[7] >= 0
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[4], indices[5], indices[6]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[4], indices[6], indices[7]
                        )
                    );
                } // end if (... BOTTOM_FACE
            } // end for (ix...
        } // end for (iy...
    } // end for (slice...
}

void computeTriangulationOfVoxelSet_MY(
    short* pointer,
    int threshold,
    Triangulation& triangulation,
    const VoxelSet& voxelSet,
    const R3Point& origin,
    double dx, double dy, double dz
) {
    const VoxelBox& voxelBox = voxelSet.voxelBox;
    int sliceStart = voxelBox.origin.slice - 2;
    if (sliceStart < 0)
        sliceStart = 1;

    int sliceFinish = voxelBox.origin.slice + voxelBox.height + 1;
    if (sliceFinish >= voxelSet.maxSlices)
        sliceFinish = voxelSet.maxSlices - 1;

    int ixmin = voxelBox.origin.point.x - 2;
    if (ixmin < 0)
        ixmin = 1;

    int ixmax = voxelBox.origin.point.x + voxelBox.width + 1;
    if (ixmax >= voxelSet.xMax)
        ixmax = voxelSet.xMax - 1;

    int iymin = voxelBox.origin.point.y - 2;
    if (iymin < 0)
        iymin = 1;

    int iymax = voxelBox.origin.point.y + voxelBox.depth + 1;
    if (iymax >= voxelSet.yMax)
        iymax = voxelSet.yMax - 1;

    triangulation.clear();

    Voxel ImageCentralVoxel(voxelSet.maxSlices / 2, voxelSet.xMax / 2, voxelSet.yMax / 2);
    R3Point ImageCenter = voxel3DCoord(ImageCentralVoxel, origin, dx, dy, dz); 
    triangulation.imageCenter = ImageCenter;

    //for (int slice = sliceStart; slice <= sliceFinish; ++slice) {
    //    for (int iy = iymin; iy <= iymax; ++iy) {
    //        for (int ix = ixmin; ix <= ixmax; ++ix) {
    for (int slice = sliceFinish/2; slice <= sliceFinish*3/4; ++slice) {
        for (int iy = iymax/2; iy <= iymax * 3 / 4; ++iy) {
            for (int ix = ixmax/2; ix <= ixmax * 3 / 4; ++ix) {
                if (voxelSet.voxelAt(slice, ix, iy) == 0) // if(voxel is out of ROI)
                    continue;
                //if (triangulation.triangles.size() > 2'000'000)
                //    return;
                // Enumeration of cube vertices and faces:
                //        7         6
                //       +---------+          z
                //      /   top   /|        ^
                //   4 / |       / |        |
                //    +---------+5 |        |
                //Left|         |  | Right  |    ^ y
                //    |  |  back|  |        |   /
                //    |  + - - -|- +        |  /
                //    |  3      | /2        | /
                //    |/  Front |/          |/
                //    +---------+           ------> x
                //   0  bottom   1

                R3Point cubeVertices[8];

                Voxel cube(slice, ix, iy);
                R3Point cubeCenter = voxel3DCoord(
                    cube, origin, dx, dy, dz
                );

                // Bottom
                Voxel BottomVoxel =
                    Voxel(cube.slice - 1, cube.point.x, cube.point.y);
                R3Point bottomCenter = voxel3DCoord(
                    BottomVoxel,
                    origin, dx, dy, dz
                );

                // Top
                Voxel TopVoxel =
                    Voxel(cube.slice + 1, cube.point.x, cube.point.y);
                R3Point topCenter = voxel3DCoord(
                    TopVoxel,
                    origin, dx, dy, dz
                );

                // Left
                Voxel LeftVoxel =
                    Voxel(cube.slice, cube.point.x - 1, cube.point.y);
                R3Point leftCenter = voxel3DCoord(
                    LeftVoxel,
                    origin, dx, dy, dz
                );

                // Right
                Voxel RightVoxel =
                    Voxel(cube.slice, cube.point.x + 1, cube.point.y);
                R3Point rightCenter = voxel3DCoord(
                    RightVoxel,
                    origin, dx, dy, dz
                );

                // Front
                Voxel FrontVoxel =
                    Voxel(cube.slice, cube.point.x, cube.point.y - 1);
                R3Point frontCenter = voxel3DCoord(
                    FrontVoxel,
                    origin, dx, dy, dz
                );

                // Back
                Voxel BackVoxel =
                    Voxel(cube.slice, cube.point.x, cube.point.y + 1);
                R3Point backCenter = voxel3DCoord(
                    BackVoxel,
                    origin, dx, dy, dz
                );



                cubeVertices[0] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;
                cubeVertices[1] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;
                cubeVertices[2] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;
                cubeVertices[3] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;

                cubeVertices[4] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;
                cubeVertices[5] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;
                cubeVertices[6] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;
                cubeVertices[7] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;

                R3Vector cubeGradient = cube.VoxelGradient(pointer, dx, dy, dz);// HU per mm
                short cubeDensity = cube.VoxelDensity(pointer);
                cubeDensity -= threshold; // new level of zero density

                // ThresholdFunction(x,y,z) = Density(x,y,z) - threshold
                // f(x,y,z) = f(x0,y0,z0) + {x-x0, y-y0, z-z0} * grad f
                double ThresholdFunction[8];
                for (int i = 0; i < 8; ++i) // for cube vertices
                    ThresholdFunction[i] = cubeDensity + (cubeVertices[i] - cubeVertices[0]).scalarProduct(cubeGradient);
                // for cube center
                double ThresholdFunctionCubeCenter = cubeDensity + (cubeCenter - cubeVertices[0]).scalarProduct(cubeGradient);

                std::pair<R3Point, double> CubeCenterPair(cubeCenter, ThresholdFunctionCubeCenter);
                std::pair<R3Point, double> CubeVertexPairs[8];
                for (int i = 0; i < 8; ++i)
                    CubeVertexPairs[i] = { cubeVertices[i], ThresholdFunction[i] };



                // Front face
                // if (front_neihbour is out of ROI or out of maxVolume)
                if (voxelSet.faceOpen(cube, Voxel::FACE_FRONT)) {

                    R3Vector FrontVoxelGradient =
                        FrontVoxel.VoxelGradient(pointer, dx, dy, dz); // HU per mm
                    short FrontDensity = FrontVoxel.VoxelDensity(pointer);
                    FrontDensity -= threshold;

                    double ThresholdFunctionFrontCenter = FrontDensity +
                        (cubeCenter - cubeVertices[0]).scalarProduct(FrontVoxelGradient);

                    // process tetrahedrons with 0, 1, 5, 4, cubeCenter, frontCenter vertecies
                    std::pair<R3Point, double> FrontCenterPair(frontCenter, ThresholdFunctionFrontCenter);

                    std::pair<R3Point, double> AdjacentCubeVertexPairs[4];
                    AdjacentCubeVertexPairs[0] = CubeVertexPairs[0];
                    AdjacentCubeVertexPairs[1] = CubeVertexPairs[1];
                    AdjacentCubeVertexPairs[2] = CubeVertexPairs[5];
                    AdjacentCubeVertexPairs[3] = CubeVertexPairs[4];

                    for (int i = 0; i < 4; ++i) {
                        // cube edges: 0-1, 1-2, 2-3, (!) 3-0 (!)
                        Tetrahedron CurrentTetrahedron(CubeCenterPair, FrontCenterPair,
                            AdjacentCubeVertexPairs[i], AdjacentCubeVertexPairs[(i + 1) % 4]);
                        triangulation.TriangulationOfTetrahedron(CurrentTetrahedron);
                    }
                } // end if (... FRONT_FACE

                // Back face
                // if (back_neihbour is out of ROI or out of maxVolume)
                if (voxelSet.faceOpen(cube, Voxel::FACE_BACK)) {

                    R3Vector BackVoxelGradient =
                        BackVoxel.VoxelGradient(pointer, dx, dy, dz); // HU per mm
                    short BackDensity = BackVoxel.VoxelDensity(pointer);
                    BackDensity -= threshold;

                    double ThresholdFunctionBackCenter = BackDensity +
                        (cubeCenter - cubeVertices[0]).scalarProduct(BackVoxelGradient);

                    // process tetrahedrons with 2, 6, 7, 3, cubeCenter, backCenter vertecies
                    std::pair<R3Point, double> BackCenterPair(backCenter, ThresholdFunctionBackCenter);

                    std::pair<R3Point, double> AdjacentCubeVertexPairs[4];
                    AdjacentCubeVertexPairs[0] = CubeVertexPairs[2];
                    AdjacentCubeVertexPairs[1] = CubeVertexPairs[6];
                    AdjacentCubeVertexPairs[2] = CubeVertexPairs[7];
                    AdjacentCubeVertexPairs[3] = CubeVertexPairs[3];

                    for (int i = 0; i < 4; ++i) {
                        Tetrahedron CurrentTetrahedron(CubeCenterPair, BackCenterPair,
                            AdjacentCubeVertexPairs[i], AdjacentCubeVertexPairs[(i + 1) % 4]);
                        triangulation.TriangulationOfTetrahedron(CurrentTetrahedron);
                    }

                } // end if (... BACK_FACE

                // Left face
                // if (left_neihbour is out of ROI or out of maxVolume)
                if (voxelSet.faceOpen(cube, Voxel::FACE_LEFT)) {

                    R3Vector LeftVoxelGradient =
                        LeftVoxel.VoxelGradient(pointer, dx, dy, dz); // HU per mm
                    short LeftDensity = LeftVoxel.VoxelDensity(pointer);
                    LeftDensity -= threshold;

                    double ThresholdFunctionLeftCenter = LeftDensity +
                        (cubeCenter - cubeVertices[0]).scalarProduct(LeftVoxelGradient);

                    // process tetrahedrons with 0, 3, 7, 4, cubeCenter, leftCenter vertecies
                    std::pair<R3Point, double> LeftCenterPair(leftCenter, ThresholdFunctionLeftCenter);

                    std::pair<R3Point, double> AdjacentCubeVertexPairs[4];
                    AdjacentCubeVertexPairs[0] = CubeVertexPairs[0];
                    AdjacentCubeVertexPairs[1] = CubeVertexPairs[3];
                    AdjacentCubeVertexPairs[2] = CubeVertexPairs[7];
                    AdjacentCubeVertexPairs[3] = CubeVertexPairs[4];

                    for (int i = 0; i < 4; ++i) {
                        Tetrahedron CurrentTetrahedron(CubeCenterPair, LeftCenterPair,
                            AdjacentCubeVertexPairs[i], AdjacentCubeVertexPairs[(i + 1) % 4]);
                        triangulation.TriangulationOfTetrahedron(CurrentTetrahedron);
                    }

                } // end if (... LEFT_FACE

                // Right face
                // if (right_neihbour is out of ROI or out of maxVolume)
                if (voxelSet.faceOpen(cube, Voxel::FACE_RIGHT)) {

                    R3Vector RightVoxelGradient =
                        RightVoxel.VoxelGradient(pointer, dx, dy, dz); // HU per mm
                    short RightDensity = RightVoxel.VoxelDensity(pointer);
                    RightDensity -= threshold;

                    double ThresholdFunctionRightCenter = RightDensity +
                        (cubeCenter - cubeVertices[0]).scalarProduct(RightVoxelGradient);

                    // process tetrahedrons with 1, 2, 6, 5, cubeCenter, rightCenter vertecies
                    std::pair<R3Point, double> RightCenterPair(rightCenter, ThresholdFunctionRightCenter);

                    std::pair<R3Point, double> AdjacentCubeVertexPairs[4];
                    AdjacentCubeVertexPairs[0] = CubeVertexPairs[1];
                    AdjacentCubeVertexPairs[1] = CubeVertexPairs[2];
                    AdjacentCubeVertexPairs[2] = CubeVertexPairs[6];
                    AdjacentCubeVertexPairs[3] = CubeVertexPairs[5];

                    for (int i = 0; i < 4; ++i) {
                        Tetrahedron CurrentTetrahedron(CubeCenterPair, RightCenterPair,
                            AdjacentCubeVertexPairs[i], AdjacentCubeVertexPairs[(i + 1) % 4]);
                        triangulation.TriangulationOfTetrahedron(CurrentTetrahedron);
                    }

                } // end if (... RIGHT_FACE

                // Bottom face
                // if (bottom_neihbour is out of ROI or out of maxVolume)
                if (voxelSet.faceOpen(cube, Voxel::FACE_BOTTOM)) {

                    R3Vector BottomVoxelGradient =
                        BottomVoxel.VoxelGradient(pointer, dx, dy, dz); // HU per mm
                    short BottomDensity = BottomVoxel.VoxelDensity(pointer);
                    BottomDensity -= threshold;

                    double ThresholdFunctionBottomCenter = BottomDensity +
                        (cubeCenter - cubeVertices[0]).scalarProduct(BottomVoxelGradient);

                    // process tetrahedrons with 0, 1, 2, 3, cubeCenter, bottomCenter vertecies
                    std::pair<R3Point, double> BottomCenterPair(bottomCenter, ThresholdFunctionBottomCenter);

                    std::pair<R3Point, double> AdjacentCubeVertexPairs[4];
                    AdjacentCubeVertexPairs[0] = CubeVertexPairs[0];
                    AdjacentCubeVertexPairs[1] = CubeVertexPairs[1];
                    AdjacentCubeVertexPairs[2] = CubeVertexPairs[2];
                    AdjacentCubeVertexPairs[3] = CubeVertexPairs[3];

                    for (int i = 0; i < 4; ++i) {
                        Tetrahedron CurrentTetrahedron(CubeCenterPair, BottomCenterPair,
                            AdjacentCubeVertexPairs[i], AdjacentCubeVertexPairs[(i + 1) % 4]);
                        triangulation.TriangulationOfTetrahedron(CurrentTetrahedron);
                    }

                } // end if (... BOTTOM_FACE

                // Top face
                // if (top_neihbour is out of ROI or out of maxVolume)
                if (voxelSet.faceOpen(cube, Voxel::FACE_TOP)) {

                    R3Vector TopVoxelGradient =
                        TopVoxel.VoxelGradient(pointer, dx, dy, dz); // HU per mm
                    short TopDensity = TopVoxel.VoxelDensity(pointer);
                    TopDensity -= threshold;

                    double ThresholdFunctionTopCenter = TopDensity +
                        (cubeCenter - cubeVertices[0]).scalarProduct(TopVoxelGradient);

                    // process tetrahedrons with 4, 5, 6, 7, cubeCenter, topCenter vertecies
                    std::pair<R3Point, double> TopCenterPair(topCenter, ThresholdFunctionTopCenter);

                    std::pair<R3Point, double> AdjacentCubeVertexPairs[4];
                    AdjacentCubeVertexPairs[0] = CubeVertexPairs[4];
                    AdjacentCubeVertexPairs[1] = CubeVertexPairs[5];
                    AdjacentCubeVertexPairs[2] = CubeVertexPairs[6];
                    AdjacentCubeVertexPairs[3] = CubeVertexPairs[7];

                    for (int i = 0; i < 4; ++i) {
                        Tetrahedron CurrentTetrahedron(CubeCenterPair, TopCenterPair,
                            AdjacentCubeVertexPairs[i], AdjacentCubeVertexPairs[(i + 1) % 4]);
                        triangulation.TriangulationOfTetrahedron(CurrentTetrahedron);
                    }

                } // end if (... TOP_FACE

            } // end for (ix...
        } // end for (iy...
    } // end for (slice...
    //triangulation.orientate();
}

