//#include <QFile>
//#include <QDataStream>
#include "voxelset.h"
#include "Triangulation.h"

#include <map>
#include <deque>

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
    std::map<int, std::set<int>>& VerxteNeighbours,
    Triangulation& triangulation,
    const VoxelSet& voxelSet,
    const R3Point& origin,
    double dx, double dy, double dz
) {
    int sliceStart = 0, sliceFinish = 0, ixmin = 0, ixmax = 0, iymin = 0, iymax = 0;
    BoxBorders(voxelSet, sliceStart, sliceFinish, ixmin, ixmax, iymin, iymax);

    triangulation.clear();
    
    // Indices of extended voxels vertices in array
    // Each voxel produce 8 vertices == extended voxels
    // voxel (slice, x, y) -> 8 extended voxels:
    //       (2*slice + s0, 2*x + s1, 2*y + s2), where
    //                                           si = +-1
    std::map<Voxel, int> vertexIndices;
    std::vector<int> ind;
    //for (int slice = sliceStart; slice <= sliceFinish; ++slice) {
    //    for (int iy = iymin; iy <= iymax; ++iy) {
    //        for (int ix = ixmin; ix <= ixmax; ++ix) {
    for (int slice = sliceFinish / 2; slice <= sliceFinish * 3 / 4; ++slice) {
        for (int iy = iymax / 2; iy <= iymax * 3 / 4; ++iy) {
            for (int ix = ixmax / 2; ix <= ixmax * 3 / 4; ++ix) {
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

                const Voxel cube(slice, ix, iy);
                R3Point cubeCenter, bottomCenter, topCenter, leftCenter, rightCenter, frontCenter, backCenter;
                InitializeNeighboursCentres(cube, origin, dx, dy,dz,
                    cubeCenter, bottomCenter, topCenter, leftCenter, rightCenter, frontCenter, backCenter);

                R3Point cubeVertices[8];
                InitializeCubeVerticies(cubeVertices,
                    cubeCenter, bottomCenter, topCenter, leftCenter, rightCenter, frontCenter, backCenter);

                Voxel extendedVoxels[8];
                InitializeExtendedVoxels(extendedVoxels, cube);

                int indices[8]; // Indices of vertices in array
                for (int iv = 0; iv < 8; ++iv)
                    indices[iv] = (-1);

                const std::vector<Voxel::Face> faces = { Voxel::FACE_FRONT,
                        Voxel::FACE_BACK, Voxel::FACE_LEFT, Voxel::FACE_RIGHT, Voxel::FACE_BOTTOM, Voxel::FACE_TOP };

                for (const auto& face : faces)
                {
                    FaceTriangulation(voxelSet, cube, face,
                        vertexIndices, extendedVoxels, triangulation,
                        indices, cubeVertices, VerxteNeighbours);
                }

            }
        }
    }
}

void FillNeighbours(std::map<int, std::set<int>>& neighbours, std::vector<int>& indicies)
{
    for (int i = 0; i < indicies.size(); ++i)
    {
        neighbours[indicies[i]].insert(indicies[(i + 1) % 3]);
        neighbours[indicies[i]].insert(indicies[(i + 2) % 3]);
    }
    std::fill(indicies.begin(), indicies.end(), 0);
}

void BoxBorders(const VoxelSet& voxelSet,
    int& sliceStart, int& sliceFinish, int& ixmin, int& ixmax, int& iymin, int& iymax)
{
    const VoxelBox& voxelBox = voxelSet.voxelBox;
    sliceStart = voxelBox.origin.slice - 2;
    if (sliceStart < 0)
        sliceStart = 0;

    sliceFinish = voxelBox.origin.slice + voxelBox.height + 1;
    if (sliceFinish >= voxelSet.maxSlices)
        sliceFinish = voxelSet.maxSlices - 1;

    ixmin = voxelBox.origin.point.x - 2;
    if (ixmin < 0)
        ixmin = 0;

    ixmax = voxelBox.origin.point.x + voxelBox.width + 1;
    if (ixmax >= voxelSet.xMax)
        ixmax = voxelSet.xMax - 1;

    iymin = voxelBox.origin.point.y - 2;
    if (iymin < 0)
        iymin = 0;

    iymax = voxelBox.origin.point.y + voxelBox.depth + 1;
    if (iymax >= voxelSet.yMax)
        iymax = voxelSet.yMax - 1;
}

void InitializeNeighboursCentres(const Voxel& cube, const R3Point& origin, double dx, double dy, double dz, R3Point& cubeCenter,
    R3Point& bottomCenter, R3Point& topCenter,
    R3Point& leftCenter, R3Point& rightCenter,
    R3Point& frontCenter, R3Point& backCenter)
{
    cubeCenter = voxel3DCoord(cube, origin, dx, dy, dz);

    Voxel neighborVoxel = Voxel(cube.slice - 1, cube.point.x, cube.point.y);
    bottomCenter = voxel3DCoord(
        neighborVoxel,
        origin, dx, dy, dz);

    neighborVoxel.slice += 2;
    topCenter = voxel3DCoord(
        neighborVoxel,
        origin, dx, dy, dz);

    neighborVoxel = Voxel(cube.slice, cube.point.x - 1, cube.point.y);
    leftCenter = voxel3DCoord( 
        neighborVoxel,
        origin, dx, dy, dz);

    neighborVoxel.point.x += 2;
    rightCenter = voxel3DCoord(
        neighborVoxel,
        origin, dx, dy, dz);

    neighborVoxel = Voxel(cube.slice, cube.point.x, cube.point.y - 1);
    frontCenter = voxel3DCoord(
        neighborVoxel,
        origin, dx, dy, dz);

    neighborVoxel.point.y += 2;
    backCenter = voxel3DCoord(
        neighborVoxel,
        origin, dx, dy, dz);
}

void InitializeCubeVerticies(R3Graph::R3Point cubeVertices[8],
    R3Graph::R3Point& cubeCenter, R3Graph::R3Point& bottomCenter,
    R3Graph::R3Point& topCenter, R3Graph::R3Point& leftCenter,
    R3Graph::R3Point& rightCenter, R3Graph::R3Point& frontCenter,
    R3Graph::R3Point& backCenter)
{
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
}

void InitializeExtendedVoxels(Voxel extendedVoxels[8], const Voxel& cube)
{
    for (int iv = 0; iv < 8; ++iv) {
        extendedVoxels[iv] = Voxel(
            cube.slice * 2, cube.point.x * 2, cube.point.y * 2
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
}

void FaceTriangulation(const VoxelSet& voxelSet, const Voxel& cube, const Voxel::Face& face,
    std::map<Voxel, int>& vertexIndices, Voxel extendedVoxels[8], Triangulation& triangulation,
    int indices[8], R3Graph::R3Point cubeVertices[8], std::map<int, std::set<int>>& VerxteNeighbours)
{
    std::vector<int> ind;
    ind.reserve(3);
    int i = 0, j = 0, k = 0, l = 0;
    InitializeVertexNumbers(face, i, j, k, l);

    if (!voxelSet.faceOpen(cube, face))
        return;
    else{
        // Add vertices
        AddVertex(vertexIndices, extendedVoxels, triangulation, cubeVertices, indices, i);
        AddVertex(vertexIndices, extendedVoxels, triangulation, cubeVertices, indices, j);
        AddVertex(vertexIndices, extendedVoxels, triangulation, cubeVertices, indices, k);
        AddVertex(vertexIndices, extendedVoxels, triangulation, cubeVertices, indices, l);

        // Add triangles for this face
        assert(
            indices[i] >= 0 &&
            indices[j] >= 0 &&
            indices[k] >= 0 &&
            indices[l] >= 0
        );

        // triangle verticies in clockwise direction
        triangulation.triangles.push_back(
            Triangulation::Triangle(
                indices[i], indices[j], indices[k]
            )
        );
        InitializeNormal(face, triangulation.triangles.back().Normal);      
        ind = { indices[i], indices[j], indices[k] };
        FillNeighbours(VerxteNeighbours, ind);

        triangulation.triangles.push_back(
            Triangulation::Triangle(
                indices[k], indices[l], indices[i]
            )
        );
        InitializeNormal(face, triangulation.triangles.back().Normal);
        ind = { indices[k], indices[l], indices[i] };
        FillNeighbours(VerxteNeighbours, ind);
    }
}

void AddVertex(std::map<Voxel, int>& vertexIndices, Voxel extendedVoxels[8], Triangulation& triangulation,
    R3Graph::R3Point cubeVertices[8], int indices[8], const int i)
{
    if (
        vertexIndices.count(
            extendedVoxels[i]
        ) == 0
        ) {
        // Add vertex to triangulation
        triangulation.vertices.push_back(
            cubeVertices[i]
        );
        indices[i] = (int)triangulation.vertices.size() - 1;
        vertexIndices[extendedVoxels[i]] = indices[i];
    }
    else {
        // Point is already in the array
        indices[i] = vertexIndices[extendedVoxels[i]];
    }
}

void InitializeVertexNumbers(const Voxel::Face& face, int& i, int& j, int& k, int& l)
{   // Look at cube illustration in computeTriangulationOfVoxelSet()
    // numeration in clockwise direction
    if (face == Voxel::Face::FACE_FRONT)
        i = 0, j = 1, k = 5, l = 4;
    
    else if (face == Voxel::Face::FACE_BACK)
        i = 2, j = 3, k = 7, l = 6;
    
    else if (face == Voxel::Face::FACE_LEFT)
        i = 0, j = 3, k = 7, l = 4;
    
    else if (face == Voxel::Face::FACE_RIGHT)
        i = 1, j = 2, k = 6, l = 5;
    
    else if (face == Voxel::Face::FACE_BOTTOM)
        i = 0, j = 1, k = 2, l = 3;
    
    else if (face == Voxel::Face::FACE_TOP)
        i = 4, j = 5, k = 6, l = 7;
    
}

void InitializeNormal(const Voxel::Face& face, R3Graph::R3Vector& Normal)
{
    if (face == Voxel::Face::FACE_FRONT)
        Normal = { 0., -1., 0. };
 
    else if (face == Voxel::Face::FACE_BACK)
        Normal = { 0., 1., 0. };
    
    else if (face == Voxel::Face::FACE_LEFT)
        Normal = { -1., 0., 0. };
    
    else if (face == Voxel::Face::FACE_RIGHT)
        Normal = { 1., 0., 0. };
    
    else if (face == Voxel::Face::FACE_BOTTOM)
        Normal = { 0., 0., -1. };
    
    else if (face == Voxel::Face::FACE_TOP)
        Normal = { 0., 0., 1. };
}

void computeTriangulationOfVoxelSet_MY(
    short* pointer,
    int threshold,
    Triangulation& triangulation,
    const VoxelSet& voxelSet,
    const R3Point& origin,
    double dx, double dy, double dz
) {
    //const VoxelBox& voxelBox = voxelSet.voxelBox;
    int sliceStart = 0, sliceFinish = 0, ixmin = 0, ixmax = 0, iymin = 0, iymax = 0;
    BoxBorders(voxelSet, sliceStart, sliceFinish, ixmin, ixmax, iymin, iymax);

    triangulation.clear();

    std::set<Voxel> Borderlayer, IsCalculate;
    FillBorderLayer(voxelSet, Borderlayer);

    for (const auto& cube : Borderlayer)
    {   
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

        R3Point cubeCenter, bottomCenter, topCenter, leftCenter, rightCenter, frontCenter, backCenter;
        InitializeNeighboursCentres(cube, origin, dx, dy, dz,
            cubeCenter, bottomCenter, topCenter, leftCenter, rightCenter, frontCenter, backCenter);

        R3Point cubeVertices[8];
        InitializeCubeVerticies(cubeVertices,
            cubeCenter, bottomCenter, topCenter, leftCenter, rightCenter, frontCenter, backCenter);

        std::vector<Voxel::Face> faces = { Voxel::FACE_FRONT,
        Voxel::FACE_BACK, Voxel::FACE_LEFT, Voxel::FACE_RIGHT, Voxel::FACE_BOTTOM, Voxel::FACE_TOP };

        for (const auto& face : faces)
        {
            const Voxel neighbour (cube.slice + Voxel::FACE_DIRECTIONS[face][0],
                cube.point.x + Voxel::FACE_DIRECTIONS[face][1],
                cube.point.y + Voxel::FACE_DIRECTIONS[face][2]);
            
            if(IsCalculate.count(neighbour) == 0)
                TriangulateBySkala(cube, face, Borderlayer, triangulation, origin, dx, dy, dz, pointer, cubeVertices, threshold);
        }
            

        IsCalculate.insert(cube);
    }
}

void TriangulateBySkala(const Voxel& cube, const Voxel::Face& face, const std::set<Voxel>& Borderlayer, Triangulation& triangulation,
    const R3Graph::R3Point& origin, double dx, double dy, double dz, short* pointer, R3Graph::R3Point cubeVertices[8], const short threshold)
{
    const Voxel Neighbour(cube.slice + Voxel::FACE_DIRECTIONS[face][0],
        cube.point.x + Voxel::FACE_DIRECTIONS[face][1],
        cube.point.y + Voxel::FACE_DIRECTIONS[face][2]);

    if (Borderlayer.count(Neighbour) == 0)
        return;
    else 
    {
        R3Point cubeCenter = voxel3DCoord(cube, origin, dx, dy, dz);
        R3Point NeighbourCenter = voxel3DCoord(Neighbour, origin, dx, dy, dz);

        R3Vector cubeGradient = cube.VoxelGradient(pointer, dx, dy, dz);// HU per mm
        double cubeDensity = (double)cube.VoxelDensity(pointer) - threshold;
        R3Vector NeighbourGradient = Neighbour.VoxelGradient(pointer, dx, dy, dz);// HU per mm
        double NeighbourDensity = (double)Neighbour.VoxelDensity(pointer) - threshold;

        R3Vector voxelDiagonalHalf = cubeCenter - cubeVertices[0];
        double ThresholdFunctionCubeCenter = cubeDensity + voxelDiagonalHalf.scalarProduct(cubeGradient);
        double ThresholdFunctionNeighbourCenter = NeighbourDensity + voxelDiagonalHalf.scalarProduct(NeighbourGradient);

        std::pair<R3Point, double> CubeCenterPair(cubeCenter, ThresholdFunctionCubeCenter);
        std::pair<R3Point, double> NeighbourCenterPair(NeighbourCenter, ThresholdFunctionNeighbourCenter);

        int i = 0, j = 0, k = 0, l = 0;
        InitializeVertexNumbers(face, i, j, k, l);
        const std::vector<int> indicies = { i, j, k, l };

        std::pair<R3Point, double> CubeVertexPairs[4];
        for (int i = 0; i < 4; ++i)
        {
            int ind = indicies[i];
            double ThresholdFunction = cubeDensity + (cubeVertices[ind] - cubeVertices[0]).scalarProduct(cubeGradient);
            CubeVertexPairs[i] = { cubeVertices[ind], ThresholdFunction };
        }

        for (int i = 0; i < 4; ++i) 
        {
            // cube edges: 0-1, 1-2, 2-3, (!) 3-0 (!)
            DensityTetrahedron CurrentTetrahedron(CubeCenterPair, NeighbourCenterPair,
                CubeVertexPairs[i], CubeVertexPairs[(i + 1) % 4]);
            triangulation.TriangulationOfTetrahedron(CurrentTetrahedron);
        }
    }
}

void FillBorderLayer(const VoxelSet& voxelSet, std::set<Voxel>& Borderlayer)
{
    const VoxelBox& voxelBox = voxelSet.voxelBox;
    int sliceStart = 0, sliceFinish = 0, ixmin = 0, ixmax = 0, iymin = 0, iymax = 0;
    BoxBorders(voxelSet, sliceStart, sliceFinish, ixmin, ixmax, iymin, iymax);
    
    std::vector<Voxel::Face> faces = {  Voxel::FACE_FRONT,
        Voxel::FACE_BACK, Voxel::FACE_LEFT, Voxel::FACE_RIGHT, Voxel::FACE_BOTTOM, Voxel::FACE_TOP };

    for (int slice = sliceFinish / 2; slice <= sliceFinish * 3 / 4; ++slice) 
    {
        for (int iy = iymax / 2; iy <= iymax * 3 / 4; ++iy) 
        {
            for (int ix = ixmax / 2; ix <= ixmax * 3 / 4; ++ix) 
            {
                if (voxelSet.voxelAt(slice, ix, iy) == 0)
                    continue;

                Voxel current{ slice, ix, iy };
                for (const auto& face : faces)
                    FillBorderNeighbour(current, face, Borderlayer, voxelSet);
            }
        }
    }
}

void FillBorderNeighbour(const Voxel& cube, const Voxel::Face& face, std::set<Voxel>& Borderlayer, const VoxelSet& voxelSet)
{
    if (!voxelSet.faceOpen(cube, face))
        return;
    else
    {
        Borderlayer.insert(cube);
        Voxel neighbour{
            cube.slice + Voxel::FACE_DIRECTIONS[face][0],
            cube.point.x + Voxel::FACE_DIRECTIONS[face][1],
            cube.point.y + Voxel::FACE_DIRECTIONS[face][2]
        };
        Borderlayer.insert(neighbour);

        // Add voxels around cube and its neighbour
        for (int i = -2; i < 3; ++i)
        {
            for (int j = -2; j < 3; ++j)
            {
                Voxel cubeSurrounder, neighbourSurrounder;
                if (face == Voxel::FACE_FRONT || face == Voxel::FACE_BACK)
                {
                    cubeSurrounder = Voxel{ cube.slice + i, cube.point.x + j, cube.point.y };
                    neighbourSurrounder = Voxel{ neighbour.slice + i, neighbour.point.x + j, neighbour.point.y };
                }
                else if (face == Voxel::FACE_LEFT || face == Voxel::FACE_RIGHT)
                {
                    cubeSurrounder = Voxel{  cube.slice + i, cube.point.x, cube.point.y + j };
                    neighbourSurrounder = Voxel{ neighbour.slice + i, neighbour.point.x, neighbour.point.y + j };
                }
                else if (face == Voxel::FACE_TOP || face == Voxel::FACE_BOTTOM)
                {
                    cubeSurrounder = Voxel{ cube.slice, cube.point.x + i, cube.point.y + j };
                    neighbourSurrounder = Voxel{ neighbour.slice, neighbour.point.x + i, neighbour.point.y + j };
                }
                Borderlayer.insert(cubeSurrounder);
                Borderlayer.insert(neighbourSurrounder);
            }
        }
    }
}

void InitializeVoxels(const Voxel& cube, Voxel& BottomVoxel, Voxel& TopVoxel, Voxel& LeftVoxel, Voxel& RightVoxel,
    Voxel& FrontVoxel, Voxel& BackVoxel)
{
    BottomVoxel = Voxel(cube.slice - 1, cube.point.x, cube.point.y);
    TopVoxel = Voxel(cube.slice + 1, cube.point.x, cube.point.y);
    LeftVoxel = Voxel(cube.slice, cube.point.x - 1, cube.point.y);
    RightVoxel = Voxel(cube.slice, cube.point.x + 1, cube.point.y);
    FrontVoxel = Voxel(cube.slice, cube.point.x, cube.point.y - 1);
    BackVoxel = Voxel(cube.slice, cube.point.x, cube.point.y + 1);
}

void SkalaTriangulation(const Voxel& cube, const Voxel::Face& face, short* pointer, double dx, double dy, double dz,
    int threshold, R3Graph::R3Point cubeVertices[8], std::pair<R3Graph::R3Point, double>& CubeCenterPair, Triangulation& triangulation,
    const VoxelSet& voxelSet, const R3Graph::R3Point& origin, const std::pair<R3Graph::R3Point, double> CubeVertexPairs[8])
{
    if (!voxelSet.faceOpen(cube, face))
        return;
    else
    {
        Voxel Neighbour(
            cube.slice + Voxel::FACE_DIRECTIONS[face][0],
            cube.point.x + Voxel::FACE_DIRECTIONS[face][1],
            cube.point.y + Voxel::FACE_DIRECTIONS[face][2]
        );
        R3Vector NeighbourVoxelGradient =
            Neighbour.VoxelGradient(pointer, dx, dy, dz); // HU per mm
        short FrontDensity = Neighbour.VoxelDensity(pointer);
        FrontDensity -= threshold;

        double ThresholdFunctionFrontCenter = FrontDensity +
            (CubeCenterPair.first - cubeVertices[0]).scalarProduct(NeighbourVoxelGradient);

        R3Point NeighbourCenter = voxel3DCoord(cube, origin, dx, dy, dz);
        std::pair<R3Point, double> NeighbourCenterPair(NeighbourCenter, ThresholdFunctionFrontCenter);

        int i = 0, j = 0, k = 0, l = 0;
        InitializeVertexNumbers(face, i, j, k, l);
        std::pair<R3Point, double> AdjacentCubeVertexPairs[4];
        AdjacentCubeVertexPairs[0] = CubeVertexPairs[i];
        AdjacentCubeVertexPairs[1] = CubeVertexPairs[j];
        AdjacentCubeVertexPairs[2] = CubeVertexPairs[k];
        AdjacentCubeVertexPairs[3] = CubeVertexPairs[l];

        for (int i = 0; i < 4; ++i) {
            // cube edges: 0-1, 1-2, 2-3, (!) 3-0 (!)
            DensityTetrahedron CurrentTetrahedron(CubeCenterPair, NeighbourCenterPair,
                AdjacentCubeVertexPairs[i], AdjacentCubeVertexPairs[(i + 1) % 4]);
            triangulation.TriangulationOfTetrahedron(CurrentTetrahedron);
        }
    }
}

void Triangulate_Custom(
    std::map<int, std::set<int>>& VerxteNeighbours,
    Triangulation& triangulation,
    const VoxelSet& voxelSet,
    const R3Point& origin,
    double dx, double dy, double dz
) {
    int sliceStart = 0, sliceFinish = 0, ixmin = 0, ixmax = 0, iymin = 0, iymax = 0;
    BoxBorders(voxelSet, sliceStart, sliceFinish, ixmin, ixmax, iymin, iymax);

    triangulation.clear();

    // Indices of extended voxels vertices in array
    // Each voxel produce 8 vertices == extended voxels
    // voxel (slice, x, y) -> 8 extended voxels:
    //       (2*slice + s0, 2*x + s1, 2*y + s2), where
    //                                           si = +-1
    std::map<Voxel, int> vertexIndices;

    //for (int slice = sliceFinish * 1 / 2; slice <= sliceFinish; ++slice)
    //{
    //    for (int iy = iymax * 1 / 3; iy <= iymax; ++iy) 
    //    {
    //        for (int ix = ixmax * 1 / 2; ix <= ixmax * 3 / 4; ++ix) 
    //        {
    for (int slice = sliceStart; slice <= sliceFinish; ++slice)
    {
        for (int iy = iymin; iy <= iymax; ++iy)
        {
            for (int ix = ixmin; ix <= ixmax; ++ix)
            {
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

                const Voxel cube(slice, ix, iy);

                R3Point cubeCenter, bottomCenter, topCenter, leftCenter, rightCenter, frontCenter, backCenter;
                InitializeNeighboursCentres(cube, origin, dx, dy, dz,
                    cubeCenter, bottomCenter, topCenter, leftCenter, rightCenter, frontCenter, backCenter);

                R3Point cubeVertices[8];
                InitializeCubeVerticies(cubeVertices,
                    cubeCenter, bottomCenter, topCenter, leftCenter, rightCenter, frontCenter, backCenter);

                Voxel extendedVoxels[8];
                InitializeExtendedVoxels(extendedVoxels, cube);

                int indices[9]; // Indices of vertices in array
                for (int iv = 0; iv < 9; ++iv)
                    indices[iv] = (-1);

                const std::vector<Voxel::Face> faces = { Voxel::FACE_FRONT,
                        Voxel::FACE_BACK, Voxel::FACE_LEFT, Voxel::FACE_RIGHT, Voxel::FACE_BOTTOM, Voxel::FACE_TOP };

                for (const auto& face : faces)
                {
                    int i = 0, j = 0, k = 0, l = 0;
                    InitializeVertexNumbers(face, i, j, k, l);

                    if (!voxelSet.faceOpen(cube, face))
                        continue;
                    else 
                    {
                        
                        // Add vertices
                        AddVertex(vertexIndices, extendedVoxels, triangulation, cubeVertices, indices, i);
                        AddVertex(vertexIndices, extendedVoxels, triangulation, cubeVertices, indices, j);
                        AddVertex(vertexIndices, extendedVoxels, triangulation, cubeVertices, indices, k);
                        AddVertex(vertexIndices, extendedVoxels, triangulation, cubeVertices, indices, l);

                        // Add triangles for this face
                        assert(
                            indices[i] >= 0 &&
                            indices[j] >= 0 &&
                            indices[k] >= 0 &&
                            indices[l] >= 0
                        );


                        Voxel extendedNeighbour(
                            2 * (cube.slice + Voxel::FACE_DIRECTIONS[face][0]),
                            2 * (cube.point.x + Voxel::FACE_DIRECTIONS[face][1]),
                            2 * (cube.point.y + Voxel::FACE_DIRECTIONS[face][2])
                        );
                        R3Point neighbourCenter;
                        if (face == Voxel::FACE_FRONT)
                        {
                            neighbourCenter = frontCenter;
                        }
                        else if(face == Voxel::FACE_BACK)
                        {
                            neighbourCenter = backCenter;
                        }
                        else if (face == Voxel::FACE_LEFT)
                        {
                            neighbourCenter = leftCenter;
                        }
                        else if (face == Voxel::FACE_RIGHT)
                        {
                            neighbourCenter = rightCenter;
                        }
                        else if (face == Voxel::FACE_BOTTOM)
                        {
                            neighbourCenter = bottomCenter;
                        }
                        else if (face == Voxel::FACE_TOP)
                        {
                            neighbourCenter = topCenter;
                        }

                        if (
                            vertexIndices.count(extendedNeighbour) == 0) {
                            // Add vertex to triangulation
                            triangulation.vertices.push_back(
                                neighbourCenter
                            );
                            indices[8] = (int)triangulation.vertices.size() - 1;
                            vertexIndices[extendedNeighbour] = indices[8];
                        }
                        else {
                            // Point is already in the array
                            indices[8] = vertexIndices[extendedNeighbour];
                        }

                        // 1
                        triangulation.triangles.push_back(
                            Triangulation::Triangle(
                                indices[i], indices[j], indices[8]
                            )
                        );
                        InitializeNormal_Custom(face, triangulation.triangles.back().Normal);
                        {
                            std::vector<int> ind = { indices[i], indices[j], indices[8] };
                            FillNeighbours(VerxteNeighbours, ind);
                        }
                        // 2
                        triangulation.triangles.push_back(
                            Triangulation::Triangle(
                                indices[j], indices[k], indices[8]
                            )
                        );
                        InitializeNormal_Custom(face, triangulation.triangles.back().Normal);
                        {
                            std::vector<int> ind = { indices[j], indices[k], indices[8] };
                            FillNeighbours(VerxteNeighbours, ind);
                        }

                        // 3
                        triangulation.triangles.push_back(
                            Triangulation::Triangle(
                                indices[k], indices[l], indices[8]
                            )
                        );
                        InitializeNormal_Custom(face, triangulation.triangles.back().Normal);
                        {
                            std::vector<int> ind = { indices[k], indices[l], indices[8] };
                            FillNeighbours(VerxteNeighbours, ind);
                        }

                        // 4 
                        triangulation.triangles.push_back(
                            Triangulation::Triangle(
                                indices[l], indices[i], indices[8]
                            )
                        );
                        InitializeNormal_Custom(face, triangulation.triangles.back().Normal);
                        {
                            std::vector<int> ind = { indices[l], indices[i], indices[8] };
                            FillNeighbours(VerxteNeighbours, ind);
                        }
                    }
                }
            }
        }
    }
}

void InitializeNormal_Custom(const Voxel::Face& face, R3Graph::R3Vector& Normal)
{//TODO Repair Normals
    if (face == Voxel::Face::FACE_FRONT)
        Normal = { 0., -1., 0. };

    else if (face == Voxel::Face::FACE_BACK)
        Normal = { 0., 1., 0. };

    else if (face == Voxel::Face::FACE_LEFT)
        Normal = { -1., 0., 0. };

    else if (face == Voxel::Face::FACE_RIGHT)
        Normal = { 1., 0., 0. };

    else if (face == Voxel::Face::FACE_BOTTOM)
        Normal = { 0., 0., -1. };

    else if (face == Voxel::Face::FACE_TOP)
        Normal = { 0., 0., 1. };
}
