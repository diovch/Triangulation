#ifndef VOXELSET_H
#define VOXELSET_H


//#include <vector>
//#include <QString>
#include <string>
#include <set>
#include "r2geom.h"
#include "roi.h"        /* For Bitmask */
#include "R3Graph.h"

class Triangulation;

class Voxel {
public:
    int slice;
    I2Point point;

    Voxel() :
        slice(0),
        point()
    {}

    Voxel(int s, const I2Point& p) :
        slice(s),
        point(p)
    {}

    Voxel(int s, int x, int y) :
        slice(s),
        point(x, y)
    {}

    bool operator==(const Voxel& v) const {
        return (
            slice == v.slice && point == v.point
            );
    }

    bool operator!=(const Voxel& v) const {
        return !(*this == v);
    }

    bool operator<(const Voxel& v) const {
        return (
            slice < v.slice || (
                slice <= v.slice && point < v.point
                )
            );
    }

    bool operator<=(const Voxel& v) const {
        return (
            slice < v.slice || (
                slice <= v.slice && point <= v.point
                )
            );
    }

    bool operator>(const Voxel& v) const {
        return !(*this <= v);
    }

    bool operator>=(const Voxel& v) const {
        return !(*this < v);
    }

    Voxel operator+(const Voxel& v) const {
        return Voxel(
            slice + v.slice, point + v.point
        );
    }

    // Enumeration of voxel faces
    enum Face{
        FACE_LEFT = 0,
        FACE_RIGHT = 1,
        FACE_FRONT = 2,
        FACE_BACK = 3,
        FACE_BOTTOM = 4,
        FACE_TOP = 5
    };

    // For given face, contains shifts in slice, x, y = -1, 0, 1
    static const int FACE_DIRECTIONS[6][3]; // Defined in "VoxelSet.cpp"
    short VoxelDensity(short* p) {
        auto VoxelPointer = p + (point.x + point.y * 512 + slice * 512 * 512);
        return *VoxelPointer;
    }

    R3Graph::R3Vector VoxelGradient(short* pointer, double dx, double dy, double dz){// HU per mm
    
        Voxel RightVoxel = Voxel(slice, point.x + 1, point.y);
        Voxel BackVoxel = Voxel(slice, point.x, point.y + 1);
        Voxel TopVoxel = Voxel(slice + 1, point.x, point.y);
        double density = (double)(*this).VoxelDensity(pointer);

        R3Graph::R3Vector result = R3Graph::R3Vector(
            ((double)(RightVoxel.VoxelDensity(pointer)) - density) / dx,
            ((double)(BackVoxel.VoxelDensity(pointer)) - density) / dy,
            ((double)(TopVoxel.VoxelDensity(pointer)) - density) / dz
        );
         return result;
    }
};

class VoxelBox {
public:
    Voxel origin;
    int width;
    int depth;
    int height;

public:
    VoxelBox() :
        origin(),
        width(0),
        depth(0),
        height(0)
    {}

    VoxelBox(const Voxel& v, int w, int d, int h) :
        origin(v),
        width(w),
        depth(d),
        height(h)
    {}

    int xMin() const { return origin.point.x; }
    int yMin() const { return origin.point.y; }
    int sliceMin() const { return origin.slice; }
    int xMax() const { return origin.point.x + width; }
    int yMax() const { return origin.point.y + depth; }
    int sliceMax() const { return origin.slice + height; }
    void setXMin(int x_min) { origin.point.x = x_min; }
    void setXMax(int x_max) { width = x_max - origin.point.x; }
    void setYMin(int y_min) { origin.point.y = y_min; }
    void setYMax(int y_max) { depth = y_max - origin.point.y; }
    void setSliceMin(int s) { origin.slice = s; }
    void setSliceMax(int s) { height = s - origin.slice; }
};

class VoxelSet {
public:
    int xMax;
    int yMax;
    int maxSlices;
    std::vector<Bitmask> bitmasks;
    double numVoxels;
    VoxelBox voxelBox;
    int threshold;
    bool detected3D;    // 3D detection is done

    VoxelSet() :
        xMax(512),
        yMax(512),
        maxSlices(0),
        bitmasks(),
        numVoxels(0.),
        voxelBox(),
        threshold(300),
        detected3D(false)
    {}

    void initialize(
        int x_max, int y_max, int max_slices
    );

    void clear();

    size_t size() const { return size_t(numVoxels); }

    const Bitmask& bitmaskAt(int slice) const {
        return bitmasks.at(slice);
    }
    Bitmask& bitmaskAt(int slice) {
        return bitmasks.at(slice);
    }

    int voxelAt(int slice, int x, int y) const {
        return bitmasks.at(slice).pixelAt(x, y);
    }
    int voxelAt(const Voxel& v) const {
        return voxelAt(v.slice, v.point.x, v.point.y);
    }

    void setVoxelValue(int slice, int x, int y, int v) {
        bitmasks.at(slice).setPixValue(x, y, v);
    }
    void setVoxelValue(const Voxel& v, int val) {
        setVoxelValue(v.slice, v.point.x, v.point.y, val);
    }

    void addVoxel(const Voxel& v, int val) {
        setVoxelValue(v, val);
        if (numVoxels == 0) {
            // Initialize a voxel box
            voxelBox.origin = v;
            voxelBox.width = 0;
            voxelBox.depth = 0;
            voxelBox.height = 0;
        }
        else {
            if (v.point.x < voxelBox.xMin()) {
                voxelBox.width += voxelBox.xMin() - v.point.x;
                voxelBox.setXMin(v.point.x);
            }
            else if (v.point.x > voxelBox.xMax()) {
                voxelBox.width += v.point.x - voxelBox.xMax();
            }

            if (v.point.y < voxelBox.yMin()) {
                voxelBox.depth += voxelBox.yMin() - v.point.y;
                voxelBox.setYMin(v.point.y);
            }
            else if (v.point.y > voxelBox.yMax()) {
                voxelBox.depth += v.point.y - voxelBox.yMax();
            }

            if (v.slice < voxelBox.sliceMin()) {
                voxelBox.height += voxelBox.sliceMin() - v.slice;
                voxelBox.setSliceMin(v.slice);
            }
            else if (v.slice > voxelBox.sliceMax()) {
                voxelBox.height += v.slice - voxelBox.sliceMax();
            }
        }
        numVoxels += 1.;
    }

    bool voxelInVolume(const Voxel& v) const {
        return (
            0 <= v.slice && v.slice < maxSlices &&
            0 <= v.point.x && v.point.x < xMax &&
            0 <= v.point.y && v.point.y < yMax
            );
    }

    bool faceOpen(const Voxel& v, int face) const; // Voxel::FACE_LEFT, etc.

    void computeVoxelBox();

    //bool save(QString path) const;
    //bool load(QString path);
    bool save(std::string path) const;
    bool load(std::string path);


};

class PackedVoxelSet {
public:
    int xMax;
    int yMax;
    int maxSlices;
    std::vector<PackedBitmask> packedBitmasks;
    double numVoxels;
    VoxelBox voxelBox;
    int threshold;
    bool detected3D;    // 3D detection is done

    PackedVoxelSet() :
        xMax(512),
        yMax(512),
        maxSlices(0),
        packedBitmasks(),
        numVoxels(0.),
        voxelBox(),
        threshold(300),
        detected3D(false)
    {}

    PackedVoxelSet(const VoxelSet& vs);
    PackedVoxelSet& operator=(const VoxelSet& vs);
    VoxelSet& unpack(VoxelSet& vs) const;
};

double detectVoxelSet(
    double (*f)(const Voxel&),
    double threshold,
    const VoxelBox& voxelBox,
    const Voxel& seed,
    VoxelSet& voxelSet
);

double detectVoxelSetFromCta(
    double (*f)(const Voxel&, short* pointer),
    double threshold,
    const VoxelBox& voxelBox,
    const Voxel& seed,
    short* pointer,
    VoxelSet& voxelSet
);

void computeTriangulationOfVoxelSet(
    std::map<int, std::set<int>>& neighbours,
    Triangulation& triangulation,
    const VoxelSet& voxelSet,
    const R3Graph::R3Point& origin,
    double dx, double dy, double dz
);

void FillNeighbours(std::map<int, std::set<int>>&, std::vector<int>&);

void computeTriangulationOfVoxelSet_MY(
    short* pointer,
    int threshold,
    Triangulation& triangulation,
    const VoxelSet& voxelSet,
    const R3Graph::R3Point& origin,
    double dx, double dy, double dz
);

inline R3Graph::R3Point voxel3DCoord(
    const Voxel& voxel,
    const R3Graph::R3Point& origin,
    double dx, double dy, double dz
) {
    return R3Graph::R3Point(
        origin.x + voxel.point.x * dx + dx / 2.,
        origin.y + voxel.point.y * dy + dy / 2.,
        origin.z + voxel.slice * dz + dz / 2.
    );
}

void BoxBorders(const VoxelBox&, const VoxelSet&, int&, int&, int&, int&, int&, int&);

void InitializeNeighboursCentres(const Voxel& cube, 
    const R3Graph::R3Point& origin, double dx, double dy,
    double dz, R3Graph::R3Point& cubeCenter,
    R3Graph::R3Point& bottomCenter, R3Graph::R3Point& topCenter,
    R3Graph::R3Point& leftCenter, R3Graph::R3Point& rightCenter,
    R3Graph::R3Point& frontCenter, R3Graph::R3Point& backCenter);

void InitializeCubeVerticies(R3Graph::R3Point cubeVertices[8],
    R3Graph::R3Point& cubeCenter, R3Graph::R3Point& bottomCenter, 
    R3Graph::R3Point& topCenter, R3Graph::R3Point& leftCenter, 
    R3Graph::R3Point& rightCenter, R3Graph::R3Point& frontCenter, 
    R3Graph::R3Point& backCenter);

void InitializeExtendedVoxels(Voxel extendedVoxels[8], const Voxel& cube);

void FaceTriangulation(const VoxelSet& voxelSet, const Voxel& cube, const Voxel::Face& face,
    std::map<Voxel, int>& vertexIndices, Voxel extendedVoxels[8], Triangulation& triangulation,
    int indices[8], R3Graph::R3Point cubeVertices[8], std::map<int, std::set<int>>& VerxteNeighbours);

void InitializeVertexNumbers(const Voxel::Face& face, int& i, int& j, int& k, int& l);

void InitializeNormal(const Voxel::Face& face, R3Graph::R3Vector& Normal);

void AddVertex(std::map<Voxel, int>& vertexIndices, Voxel extendedVoxels[8], Triangulation& triangulation, 
    R3Graph::R3Point cubeVertices[8], int indices[8], const int i);

void InitializeVoxels(Voxel& cube, Voxel& BottomVoxel, Voxel& TopVoxel, Voxel& LeftVoxel, Voxel& RightVoxel,
    Voxel& FrontVoxel, Voxel& BackVoxel);

#endif
