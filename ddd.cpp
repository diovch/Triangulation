#include <deque>
#include <utility>
#include <climits>
#include <cassert>
#include "voxelset.h"
#include "R3Graph.h"

using namespace R3Graph;

static const Voxel NEIGHBOURS6[6] = {
    Voxel(0, 1, 0),
    Voxel(0, -1, 0),
    Voxel(0, 0, 1),
    Voxel(0, 0, -1),
    Voxel(1, 0, 0),
    Voxel(-1, 0, 0)
};

double detectVoxelSet(
    double (*f)(const Voxel&),
    double threshold,
    const VoxelBox& voxelBox,
    const Voxel& seed,
    VoxelSet& voxelSet
) {
    int sliceStart = voxelBox.origin.slice;
    int sliceFinish = sliceStart + voxelBox.height - 1;
    int xmax = voxelBox.width, ymax = voxelBox.depth;
    int zmax = voxelBox.height;

    voxelSet.initialize(xmax, ymax, zmax);

    if (seed.slice < sliceStart || seed.slice > sliceFinish)
        return 0.;

    if ((*f)(seed) < threshold)
        return 0.;

    std::deque<Voxel> deq;
    voxelSet.setVoxelValue(seed, ROI_POSITIVE);
    deq.push_back(seed);
    double num = 1.;
    double maxNum = double(sliceFinish + 1 - sliceStart) *
        double(xmax) * double(ymax);
    while (!deq.empty()) {
        Voxel t = deq.front(); deq.pop_front();
        for (int i = 0; i < 6; ++i) {
            Voxel n = t + NEIGHBOURS6[i];
            if (
                n.slice < sliceStart || n.slice > sliceFinish ||
                n.point.x <= 0 || n.point.x >= xmax ||
                n.point.y <= 0 || n.point.y >= ymax
                )
                continue;
            if (voxelSet.voxelAt(n) != 0)
                continue;   // Voxel is already in the set

            double v = (*f)(n);
            if (v < threshold)
                continue;
            voxelSet.addVoxel(n, ROI_POSITIVE);
            deq.push_back(n);
            num += 1.;
        }
        if (num >= maxNum) {
            // Should not come here!
            break;
        }
    }
    voxelSet.detected3D = true;
    return num;
}

// my realization with image pointer
double detectVoxelSetFromCta(
    double (*f)(const Voxel&, short* pointer),
    double threshold,
    const VoxelBox& voxelBox,
    const Voxel& seed,
    short* pointer,
    VoxelSet& voxelSet
) {
    int sliceStart = voxelBox.origin.slice;
    int sliceFinish = sliceStart + voxelBox.height - 1;
    int xmax = voxelBox.width, ymax = voxelBox.depth;
    int zmax = voxelBox.height;
    //                   512   512   389
    voxelSet.initialize(xmax, ymax, zmax);

    if (seed.slice < sliceStart || seed.slice > sliceFinish)
        return 0.;

    if ((*f)(seed, pointer) < threshold)
        return 0.;

    std::deque<Voxel> deq;
    voxelSet.setVoxelValue(seed, ROI_POSITIVE); // Write 1 in bitmask
    deq.push_back(seed);
    double num = 1.;
    double maxNum = double(sliceFinish + 1 - sliceStart) *
        double(xmax) * double(ymax);
    while (!deq.empty()) {
        Voxel t = deq.front(); deq.pop_front();
        for (int i = 0; i < 6; ++i) {
            Voxel n = t + NEIGHBOURS6[i];
            if (
                n.slice < sliceStart || n.slice > sliceFinish ||
                n.point.x <= 0 || n.point.x >= xmax ||
                n.point.y <= 0 || n.point.y >= ymax
                )
                continue;
            if (voxelSet.voxelAt(n) != 0)
                continue;   // Voxel is already in the set

            double v = (*f)(n, pointer);
            if (v < threshold)
                continue;
            voxelSet.addVoxel(n, ROI_POSITIVE);
            deq.push_back(n);
            num += 1.;
        }
        if (num >= maxNum) {
            // Should not come here!
            break;
        }
    }
    voxelSet.detected3D = true;
    return num;
}