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
    double threshold,
    const VoxelBox& voxelBox,
    const Voxel& seed,
    short* pointer,
    unsigned char* mask_pointer,
    unsigned char maskLabel,
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

    if (VoxelDensity(seed, pointer) < threshold || VoxelType(seed, mask_pointer) != maskLabel)
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

            
            if (VoxelDensity(n, pointer) < threshold || VoxelType(n, mask_pointer) != maskLabel)
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

void FillVoids(VoxelSet& voxelSet)
{// TODO Debug FillVoids()
    int sliceStart = 0, sliceFinish = 0, ixmin = 0, ixmax = 0, iymin = 0, iymax = 0;
    BoxBorders(voxelSet, sliceStart, sliceFinish, ixmin, ixmax, iymin, iymax);

    VoxelBox& voxelBox = voxelSet.voxelBox;
    int NumVoxels = voxelBox.depth * voxelBox.width * voxelBox.height;
    int NumOutsideVoids = 0;

    std::set<Voxel> OutsideVoids;
    std::deque<Voxel> QueueOfNeighbours;

    for (int slice = sliceStart; slice <= sliceFinish; ++slice)
    {
        Voxel seed = SearchSliceSeed(voxelSet, slice);
        QueueOfNeighbours.push_back(seed);
        OutsideVoids.insert(seed);

        while (!QueueOfNeighbours.empty())
        {
            Voxel t = QueueOfNeighbours.front(); QueueOfNeighbours.pop_front();
            for (int i = 0; i < 4; ++i) {
                Voxel n = t + NEIGHBOURS6[i]; // at slice
                if (OutsideVoids.count(n) == 1 ||
                    n.slice < sliceStart || n.slice > sliceFinish ||
                    n.point.x < ixmin || n.point.x > ixmax ||
                    n.point.y < iymin || n.point.y > iymax
                    )
                    continue;

                if (voxelSet.voxelAt(n) == 0)
                {
                    QueueOfNeighbours.push_back(n);
                    OutsideVoids.insert(n);
                    NumOutsideVoids++;
                }
            }
        }

        for (int iy = iymin; iy <= iymax; ++iy)
        {
            for (int ix = ixmin; ix <= ixmax; ++ix)
            {
                Voxel current(slice, ix, iy);
                if (voxelSet.voxelAt(current) == 0 && OutsideVoids.count(current) == 0)
                {
                    voxelSet.addVoxel(current, ROI_POSITIVE);
                }
            }
        }
        OutsideVoids.clear();
    }

}

Voxel SearchSeed(VoxelSet& voxelSet)
{
    int sliceStart = 0, sliceFinish = 0, ixmin = 0, ixmax = 0, iymin = 0, iymax = 0;
    BoxBorders(voxelSet, sliceStart, sliceFinish, ixmin, ixmax, iymin, iymax);

    for (int slice = sliceStart; slice <= sliceFinish; ++slice)
    {
        for (int iy = iymin; iy <= iymax; ++iy)
        {
            for (int ix = ixmin; ix <= ixmax; ++ix)
            {
                Voxel seed(slice, ix, iy);
                if (voxelSet.voxelAt(seed) == 0)
                {
                    return seed;
                }
            }
        }
    }
    return {};
}

int CountRoiVoxels(VoxelSet& voxelSet)
{
    int num = 0;

    int sliceStart = 0, sliceFinish = 0, ixmin = 0, ixmax = 0, iymin = 0, iymax = 0;
    BoxBorders(voxelSet, sliceStart, sliceFinish, ixmin, ixmax, iymin, iymax);

    for (int slice = sliceStart; slice <= sliceFinish; ++slice)
        for (int iy = iymin; iy <= iymax; ++iy)
            for (int ix = ixmin; ix <= ixmax; ++ix)
            {
                Voxel seed(slice, ix, iy);
                if (voxelSet.voxelAt(seed) != 0)
                    ++num;
            }
    return num;
}

Voxel SearchSliceSeed(const VoxelSet & voxelSet, int slice)
{
    int sliceStart = 0, sliceFinish = 0, ixmin = 0, ixmax = 0, iymin = 0, iymax = 0;
    BoxBorders(voxelSet, sliceStart, sliceFinish, ixmin, ixmax, iymin, iymax);

    
    for (int iy = iymin; iy <= iymax; ++iy)
    {
        for (int ix = ixmin; ix <= ixmax; ++ix)
        {
            Voxel seed(slice, ix, iy);
            if (voxelSet.voxelAt(seed) == 0)
            {
                return seed;
            }
        }
    }
    return {};
}


void ContourSegmentation(VoxelSet& voxelSet, const Voxel& seed)
{   // Moore NeighborhoodTracing
    if (voxelSet.voxelAt(seed) == 0)
        return;

    int zStart = seed.slice, xStart = seed.point.x, yStart = seed.point.y;
    Voxel bug = seed;

    for (int slice = 355; slice < voxelSet.maxSlices; ++slice)
    {
        DirectionOfMovement direction = Y_NEGATIVE;
        while (voxelSet.voxelAt(bug) != ROI_POSITIVE) // first free voxel
        {
            bug.point.y--;
        }
        voxelSet.addVoxel(bug, ROI_MANUAL_NEGATIVE_BORDER);
        I2Point initialPosition{ bug.point }, currentPosition{};
        Voxel leftHandNeighbour, forwardNeighbour;


        //    _     bug has an origin which it comes from and three futher directions - left, forward, right
        //  _| |_   if left spot is free it turns left
        // |_   _|  if forward - it moves forward
        //   |_|    otherwise it turns right, so ROI area is being on the left hand 


        while (currentPosition != initialPosition)
        {
            leftHandNeighbour = LeftHandNeighbour(bug, direction);
            forwardNeighbour = bug + NEIGHBOURS6[direction];

            if (voxelSet.voxelAt(leftHandNeighbour) != ROI_POSITIVE) 
            {
                bug = leftHandNeighbour;
                TurnLeft(direction);
            }
            else if (voxelSet.voxelAt(forwardNeighbour) != ROI_POSITIVE)
            {
                bug = forwardNeighbour;
            }
            else
            {
                bug = RightHandNeighbour(bug, direction);
            }
            voxelSet.addVoxel(bug, ROI_MANUAL_NEGATIVE_BORDER);
            currentPosition = bug.point;
        }
        
    }
}

Voxel LeftHandNeighbour(const Voxel& bug, const DirectionOfMovement& direction)
{
    if (direction == Y_NEGATIVE)
    {
        return { bug.slice, bug.point.x + 1, bug.point.y };
    }
    else if (direction == Y_POSITIVE)
    {
        return { bug.slice, bug.point.x - 1, bug.point.y };
    }
    else if (direction == X_NEGATIVE)
    {
        return { bug.slice, bug.point.x, bug.point.y - 1};
    }
    else if(direction == X_POSITIVE)
    {
        return { bug.slice, bug.point.x, bug.point.y + 1};
    }
    else
    {
        return {};
    }
}

Voxel RightHandNeighbour(const Voxel& bug, const DirectionOfMovement& direction)
{
    if (direction == Y_NEGATIVE)
    {
        return { bug.slice, bug.point.x - 1, bug.point.y };
    }
    else if (direction == Y_POSITIVE)
    {
        return { bug.slice, bug.point.x + 1, bug.point.y };
    }
    else if (direction == X_NEGATIVE)
    {
        return { bug.slice, bug.point.x, bug.point.y + 1 };
    }
    else if (direction == X_POSITIVE)
    {
        return { bug.slice, bug.point.x, bug.point.y - 1 };
    }
    else
    {
        return {};
    }
}

void TurnLeft(DirectionOfMovement& direction)
{
    if (direction == Y_NEGATIVE)
    {
        direction = X_POSITIVE;
    }
    else if (direction == Y_POSITIVE)
    {
        direction = X_NEGATIVE;
    }
    else if (direction == X_NEGATIVE)
    {
        direction = Y_NEGATIVE;
    }
    else if (direction == X_POSITIVE)
    {
        direction = Y_POSITIVE;
    }
}

unsigned char VoxelType(const Voxel& v, unsigned char* p)
{
    auto VoxelPointer = p + (v.point.x + v.point.y * 512 + v.slice * 512 * 512);
    return *VoxelPointer;
}


short VoxelDensity(const Voxel& v, short* p)
{
    auto VoxelPointer = p + (v.point.x + v.point.y * 512 + v.slice * 512 * 512);
    return *VoxelPointer;
}