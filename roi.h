#ifndef ROI_H
#define ROI_H

#include <vector>
#include <cassert>
//#include <QDataStream>
#include <fstream>
#include "matrix.h"
#include "r2geom.h"

class QDataStream;

const int ROI_POSITIVE_BIT = 1; // Positive mask
const int ROI_NEGATIVE_BIT = 2; // Negative mask

const int ROI_POSITIVE = 1;
const int ROI_NEGATIVE = 2;
const int ROI_MANUAL_BORDER = 5;                // 0101
const int ROI_MANUAL_NEGATIVE_BORDER = 6;       // 0110
const int ROI_CONNECTING_VALUE = 9;             // 1001
const int ROI_SPECIAL_BORDER_VALUE = 255;

class Bitmask : public ByteMatrix {
private:
    void resize(int w, int h) {
        ByteMatrix::resize(h, w);
    }

public:
    Bitmask(int w = 1, int h = 1) :
        ByteMatrix(h, w)
    {}

    int width() const { return cols(); }
    int height() const { return rows(); }
    int getXMax() const { return width(); }
    int getYMax() const { return height(); }

    int pixelAt(int x, int y) const {
        return at(y, x);
    }
    unsigned char& pixelRef(int x, int y) {
        return at(y, x);
    }
    int pixelAt(const I2Point& p) const {
        return pixelAt(p.x, p.y);
    }
    unsigned char& pixelRef(const I2Point& p) {
        return pixelRef(p.x, p.y);
    }
    void setPixValue(int x, int y, int v) {
        pixelRef(x, y) = (unsigned char)(v);
    }
    int getPixValue(int x, int y) const {
        return pixelAt(x, y);
    }
    void setPixValue(const I2Point& p, int v) {
        pixelRef(p) = (unsigned char)(v);
    }
    int getPixValue(const I2Point& p) const {
        return pixelAt(p);
    }

    void resizeBitmask(int w, int h) {
        ByteMatrix::resize(h, w);
    }

    bool isPositive(int x, int y) const {
        return ((pixelAt(x, y) & ROI_POSITIVE) != 0);
    }
    void setPositive(int x, int y) {
        pixelRef(x, y) |= (unsigned char)(ROI_POSITIVE);
    }
    void clearPositive(int x, int y) {
        pixelRef(x, y) &= ~((unsigned char)(ROI_POSITIVE));
    }

    bool isNegative(int x, int y) const {
        return ((pixelAt(x, y) & ROI_NEGATIVE) != 0);
    }
    void setNegative(int x, int y) {
        pixelRef(x, y) |= (unsigned char)(ROI_NEGATIVE);
    }
    void clearNegative(int x, int y) {
        pixelRef(x, y) &= ~((unsigned char)(ROI_NEGATIVE));
    }

    void clear() {
        zero();
    }

    void fill(int v = 0) {
        ByteMatrix::fill((unsigned char)v);
    }

    void drawLine(
        const I2Point& p0, const I2Point& p1, int v = ROI_POSITIVE,
        bool thickLine = false
    );

    void drawContour(
        const I2Contour& c, int v = ROI_POSITIVE,
        bool closed = true, bool thickLine = false
    );
    void drawContour(
        const R2Contour& c, int v = ROI_POSITIVE
    );

    void paintContour(
        const I2Contour& c,
        int v = ROI_POSITIVE,
        int vBorder = ROI_POSITIVE,
        int numNeighbours = 4,
        const I2Point* seed = 0,
        const Bitmask* externalBitmask = 0
    );
    void paintContour(
        const R2Contour& c, int v = ROI_POSITIVE
    );
    void regionGrow4(
        const I2Point& p, int v = ROI_POSITIVE,
        const Bitmask* externalBitmask = 0
    );
    void regionGrow8(
        const I2Point& p, int v = ROI_POSITIVE,
        const Bitmask* externalBitmask = 0
    );

    void paintContourInternalArea(
        const I2Contour& c,
        int v = ROI_POSITIVE,
        int numNeighbours = 4,
        const I2Point* seed = 0,
        const Bitmask* externalBitmask = 0
    );

    bool write(std::fstream& s) const;
    bool read(std::fstream& s);
};

class PackedBitmask {
public:
    int width;
    int height;
    std::vector<unsigned char> packedMatrix;

public:
    PackedBitmask() :
        width(0),
        height(0),
        packedMatrix()
    {}
    PackedBitmask(const Bitmask& m);

    PackedBitmask& operator=(const Bitmask& m);
    Bitmask& unpack(Bitmask& m) const;
};

typedef std::vector<I2Contour> I2ContourVector;

class ROI {
public:
    int sliceMin;
    int sliceMax;
    int maxSlices;
    int xMax;
    int yMax;
    std::vector<Bitmask> bitmasks;
    std::vector<I2ContourVector> contours;
    I2Rectangle framingRect;
    int firstVertBorder;    // (-1) for undefined
    int secondVertBorder;   // (-1) for undefined

public:
    ROI() :
        sliceMin(0),
        sliceMax(0),
        maxSlices(0),
        xMax(512),
        yMax(512),
        bitmasks(),
        contours(),
        framingRect(0, 0, 512, 512),
        firstVertBorder(-1),
        secondVertBorder(-1)
    {}

    int lowerSlice() const;
    int upperSlice() const;

    void initialize(
        int max_slices,
        int w, int h
    );
    void clear();

    int pixelAt(int sliceIdx, int x, int y) const {
        return bitmasks.at(sliceIdx).pixelAt(x, y);
    }
    unsigned char& pixelRef(int sliceIdx, int x, int y) {
        return bitmasks.at(sliceIdx).pixelRef(x, y);
    }
    const Bitmask& bitmaskAt(int sliceIdx) const {
        return bitmasks.at(sliceIdx);
    }
    Bitmask& bitmaskAt(int sliceIdx) {
        return bitmasks.at(sliceIdx);
    }

    I2ContourVector& contourArrayAt(int sliceIdx) {
        return contours.at(sliceIdx);
    }
    const I2ContourVector& contourArrayAt(int sliceIdx) const {
        return contours.at(sliceIdx);
    }

    void updateROI(int sliceIdx);
};

class RLEEncoder {
private:
    //QDataStream* stream;
    std::fstream* stream;
    std::vector<unsigned char>* packedArray;
    unsigned char currentChar;
    int counter;
    int currentPos; // position in packedArray (used in reading)

public:
    RLEEncoder(
        //QDataStream* s = 0,
        std::fstream* s = 0,
        std::vector<unsigned char>* packed_array = 0
    ) :
        stream(s),
        packedArray(packed_array),
        currentChar(0),
        counter(0),
        currentPos(0)
    {}

    void initialize(
        //QDataStream* s = 0,
        std::fstream* s = 0,
        std::vector<unsigned char>* packed_array = 0
    ) {
        stream = s;
        counter = 0;
        currentChar = 0;
        packedArray = packed_array;
        if (packedArray != 0)
            packedArray->clear();
        currentPos = 0;
    }

    void writeByte(unsigned char b) {
        if (counter == 0) {
            currentChar = b;
            counter = 1;
        }
        else if (b == currentChar) {
            ++counter;
            if (counter == 255)
                flush();
        }
        else {
            flush();
            currentChar = b;
            counter = 1;
        }
    }

    void push_back(unsigned char b) {
        assert(packedArray != 0);
        if (counter == 0) {
            currentChar = b;
            counter = 1;
        }
        else if (b == currentChar) {
            ++counter;
            if (counter == 255)
                flush_array();
        }
        else {
            flush_array();
            currentChar = b;
            counter = 1;
        }
    }

    void flush() {
        assert(stream != 0);
        if (stream != 0 && counter != 0) {
            assert(counter <= 255);
            (*stream) << ((unsigned char)counter);
            (*stream) << currentChar;
        }
        counter = 0;
    }

    void flush_array() {
        assert(packedArray != 0);
        if (packedArray != 0 && counter != 0) {
            assert(counter <= 255);
            packedArray->push_back((unsigned char)counter);
            packedArray->push_back(currentChar);
            currentPos += 2;
        }
        counter = 0;
    }

    bool atEnd() const {
        if (stream != 0) {
            //return (counter == 0 && stream->atEnd());
            return (counter == 0 && stream->eof());
        }
        else {
            assert(packedArray != 0);
            return (
                packedArray == 0 ||
                (counter == 0 &&
                    currentPos >= (int)packedArray->size())
                );
        }
    }

    unsigned char readByte() {  // Reading from a file
        if (counter == 0) {
            readRunningSequence();
            if (counter == 0 && atEnd()) {
                return 0;       // Reading past end of file!
            }
        }
        // assert(counter > 0);
        if (counter == 0) {
            // Probably, incorrect format of input stream
            return 0;
        }
        --counter;
        return currentChar;
    }

    unsigned char pop_front() { // Reading from an array
        if (counter == 0) {
            extractRunningSequence();
            if (counter == 0 && atEnd()) {
                return 0;       // Reading past end of array!
            }
        }
        // assert(counter > 0);
        if (counter == 0) {
            // Probably, incorrect format of input stream
            return 0;
        }
        --counter;
        return currentChar;
    }


private:
    void readRunningSequence() { // Reading from a file
        assert(counter == 0);
        //if (stream == 0 || stream->atEnd())
        if (stream == 0 || stream->eof())
            return;
        unsigned char b;
        (*stream) >> b;
        counter = b;
        (*stream) >> currentChar;
    }

    void extractRunningSequence() { // Reading from a file
        assert(counter == 0);
        if (packedArray == 0 || currentPos >= (int)packedArray->size())
            return;
        counter = (int)packedArray->at(currentPos);
        currentChar = packedArray->at(currentPos + 1);
        currentPos += 2;
    }
};

#endif
