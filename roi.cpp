#include <deque>
#include <cassert>
#include <cstdlib>
//#include <QDataStream>
#include <fstream>
#include "roi.h"

void ROI::initialize(
    int max_slices,
    int w, int h
) {
    maxSlices = max_slices;
    sliceMin = 0;
    sliceMax = max_slices - 1;
    xMax = w;
    yMax = h;
    bitmasks.resize(max_slices);
    contours.resize(max_slices);
    for (int i = 0; i < maxSlices; ++i) {
        bitmasks.at(i).resizeBitmask(w, h);
        bitmasks.at(i).fill(ROI_POSITIVE);
        contours.at(i).clear();
    }
    framingRect = I2Rectangle(0, 0, w, h);
}

void ROI::clear() {
    for (int i = 0; i < maxSlices; ++i) {
        bitmasks.at(i).fill(ROI_POSITIVE);
        contours.at(i).clear();
    }
    framingRect = I2Rectangle(0, 0, xMax, yMax);
    sliceMin = 0;
    sliceMax = maxSlices - 1;
}

void Bitmask::drawLine(
    const I2Point& p0, const I2Point& p1, int v /* = ROI_POSITIVE */,
    bool thickLine /* = false */
) {
    assert(0 <= p0.x && p0.x < width());
    assert(0 <= p0.y && p0.y < height());
    assert(0 <= p1.x && p1.x < width());
    assert(0 <= p1.y && p1.y < height());

    if (
        (p0.x == p1.x && abs(p0.y - p1.y) <= 1) ||
        (p0.y == p1.y && abs(p0.x - p1.x) <= 1)
        ) {
        pixelRef(p0) = (unsigned char)v;
        if (thickLine) {
            if (p0.x > 0)
                pixelRef(p0.x - 1, p0.y) = (unsigned char)v;
            if (p0.x < width() - 1)
                pixelRef(p0.x + 1, p0.y) = (unsigned char)v;
            if (p0.y > 0)
                pixelRef(p0.x, p0.y - 1) = (unsigned char)v;
            if (p0.y < height() - 1)
                pixelRef(p0.x, p0.y + 1) = (unsigned char)v;
        }

        pixelRef(p1) = (unsigned char)v;
        if (thickLine) {
            if (p1.x > 0)
                pixelRef(p1.x - 1, p1.y) = (unsigned char)v;
            if (p1.x < width() - 1)
                pixelRef(p1.x + 1, p1.y) = (unsigned char)v;
            if (p1.y > 0)
                pixelRef(p1.x, p1.y - 1) = (unsigned char)v;
            if (p1.y < height() - 1)
                pixelRef(p1.x, p1.y + 1) = (unsigned char)v;
        }
        return;
    }

    pixelRef(p0.x, p0.y) = (unsigned char)v;
    if (thickLine) {
        if (p0.x > 0)
            pixelRef(p0.x - 1, p0.y) = (unsigned char)v;
        if (p0.x < width() - 1)
            pixelRef(p0.x + 1, p0.y) = (unsigned char)v;
        if (p0.y > 0)
            pixelRef(p0.x, p0.y - 1) = (unsigned char)v;
        if (p0.y < height() - 1)
            pixelRef(p0.x, p0.y + 1) = (unsigned char)v;
    }

    if (p0 == p1)
        return;
    int dx = p1.x - p0.x;
    int dy = p1.y - p0.y;
    if (abs(dx) >= abs(dy)) {
        assert(dx != 0);
        int x0, x1, y0, y1;
        if (dx > 0) {
            x0 = p0.x; y0 = p0.y;
            x1 = p1.x; y1 = p1.y;
        }
        else {
            x0 = p1.x; y0 = p1.y;
            x1 = p0.x; y1 = p0.y;
        }
        double deltaY = double(y1 - y0) / double(x1 - x0);
        assert(fabs(deltaY) <= 1.);
        int x = x0 + 1; double y = double(y0) + deltaY;
        while (true) {
            int yy = int(round(y));
            pixelRef(x, yy) = (unsigned char)v;
            if (thickLine) {
                if (x > 0)
                    pixelRef(x - 1, yy) = (unsigned char)v;
                if (x < width() - 1)
                    pixelRef(x + 1, yy) = (unsigned char)v;
                if (yy > 0)
                    pixelRef(x, yy - 1) = (unsigned char)v;
                if (yy < height() - 1)
                    pixelRef(x, yy + 1) = (unsigned char)v;
            }
            ++x;
            if (x > x1)
                break;
            y += deltaY;
        }
    }
    else {
        assert(dy != 0);
        int y0, y1, x0, x1;
        if (dy > 0) {
            y0 = p0.y; x0 = p0.x;
            y1 = p1.y; x1 = p1.x;
        }
        else {
            y0 = p1.y; x0 = p1.x;
            y1 = p0.y; x1 = p0.x;
        }
        double deltaX = double(x1 - x0) / double(y1 - y0);
        assert(fabs(deltaX) <= 1.);
        int y = y0 + 1; double x = double(x0) + deltaX;
        while (true) {
            int xx = int(round(x));
            pixelRef(xx, y) = (unsigned char)v;
            if (thickLine) {
                if (xx > 0)
                    pixelRef(xx - 1, y) = (unsigned char)v;
                if (xx < width() - 1)
                    pixelRef(xx + 1, y) = (unsigned char)v;
                if (y > 0)
                    pixelRef(xx, y - 1) = (unsigned char)v;
                if (y < height() - 1)
                    pixelRef(xx, y + 1) = (unsigned char)v;
            }
            ++y;
            if (y > y1)
                break;
            x += deltaX;
        }
    }
}

void Bitmask::drawContour(
    const I2Contour& c, int v /* = ROI_POSITIVE */,
    bool closed /* = true */, bool thickLine /* = false */
) {
    if (c.size() == 0)
        return;
    for (size_t i = 0; i < c.size() - 1; ++i) {
        drawLine(c.at(i), c.at(i + 1), v, thickLine);
    }
    if (closed)
        drawLine(c.back(), c.front(), v, thickLine);
}

void Bitmask::drawContour(
    const R2Contour& c, int v /* = ROI_POSITIVE */
) {
    if (c.size() == 0)
        return;
    I2Point p0(
        int(round(c.front().x)),
        int(round(c.front().y))
    );
    I2Point q0 = p0;
    for (size_t i = 1; i < c.size(); ++i) {
        I2Point q1(
            int(round(c.at(i).x)),
            int(round(c.at(i).y))
        );
        drawLine(q0, q1, v);
        q0 = q1;
    }
    drawLine(q0, p0, v);
}

void Bitmask::paintContour(
    const I2Contour& c,
    int v /* = ROI_POSITIVE */,
    int vBorder /* = ROI_POSITIVE */,
    int numNeighbours /* = 4 */,
    const I2Point* seed /* = 0 */,
    const Bitmask* externalBitmask /* = 0 */
) {
    drawContour(c, vBorder);
    if (c.size() <= 2)
        return;

    bool found = false;
    I2Point q;
    if (seed != 0) {
        if (c.containsStrictly(*seed)) {
            found = true;
            q = *seed;
        }
    }

    if (!found) {
        // Find a point inside a contour
        const I2Rectangle& framingRect = c.framingRect();
        q = I2Point(
            (framingRect.left() + framingRect.right()) / 2,
            (framingRect.top() + framingRect.bottom()) / 2
        );
        found = c.containsStrictly(q);
        if (!found) {
            int y = q.y;
            for (int x = q.x - 1; x > framingRect.left(); --x) {
                if (c.containsStrictly(I2Point(x, y))) {
                    q = I2Point(x, y);
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            int y = q.y;
            for (int x = q.x + 1; x < framingRect.right(); ++x) {
                if (c.containsStrictly(I2Point(x, y))) {
                    q = I2Point(x, y);
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            int x = q.x;
            for (int y = q.y - 1; y > framingRect.y; --y) {
                if (c.containsStrictly(I2Point(x, y))) {
                    q = I2Point(x, y);
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            int x = q.x;
            int ymax = framingRect.y + framingRect.height();
            for (int y = q.y + 1; y < ymax - 1; --y) {
                if (c.containsStrictly(I2Point(x, y))) {
                    q = I2Point(x, y);
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            for (size_t i = 0; i < c.size(); ++i) {
                size_t j = i + 1;
                if (j >= c.size())
                    j = 0;
                R2Point p0(c.at(i).x, c.at(i).y);
                R2Point p1(c.at(j).x, c.at(j).y);
                if (p0 == p1)
                    continue;
                R2Vector v = p1 - p0;
                R2Vector n = v.normal();
                if (fabs(n.x) >= fabs(n.y)) {
                    assert(n.x != 0.);
                    n.x /= fabs(n.x);   // n.x == +-1
                    n.y /= fabs(n.y);
                }
                else {
                    assert(n.y != 0.);
                    n.x /= fabs(n.x);   // n.x == +-1
                    n.y /= fabs(n.y);
                }
                R2Point p = p0 + v * 0.5 + n;
                q = I2Point(
                    int(round(p.x)), int(round(p.y))
                );
                if (c.containsStrictly(q)) {
                    found = true;
                    break;
                }
                p += n;
                q = I2Point(
                    int(round(p.x)), int(round(p.y))
                );
                if (c.containsStrictly(q)) {
                    found = true;
                    break;
                }
            } // end for
        }
    }

    if (!found)
        return; // Falure...

    if (numNeighbours == 4) {
        regionGrow4(q, v, externalBitmask);
    }
    else {
        assert(numNeighbours == 8);
        regionGrow8(q, v, externalBitmask);
    }
}

static const I2Vector NEIGHBOURS4[4] = {
    I2Vector(1, 0),
    I2Vector(-1, 0),
    I2Vector(0, 1),
    I2Vector(0, -1)
};

static const I2Vector NEIGHBOURS8[8] = {
    I2Vector(1, 0),
    I2Vector(-1, 0),
    I2Vector(0, 1),
    I2Vector(0, -1),
    I2Vector(1, 1),
    I2Vector(-1, 1),
    I2Vector(1, -1),
    I2Vector(-1, -1)
};

void Bitmask::regionGrow4(
    const I2Point& seed, int v /* = ROI_POSITIVE */,
    const Bitmask* externalBitmask /* = 0 */
) {
    std::deque<I2Point> deq;
    setPixValue(seed, v);
    deq.push_back(seed);
    while (!deq.empty()) {
        I2Point t = deq.front(); deq.pop_front();
        for (int i = 0; i < 4; ++i) { // For every neighbour
            I2Point n = t + NEIGHBOURS4[i];
            if (
                n.x <= 0 || n.x >= width() - 1 ||
                n.y <= 0 || n.y >= height() - 1
                )
                continue;
            if (
                externalBitmask != 0 &&
                externalBitmask->pixelAt(n) != 0
                )
                continue;

            if (pixelAt(n) == 0) {
                setPixValue(n, v);
                deq.push_back(n);
            }
        }
    }
}

void Bitmask::regionGrow8(
    const I2Point& seed, int v /* = ROI_POSITIVE */,
    const Bitmask* externalBitmask /* = 0 */
) {
    std::deque<I2Point> deq;
    setPixValue(seed, v);
    deq.push_back(seed);
    while (!deq.empty()) {
        I2Point t = deq.front(); deq.pop_front();
        for (int i = 0; i < 8; ++i) { // For every neighbour
            I2Point n = t + NEIGHBOURS8[i];
            if (
                n.x <= 0 || n.x >= width() - 1 ||
                n.y <= 0 || n.y >= height() - 1
                )
                continue;
            if (
                externalBitmask != 0 &&
                externalBitmask->pixelAt(n) != 0
                ) {
                int m = externalBitmask->pixelAt(n);
                if (
                    //... v != ROI_MANUAL_NEGATIVE_BORDER ||
                    m == ROI_MANUAL_NEGATIVE_BORDER
                    )
                    continue;
            }

            if (pixelAt(n) == 0) {
                setPixValue(n, v);
                deq.push_back(n);
            }
        }
    }
}

void Bitmask::paintContourInternalArea(
    const I2Contour& c,
    int v /* = ROI_POSITIVE */,
    int numNeghbours /* = 4 */,
    const I2Point* seed /* = 0 */,
    const Bitmask* externalBitmask /* = 0 */
) {
    Bitmask m(width(), height());
    // m.clear();
    m.paintContour(
        c, v, ROI_SPECIAL_BORDER_VALUE, numNeghbours,
        seed, externalBitmask
    );
    for (int y = 0; y < height(); ++y) {
        for (int x = 0; x < width(); ++x) {
            if (m.pixelAt(x, y) == v) {
                setPixValue(x, y, v);
            }
        }
    }
}

bool Bitmask::write(std::fstream& out) const {
    out << width();
    if (out.fail())
        return false;
    out << height();

    /*... Now we use RLE encoding
    for (int y = 0; y < height(); ++y) {
        for (int x = 0; x < width(); ++x) {
            out << (unsigned char) pixelAt(x, y);
        }
        if (out.status() == QDataStream::WriteFailed)
            return false;
    }
    ...*/

    RLEEncoder encoder(&out);
    for (int y = 0; y < height(); ++y) {
        for (int x = 0; x < width(); ++x) {
            encoder.writeByte((unsigned char)pixelAt(x, y));
        }
        if (out.fail())
            return false;
    }
    encoder.flush();
    if (out.fail())
        return false;

    return true;
}

bool Bitmask::read(std::fstream& in) {
    int w, h;
    in >> w >> h;
    if (!in.good())
        return false;
    resizeBitmask(w, h);
    clear();

    /*... Now we use RLE encoding
    for (int y = 0; y < height(); ++y) {
        for (int x = 0; x < width(); ++x) {
            unsigned char b;
            in >> b;
            pixelRef(x, y) = b;
        }
        if (in.status() != QDataStream::Ok)
            return false;
    }
    ...*/

    RLEEncoder encoder(&in);
    for (int y = 0; y < height(); ++y) {
        for (int x = 0; x < width(); ++x) {
            unsigned char b = encoder.readByte();
            pixelRef(x, y) = b;
        }
    }

    return true;
}

PackedBitmask::PackedBitmask(const Bitmask& m) :
    width(m.width()),
    height(m.height()),
    packedMatrix()
{
    RLEEncoder encoder(0, &packedMatrix);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            encoder.push_back(
                (unsigned char)m.pixelAt(x, y)
            );
        }
    }
    encoder.flush_array();
}

PackedBitmask& PackedBitmask::operator=(const Bitmask& m) {
    width = m.width();
    height = m.height();
    packedMatrix.clear();

    RLEEncoder encoder(0, &packedMatrix);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            encoder.push_back(
                (unsigned char)m.pixelAt(x, y)
            );
        }
    }
    encoder.flush_array();
    return *this;
}

Bitmask& PackedBitmask::unpack(Bitmask& m) const {
    m.resizeBitmask(width, height);
    RLEEncoder encoder(
        0,
        const_cast<std::vector<unsigned char>*>(&packedMatrix)
    );
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            m.pixelRef(x, y) = encoder.pop_front();
        }
    }
    return m;
}
