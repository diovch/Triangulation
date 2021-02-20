#ifndef R2GEOM_H
#define R2GEOM_H

#include <vector>
#include <utility>
#include <cmath>
#include "pi.h"

using std::size_t;

// const double R2_EPSILON = 1e-8;
const int MAX_INTERPOLATION_RAYS = 256;

// Predefinitions
class Bitmask;
class I2Contour;
class R2Contour;

template <class V = double> class G2Vector {
public:
    V x;
    V y;

    G2Vector(V px = V(), V py = V()):
        x(px),
        y(py)
    {}

    G2Vector<V> operator+(const G2Vector<V>& v) const {
        return G2Vector<V>(x + v.x, y + v.y);
    }

    G2Vector<V>& operator+=(const G2Vector<V>& v) {
        x += v.x;
        y += v.y;
        return *this;
    }

    G2Vector<V> operator-(const G2Vector<V>& v) const {
        return G2Vector<V>(x - v.x, y - v.y);
    }

    G2Vector<V>& operator-=(const G2Vector<V>& v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    G2Vector<V> operator*(V c) const {
        return G2Vector<V>(x*c, y*c);
    }

    V operator*(const G2Vector<V>& v) const { // Dot product
        return (x*v.x + y*v.y);
    }

    V signedArea(const G2Vector<V>& v) const {
        return (x*v.y - v.x*y);
    }

    V area(const G2Vector<V>& v) const {
        V s = x*v.y - v.x*y;
        if (s < 0)
            s = (-s);
        return s;
    }

    bool operator==(const G2Vector<V>& v) const {
        return (x == v.x && y == v.y);
    }

    bool operator!=(const G2Vector<V>& v) const {
        return !(*this == v);
    }

    bool operator<=(const G2Vector<V>& v) const {
        return (y < v.y || (y <= v.y && x <= v.x));
    }

    bool operator>(const G2Vector<V>& v) const {
        return !(*this <= v);
    }

    bool operator<(const G2Vector<V>& v) const {
        return (y < v.y || (y <= v.y && x < v.x));
    }

    bool operator>=(const G2Vector<V>& v) const {
        return !(*this < v);
    }

    double norm() const {
        return sqrt(double(x)*double(x) + double(y)*double(y));
    }

    G2Vector<V> normal() const {
        return G2Vector<V>(-y, x);
    }

    G2Vector<V> normalized() const {
        double l = norm();
        if (l <= 0.)
            return *this;
        return G2Vector<V>(V(x/l), V(y/l));
    }

    G2Vector<V>& normalize() const {
        double l = norm();
        if (l > 0.) {
            x = V(x/l);
            y = V(y/l);
        }
        return *this;
    }

    double angle(const G2Vector<V>& v) const {
        G2Vector<V> n = normal();
        double x = double(v*(*this));
        double y = double(v*n);
        return atan2(y, x);
    }
};

template <class V = double> class G2Point {
public:
    V x;
    V y;

    G2Point(V px = V(), V py = V()):
        x(px),
        y(py)
    {}

    G2Point<V> operator+(const G2Vector<V>& v) const {
        return G2Point<V>(x + v.x, y + v.y);
    }

    G2Point<V>& operator+=(const G2Vector<V>& v) {
        x += v.x;
        y += v.y;
        return *this;
    }

    G2Point<V> operator+(const G2Point<V>& v) const {
        return G2Point<V>(x + v.x, y + v.y);
    }

    G2Point<V>& operator+=(const G2Point<V>& v) {
        x += v.x;
        y += v.y;
        return *this;
    }

    G2Vector<V> operator-(const G2Point<V>& v) const {
        return G2Vector<V>(x - v.x, y - v.y);
    }

    G2Point<V> operator-(const G2Vector<V>& v) const {
        return G2Point<V>(x - v.x, y - v.y);
    }

    G2Point<V>& operator-=(const G2Vector<V>& v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    G2Point<V> operator*(V c) const {
        return G2Point<V>(x*c, y*c);
    }

    bool operator==(const G2Point<V>& v) const {
        return (x == v.x && y == v.y);
    }

    bool operator!=(const G2Point<V>& v) const {
        return !(*this == v);
    }

    bool operator<=(const G2Point<V>& v) const {
        return (y < v.y || (y <= v.y && x <= v.x));
    }

    bool operator>(const G2Point<V>& v) const {
        return !(*this <= v);
    }

    bool operator<(const G2Point<V>& v) const {
        return (y < v.y || (y == v.y && x < v.x));
    }

    bool operator>=(const G2Point<V>& v) const {
        return !(*this < v);
    }

    static double signedArea(
        const G2Point<V>& a,
        const G2Point<V>& b,
        const G2Point<V>& c
    ) {
        G2Vector<V> ab = b - a;
        G2Vector<V> ac = c - a;
        return double(ab.signedArea(ac))/2.;
    }

    static double area(
        const G2Point<V>& a,
        const G2Point<V>& b,
        const G2Point<V>& c
    ) {
        return fabs(G2Point<V>::signedArea(a, b, c));
    }

    double distance(const G2Point& p) const {
        return (p - *this).norm();
    }

    static double distance(const G2Point& p0, const G2Point& p1) {
        return p0.distance(p1);
    }
};

template <class V = double> class G2Rectangle {
public:
    V x;
    V y;
    V w;
    V h;

    G2Rectangle(
        V xx = V(), V yy = V(),
        V ww = V(), V hh = V()
    ):
        x(xx),
        y(yy),
        w(ww),
        h(hh)
    {}

    G2Rectangle(
        const G2Point<V>& p,
        V ww, V hh
    ):
        x(p.x),
        y(p.y),
        w(ww),
        h(hh)
    {}

    G2Rectangle(
        const G2Point<V>& left_bottom,
        const G2Point<V>& right_top
    ):
        x(left_bottom.x),
        y(left_bottom.y),
        w(right_top.x - left_bottom.x),
        h(right_top.y - left_bottom.y)
    {}

    G2Point<V> leftBottom() const {
        return G2Point<V>(x, y);
    }

    G2Point<V> rightTop() const {
        return G2Point<V>(x + w, y + h);
    }

    V width() const { return w; }
    V height() const { return h; }

    V left() const { return x; }
    V bottom() const { return y; }
    V right() const { return x + w; }
    V top() const { return y + h; }

    void setLeft(V l) { x = l; }
    void setBottom(V b) { y = b; }
    void setLeftBottom(const G2Point<V> p) {
        x = p.x; y = p.y;
    }
    void setWidth(V ww) { w = ww; }
    void setHeight(V hh) { h = hh; }

    bool operator==(const G2Rectangle<V>& r) const {
        return (x == r.x && y == r.y && w == r.w && h == r.h);
    }

    bool operator!=(const G2Rectangle<V>& r) const {
        return !(*this == r);
    }

    bool operator<(const G2Rectangle<V>& r) const {
        return (
            y < r.y || (
                y <= r.y && (
                    x < r.x || (
                        x <= r.x && (
                             w < r.w || (
                                w <= r.w && h < r.h
                             )
                        )
                    )
                )
            )
        );
    }

    bool operator<=(const G2Rectangle<V>& r) const {
        return (
            y < r.y || (
                y <= r.y && (
                    x < r.x || (
                        x <= r.x && (
                             w < r.w || (
                                w <= r.w && h <= r.h
                             )
                        )
                    )
                )
            )
        );
    }

    bool operator>(const G2Rectangle<V>& r) const {
        return !(*this <= r);
    }

    bool operator>=(const G2Rectangle<V>& r) const {
        return !(*this < r);
    }

    bool contains(const G2Point<V>& p) const {
        return (
            x <= p.x && p.x <= x + w &&
            y <= p.y && p.y <= y + h
        );
    }

    bool empty() const {
        return (w < V(0) || h < V(0));
    }

    V area() const {
        if (empty())
            return V(0);
        return w*h;
    }

    G2Rectangle<V> intersection(const G2Rectangle<V>& r) const {
        V x0 = left();
        if (r.left() > x0)
            x0 = r.left();
        V x1 = right();
        if (r.right() < x1)
            x1 = r.right();

        V y0 = bottom();
        if (r.bottom() > y0)
            y0 = r.bottom();
        V y1 = top();
        if (r.top() < y1)
            y1 = r.top();
        return G2Rectangle<V>(x0, y0, x1 - x0, y1 - y0);
    }

    G2Rectangle<V> combination(const G2Rectangle<V>& r) const {
        V x0 = left();
        if (r.left() < x0)
            x0 = r.left();
        V x1 = right();
        if (r.right() > x1)
            x1 = r.right();

        V y0 = bottom();
        if (r.bottom() < y0)
            y0 = r.bottom();
        V y1 = top();
        if (r.top() > y1)
            y1 = r.top();
        return G2Rectangle<V>(x0, y0, x1 - x0, y1 - y0);
    }
};

typedef G2Vector<int> I2Vector;
typedef G2Point<int> I2Point;
typedef G2Vector<double> R2Vector;
typedef G2Point<double> R2Point;
typedef G2Rectangle<int> I2Rectangle;
typedef G2Rectangle<double> R2Rectangle;

template <class T = double> class G2Contour:
    public std::vector< G2Point<T> > {
public:
    mutable bool areaComputed;
    mutable double signed_area;
    mutable bool rectComputed;
    mutable G2Rectangle<T> framing_rect;

    G2Contour():
        std::vector< G2Point<T> >(),
        areaComputed(false),
        signed_area(0.),
        rectComputed(false),
        framing_rect()
    {}

    G2Contour(int n):
        std::vector< G2Point<T> >(n),
        areaComputed(false),
        signed_area(0.),
        rectComputed(false),
        framing_rect()
    {}

    size_t size() const {
        return std::vector< G2Point<T> >::size();
    }

    void clear() {
        std::vector< G2Point<T> >::clear();
        areaComputed = false;
        rectComputed = false;
        signed_area = 0.;
    }

    void push_back(const G2Point<T>& p) {
        std::vector< G2Point<T> >::push_back(p);
        areaComputed = false;
        rectComputed = false;
    }

    void computeArea() const {
        G2Point<T> center;
        double s = 0.;
        for (size_t i = 0; i < size(); ++i) {
            size_t j = i+1;
            if (j >= size())
                j = 0;
            s += G2Point<T>::signedArea(center, (*this)[i], (*this)[j]);
        }
        signed_area = s;
        areaComputed = true;
    }

    double signedArea() const {
        if (!areaComputed)
            computeArea();
        return signed_area;
    }

    double area() const {
        return fabs(signedArea());
    }

    int orientation() const {
        double sa = signedArea();
        if (sa > 0.)
            return 1;
        else if (sa < 0.)
            return (-1);
        else
            return 0;
    }

    G2Contour<T>& orientate() {
        if (orientation() < 0) {
            size_t n2 = size()/2;
            for (int i = 0; i < int(n2); ++i) {
                std::swap(this->at(i), this->at(size() - 1 - i));
            }
        }
        return *this;
    }

    void computeFramingRect() const {
        if (size() == 0)
            return;
        T xmin = (*this)[0].x;
        T xmax = xmin;
        T ymin = (*this)[0].y;
        T ymax = ymin;

        for (size_t i = 0; i < size(); ++i) {
            T x = (*this)[i].x;
            T y = (*this)[i].y;
            if (x < xmin)
                xmin = x;
            if (x > xmax)
                xmax = x;
            if (y < ymin)
                ymin = y;
            if (y > ymax)
                ymax = y;
        }
        framing_rect.setLeft(xmin);
        framing_rect.setBottom(ymin);
        framing_rect.setWidth(xmax - xmin);
        framing_rect.setHeight(ymax - ymin);
        rectComputed = true;
    }

    G2Rectangle<T>& framingRect() {
        if (!rectComputed)
            computeFramingRect();
        return framing_rect;
    }

    const G2Rectangle<T>& framingRect() const {
        if (!rectComputed)
            computeFramingRect();
        return framing_rect;
    }

    bool onBorder(const G2Point<T>& q) const {
        for (int i = 0; i < int(size()); ++i) {
            int j = i + 1;
            if (j >= int(size()))
                j = 0;
            const G2Point<T>& p0 = this->at(i);
            const G2Point<T>& p1 = this->at(j);
            if (q == p0 || q == p1)
                return true;
            if (p0 != p1) {
                G2Vector<T> v = p1 - p0;
                G2Vector<T> n = v.normal();
                G2Vector<T> w = q - p0;
                if (
                    n*w == T(0) &&
                    w*v >= T(0) &&
                    w*w <= v*v
                ) {
                    return true;
                }
            }
        }
        return false;
    }

    bool contains(const G2Point<T>& p) const {
        if (size() <= 1)
            return false;
        if (!framingRect().contains(p))
            return false;

        double alpha = 0.;
        G2Vector<T> v0 = (*this)[0] - p;
        G2Vector<T> v_prev = v0;
        for (size_t i = 1; i < size(); ++i) {
            G2Vector<T> v = (*this)[i] - p;
            alpha += v_prev.angle(v);
            v_prev = v;
        }
        alpha += v_prev.angle(v0);

        return (fabs(alpha) >= PI); // Must be either 2*pi or 0
    }

    bool containsStrictly(const G2Point<T>& p) const {
        return (contains(p) && !onBorder(p));
    }

    R2Point centroid() const {
        if (size() == 0)
            return R2Point();
        double x = 0., y = 0.;
        for (size_t i = 0; i < size(); ++i) {
            x += double(this->at(i).x);
            y += double(this->at(i).y);
        }
        return R2Point(x/double(size()), y/double(size()));
    }
};

// typedef G2Contour<int> I2Contour;
// typedef G2Contour<double> R2Contour;

class R2Contour: public G2Contour<double> {
public:
    int contourType;
    mutable bool rayDistancesCalculated;
    mutable std::vector<double> rayDistances;

    static R2Vector ray_directions[MAX_INTERPOLATION_RAYS];
    static bool ray_directions_calculated;

public:
    R2Contour():
        contourType(0),
        rayDistancesCalculated(false)
    {}

    R2Contour(int n):
        G2Contour<double>(n),
        contourType(0),
        rayDistancesCalculated(false)
    {}

    R2Contour(const I2Contour& c);
    R2Contour& operator=(const I2Contour& c);

    // Interpolate contours that have star-like shapes.
    // 0 <= t <= 1. For t=0 the result equals *this,
    //              for t=1 the result equals c.
    R2Contour starInterpolation(const R2Contour& c, double t) const;
    void calculateRayDistances() const;

    static void calculate_ray_directions();
};

class I2Contour: public G2Contour<int> {
public:
    int contourType;

public:
    I2Contour():
        contourType(0)
    {}

    I2Contour(const R2Contour& c);
    I2Contour& operator=(const R2Contour& c);

    bool canAdd(const I2Point& p) const;
    bool canClose() const;
};

// Global functions

bool intersectStraightLines(
    const R2Point& p, const R2Vector& v,    // First line
    const R2Point& q, const R2Vector& w,    // Second line
    R2Point& intersection                   // Result
);

bool intersectLineSegmentAndLine(
    const R2Point& p0, const R2Point& p1,   // Line segment
    const R2Point& q, const R2Vector& v,    // Straight line
    R2Point& intersection                   // Result
);

bool intersectLineSegments(
    const R2Point& p0, const R2Point& p1,   // First line segment
    const R2Point& q0, const R2Point& q1,   // Second line segmant
    R2Point& intersection                   // Result
);

#endif
