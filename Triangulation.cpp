#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <clocale>
#include <cassert>
#include <set>
#include <deque>
#include "Triangulation.h"
#include "BinHeap.h"

using namespace R3Graph;

static const int TEST_MAX_TRIANGLES = 0;
static const int TEST_MAX_EDGES = 1000;

// Parsing of XML-file
static bool extractTag(FILE* f, int& tag, bool& closing);
static bool extractInt(FILE* f, int& value);
static bool extractDouble(FILE* f, double& value);
static int findTag(const char* tag);

static const char* const TAGS[] = {
    "triangulation",    // 0
    "vertices",         // 1
    "vertex",           // 2
    "coord",            // 3
    "normal",           // 4
    "triangles",        // 5
    "triangle",         // 6
    0
};

static const int TAG_TRIANGULATION = 0;
static const int TAG_VERTICES = 1;
static const int TAG_VERTEX = 2;
static const int TAG_COORD = 3;
static const int TAG_NORMAL = 4;
static const int TAG_TRIANGLES = 5;
static const int TAG_TRIANGLE = 6;

static bool writeTag(
    FILE* f, int tagID, int level, bool closing = false
);
static bool write3Int(
    FILE* f, int x, int y, int z
);
static bool write3Double(
    FILE* f, double x, double y, double z
);
static bool writeIndent(
    FILE* f, int level
);

void Triangulation::computeFramingBox() {
    if (vertices.size() == 0)
        return;
    R3Point minPoint = vertices[0].point;
    R3Point maxPoint = vertices[0].point;
    for (size_t i = 0; i < vertices.size(); ++i) {
        if (vertices[i].point.x < minPoint.x)
            minPoint.x = vertices[i].point.x;
        if (vertices[i].point.x > maxPoint.x)
            maxPoint.x = vertices[i].point.x;
        if (vertices[i].point.y < minPoint.y)
            minPoint.y = vertices[i].point.y;
        if (vertices[i].point.y > maxPoint.y)
            maxPoint.y = vertices[i].point.y;
        if (vertices[i].point.z < minPoint.z)
            minPoint.z = vertices[i].point.z;
        if (vertices[i].point.z > maxPoint.z)
            maxPoint.z = vertices[i].point.z;
    }
    box.origin = minPoint;
    box.size = maxPoint - minPoint;
}


void Triangulation::clear() {
    vertices.clear();
    triangles.clear();
    adjacentTrianglesCalculated = false;
    delete adjacentTriangles; adjacentTriangles = 0;
    trianglesOfVerticesCalculated = false;
    trianglesOfVertices.clear();
    linkedComponentsCalculated = false;
    delete linkedComponents; linkedComponents = 0;
}

void Triangulation::orientate() {
    for (size_t i = 0; i < triangles.size(); ++i) {
        Triangle& t = triangles.at(i);

        /*???
        Vertex v0 = vertices.at(t.indices[0]);
        Vertex v1 = vertices.at(t.indices[1]);
        Vertex v2 = vertices.at(t.indices[2]);
        R3Point p0 = v0.point;
        R3Point p1 = v1.point;
        R3Point p2 = v2.point;
        R3Vector meanN = (v0.normal + v1.normal + v2.normal)*(1./3.);
        R3Vector N = (p1 - p0).vectorProduct(p2 - p0);
        if (N.scalarProduct(meanN) < 0.) {
            t.invert();
        }
        ???*/

        t.invert();     //??? Nessary because of the error in Skala code:
                        //??? must be opposite orientation of triangles
    }
}


static const char * const xmlFirstLine =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";

bool Triangulation::save(const char* path) const {
    setlocale(LC_ALL, "C");

    FILE* f = fopen(path, "w");
    if (f == NULL) {
        perror("Cannot open a triangulationfile for writing");
        return false;
    }
    if (fprintf(f, "%s\n", xmlFirstLine) <= 0) {
        LError: ;
        perror("Write error");
        fclose(f);
        return false;
    }
    int level = 0;
    if (!writeTag(f, TAG_TRIANGULATION, level, false))
        goto LError;
    ++level;

    if (!writeTag(f, TAG_VERTICES, level, false))
        goto LError;
    ++level;
    for (size_t i = 0; i < vertices.size(); ++i) {
        const Vertex& v = vertices[i];

        if (!writeTag(f, TAG_VERTEX, level, false))
            goto LError;
        ++level;

        if (
            !writeIndent(f, level) ||
            !writeTag(f, TAG_COORD, (-1), false)
        )
            goto LError;
        if (!write3Double(f, v.point.x, v.point.y, v.point.z))
            goto LError;
        if (!writeTag(f, TAG_COORD, (-1), true))
            goto LError;

        if (
            !writeIndent(f, level) ||
            !writeTag(f, TAG_NORMAL, (-1), false)
        )
            goto LError;
        if (!write3Double(f, v.normal.x, v.normal.y, v.normal.z))
            goto LError;
        if (!writeTag(f, TAG_NORMAL, (-1), true))
            goto LError;

        --level;
        if (!writeTag(f, TAG_VERTEX, level, true))
            goto LError;
    }
    --level;
    if (!writeTag(f, TAG_VERTICES, level, true))
        goto LError;

    if (!writeTag(f, TAG_TRIANGLES, level, false))
        goto LError;
    ++level;
    for (size_t i = 0; i < triangles.size(); ++i) {
        const Triangle& t = triangles[i];

        if (
            !writeIndent(f, level) ||
            !writeTag(f, TAG_TRIANGLE, (-1), false)
        )
            goto LError;
        if (!write3Int(f, t.indices[0], t.indices[1], t.indices[2]))
            goto LError;
        if (!writeTag(f, TAG_TRIANGLE, (-1), true))
            goto LError;
    }
    --level;
    if (!writeTag(f, TAG_TRIANGLES, level, true))
        goto LError;

    --level;
    if (!writeTag(f, TAG_TRIANGULATION, level, true))
        goto LError;
    assert(level == 0);
    fclose(f);
    return true;
}

// level defines identation
// if level < 0, then continue the current line and
//               write '\n' for closing tag
static bool writeTag(
    FILE* f, int tagID, int level, bool closing // = false
) {
    if (level >= 0) {
        if (!writeIndent(f, level))
            return false;
    }
    if (fprintf(f, "<") <= 0)
        return false;
    if (closing) {
        if (fprintf(f, "/") <= 0)
            return false;
    }
    if (fprintf(f, "%s>", TAGS[tagID]) <= 0)
        return false;
    if (level >= 0 || closing) {
        if (fprintf(f, "\n") <= 0)
            return false;
    }
    return true;
}

static bool write3Int(
    FILE* f, int x, int y, int z
) {
    if (fprintf(f, "%d %d %d", x, y, z) <= 0)
        return false;
    return true;
}

static bool write3Double(
    FILE* f, double x, double y, double z
) {
    if (fprintf(f, "%f %f %f", x, y, z) <= 0)
        return false;
    return true;
}

static bool writeIndent(
    FILE* f, int level
) {
    // Indentation
    for (int i = 0; i < level; ++i) {
        if (fprintf(f, " ") <= 0)
            return false;
    }
    return true;
}

static int findTag(const char* tag) {
    const char* line = TAGS[0];
    int id = 0;
    while (line != 0) {
        if (strcmp(line, tag) == 0)
            return id;
        ++id;
        line = TAGS[id];
    }
    return (-1);
}

static bool extractTag(FILE* f, int& tagID, bool& closing) {
    // Using Finite State Machine
    char tag[64];
    int tagLen = 0;
    int state = 0;
    bool tagEndFound = false;
    closing = false;
    int c;
    while (!tagEndFound) {
        c = fgetc(f);
        if (c == EOF)
            return false;
        if (state == 0) { // Looking for the tag beginning
            if (c != '<') {
                continue;
            } else {
                state = 1; // Tag beginning is found
                closing = false;
                tagLen = 0;
            }
        } else if (state == 1) { // Tag beginning "<" is found
            if (c == '/') {
                closing = true;
                state = 2; // Name must follow
            } else if (isalpha(c)) {
                tag[tagLen] = (char) c;
                ++tagLen;
                state = 2; // Name beginning is found
            } else {
                state = 0; // Go to initial state
                closing = false;
                tagLen = 0;
            }
        } else if (state == 2) { // Reading tag name
            if (isalpha(c)) {
                if (tagLen < 62) {
                    tag[tagLen] = (char) c;
                    ++tagLen;
                }
            } else {
                if (c == '>') {
                    if (tagLen > 0) {
                        tagEndFound = true;
                        break;
                    } else {
                        state = 0; // Go to initial state
                        closing = false;
                        tagLen = 0;
                    }
                }
                state = 3;  // Name end is found
            }
        } else if (state == 3) { // Looking for the tag end ">"
            if (c == '>') {
                if (tagLen > 0) {
                    tagEndFound = true;
                    break;
                } else {
                    state = 0; // Go to initial state
                    closing = false;
                    tagLen = 0;
                }
            }
        }
    }
    if (!tagEndFound || tagLen == 0)
        return false;
    tag[tagLen] = 0;
    int id = findTag(tag);
    if (id < 0)
        return false;
    tagID = id;
    return true;
}

static bool extractInt(FILE* f, int& value) {
    char line[64];
    int lineLen = 0;
    int c;
    while (true) {
        c = fgetc(f);
        if (c == EOF)
            return false;
        else if (isdigit(c) || c == '+' || c == '-') {
            break;
        }
    }

    // The beginning of integer constant is found
    while (true) {
        if (lineLen < 62) {
            line[lineLen] = (char) c;
            ++lineLen;
        }
        c = fgetc(f);
        if (c == EOF)
            break;
        if (!isdigit(c)) {
            break;
        }
    }
    if (c != EOF)
        ungetc(c, f);
    line[lineLen] = 0;
    value = atoi(line);
    return true;
}

static bool extractDouble(FILE* f, double& value) {
    char line[64];
    int lineLen = 0;
    int c;
    while (true) {
        c = fgetc(f);
        if (c == EOF)
            return false;
        else if (
            isdigit(c) ||
            c == '+' || c == '-' ||
            c == '.'
        ) {
            break;
        }
    }

    // The beginning of double constant is found
    while (true) {
        if (lineLen < 62) {
            line[lineLen] = (char) c;
            ++lineLen;
        }
        c = fgetc(f);
        if (c == EOF)
            break;
        if (
            !isdigit(c) &&
            c != '+' && c != '-' &&
            c != '.' && c != 'e' && c != 'E'
        ) {
            break;
        }
    }
    if (c != EOF)
        ungetc(c, f);
    line[lineLen] = 0;
    value = atof(line);
    return true;
}

bool Triangulation::load(const char *path) {
    setlocale(LC_ALL, "C");

    clear();
    FILE* f = fopen(path, "r");
    if (f == NULL)
        return false;
    int tagID = (-1);
    bool closing = false;
    while (
        extractTag(f, tagID, closing) &&
        tagID != TAG_TRIANGULATION
    );
    if (tagID != TAG_TRIANGULATION) {
        fclose(f);
        return false;
    }

    while (
        extractTag(f, tagID, closing) &&
        tagID != TAG_VERTICES
    );
    if (tagID != TAG_VERTICES || closing) {
        fclose(f);
        return false;
    }

    // Load vertices
    vertices.clear();
    if (
        !extractTag(f, tagID, closing) ||
        tagID != TAG_VERTEX || closing
    ) {
        fclose(f);
        return false;
    }

    while (true) {
        if (
            !extractTag(f, tagID, closing) ||
            tagID != TAG_COORD || closing
        ) {
            fclose(f);
            return false;
        }

        double x, y, z;
        if (!extractDouble(f, x)) {
            fclose(f);
            return false;
        }
        if (!extractDouble(f, y)) {
            fclose(f);
            return false;
        }
        if (!extractDouble(f, z)) {
            fclose(f);
            return false;
        }
        if (
            !extractTag(f, tagID, closing) ||
            tagID != TAG_COORD || !closing
        ) {
            fclose(f);
            return false;
        }

        if (
            !extractTag(f, tagID, closing) ||
            tagID != TAG_NORMAL || closing
        ) {
            fclose(f);
            return false;
        }

        double nx, ny, nz;
        if (!extractDouble(f, nx)) {
            fclose(f);
            return false;
        }
        if (!extractDouble(f, ny)) {
            fclose(f);
            return false;
        }
        if (!extractDouble(f, nz)) {
            fclose(f);
            return false;
        }
        if (
            !extractTag(f, tagID, closing) ||
            tagID != TAG_NORMAL || !closing
        ) {
            fclose(f);
            return false;
        }

        if (
            !extractTag(f, tagID, closing) ||
            tagID != TAG_VERTEX || !closing
        ) {
            fclose(f);
            return false;
        }

        vertices.push_back(
            Triangulation::Vertex(
                R3Point(x, y, z),
                //... R3Vector(nx, ny, nz).normalize()
                R3Vector(nx, ny, nz)
            )
        );

        if (!extractTag(f, tagID, closing)) {
            fclose(f);
            return false;
        }
        if (tagID == TAG_VERTEX && !closing)
            continue;
        break;
    } // end while

    // Load triangles
    if (
        !extractTag(f, tagID, closing) ||
        tagID != TAG_TRIANGLES || closing
    ) {
        fclose(f);
        return false;
    }

    triangles.clear();
    if (
        !extractTag(f, tagID, closing) ||
        tagID != TAG_TRIANGLE || closing
    ) {
        fclose(f);
        return false;
    }

    while (true) {
        int v0, v1, v2;
        if (!extractInt(f, v0)) {
            fclose(f);
            return false;
        }
        if (!extractInt(f, v1)) {
            fclose(f);
            return false;
        }
        if (!extractInt(f, v2)) {
            fclose(f);
            return false;
        }
        if (
            !extractTag(f, tagID, closing) ||
            tagID != TAG_TRIANGLE || !closing
        ) {
            fclose(f);
            return false;
        }

        triangles.push_back(
            Triangulation::Triangle(v0, v1, v2)
        );

        if (!extractTag(f, tagID, closing)) {
            fclose(f);
            return false;
        }
        if (tagID == TAG_TRIANGLE && !closing)
            continue;
        break;
    } // end while
    fclose(f);

    computeFramingBox();
    return true;
}


//--------------------------------------------------------------------------------------------------------------
R3Point Triangulation::center() const {
    return (box.origin + box.size * 0.5);
}
//
void Triangulation::computeNormals() {
    if (!trianglesOfVerticesCalculated) {
        defineTrianglesOfVertices();
    }
    for (size_t i = 0; i < vertices.size(); ++i) {
        Vertex& v = vertices.at(i);
        const TrianglesOfVertex& triangs = trianglesOfVertices.at(i);
        R3Vector sum(0., 0., 0.);
        for (size_t t = 0; t < triangs.size(); ++t) {
            const Triangle& tr = triangles.at(triangs[t]);
            R3Point p0 = vertices.at(tr[0]).point;
            R3Point p1 = vertices.at(tr[1]).point;
            R3Point p2 = vertices.at(tr[2]).point;
            R3Vector v1 = p1 - p0; v1.normalize();
            R3Vector v2 = p2 - p0; v2.normalize();
            R3Vector n = v1.vectorProduct(v2);
            n.normalize();
            sum += n;
        }
        if (triangs.size() > 0) {
            sum *= 1. / (double)triangs.size();
        }
        v.normal = sum;
    }
}//
//void Triangulation::refine() {  // Remove double vertices
//    Triangulation t;
//    std::map<R3Point, int> pointIdx;
//    std::map<int, int> newVertexIdx;
//    int numNewVertices = 0;
//    for (int i = 0; i < (int)vertices.size(); ++i) {
//        const Vertex& v = vertices[i];
//        if (pointIdx.count(v.point) == 0) {
//            t.vertices.push_back(v);
//            pointIdx[v.point] = numNewVertices;
//            newVertexIdx[i] = numNewVertices;
//            ++numNewVertices;
//        }
//        else {
//            // Double point found
//            newVertexIdx[i] = pointIdx[v.point];
//        }
//    }
//    //Q_ASSERT(numNewVertices == (int)t.vertices.size());
//
//    //qDebug() << "refine(): vertices=" << vertices.size() <<
//    //    " different vertices=" << numNewVertices << "\n";
//
//    //if (numNewVertices == (int)vertices.size()) {
//    //    qDebug() << "refine(): no double vertices!\n";
//    //    return;
//    //}
//
//    std::set<Triangle> triangleSet;
//    for (int i = 0; i < (int)triangles.size(); ++i) {
//        const Triangle& tr = triangles[i];
//        Triangle newTriangle(
//            newVertexIdx[tr[0]],
//            newVertexIdx[tr[1]],
//            newVertexIdx[tr[2]]
//        );
//        if (triangleSet.count(newTriangle) > 0)
//            continue;   // The triangle is already presented, skip it
//        t.triangles.push_back(newTriangle);
//    }
//    /*qDebug() << "refine(): triangles=" << triangles.size() <<
//        " different triangles=" << t.triangles.size() << "\n";*/
//
//    t.computeFramingBox();
//    clear();    //???
//    *this = t;
//}

void Triangulation::defineAdjacentTriangles() const {
    invalidateAdjacentTriangles();
    if (adjacentTriangles == 0)
        adjacentTriangles = new std::vector<AdjacentTriangles>;
    adjacentTriangles->resize(triangles.size());

    // Adjacent triangles for each triangle
    std::map<Edge, int> triangleOfEdge;
    for (size_t i = 0; i < triangles.size(); ++i) {
        const Triangle& triangle = triangles[i];
        for (int j = 0; j < 3; ++j) {
            int k = j + 1;
            if (k >= 3)
                k = 0;
            Edge e(triangle[j], triangle[k]);
            if (triangleOfEdge.count(e) == 0) {
                triangleOfEdge[e] = i;
            }
            else {
                adjacentTriangles->at(triangleOfEdge[e]).add(i);
                adjacentTriangles->at(i).add(triangleOfEdge[e]);
            }
        }
    }

    adjacentTrianglesCalculated = true;
}

void Triangulation::defineTrianglesOfVertices() const {
    trianglesOfEdges.clear();

    trianglesOfVertices.resize(vertices.size());
    for (size_t i = 0; i < trianglesOfVertices.size(); ++i)
        trianglesOfVertices.at(i).clear();

    // Adjacent triangles for each vertex
    for (size_t i = 0; i < triangles.size(); ++i) {
        const Triangle& triangle = triangles.at(i);
        for (int j = 0; j < 3; ++j) {
            int v = triangle[j];
            trianglesOfVertices.at(v).push_back(i);

            int w = triangle[(j + 1) % 3];
            Edge e(v, w);
            trianglesOfEdges[e].push_back(i);
        }
    }
    trianglesOfVerticesCalculated = true;

    // starsOfVertices.clear();
    starsOfVertices.resize(vertices.size());
    for (int v = 0; v < int(vertices.size()); ++v) {
        starsOfVertices.at(v).clear();
        const TrianglesOfVertex& vertexTriangles = trianglesOfVertices.at(v);
        for (int t = 0; t < int(vertexTriangles.size()); ++t) {
            int triangleIdx = vertexTriangles[t];
            const Triangle& triangle = triangles.at(triangleIdx);
            for (int j = 0; j < 3; ++j) {
                int w = triangle[j];
                if (w != v) {
                    starsOfVertices.at(v).insert(w);
                }
            }
        }
    }

    ringsOfVertices.resize(vertices.size());
    for (size_t v = 0; v < vertices.size(); ++v) {
        computeVertexRing(v, ringsOfVertices.at(v));
    }
}

//int Triangulation::defineLinkedComponents() const {
//    if (!adjacentTrianglesCalculated)
//        defineAdjacentTriangles();
//
//    clearLinkedComponents();
//
//    if (linkedComponents == 0)
//        linkedComponents = new std::vector<LinkedComponent>;
//
//    int numTriangles = (int)triangles.size();
//    if (numTriangles == 0)
//        return 0;
//
//    BinHeapMin<char> heap(numTriangles);
//
//    int* triangleIndex = new int[numTriangles];
//    int* heapIndex = new int[numTriangles];
//    memset(heapIndex, (char)(-1), numTriangles * sizeof(int));
//
//    int infinity = 4;
//
//    int i = 0;  // Index in heap
//    for (int t = 0; t < numTriangles; ++t) {
//        heap.elements[i] = (char)infinity;    // Infinite distance
//        triangleIndex[i] = t;
//        heapIndex[t] = i;
//        ++i;
//    }
//
//    heap.indexArray = triangleIndex;
//    heap.heapIndex = heapIndex;
//    heap.numElems = numTriangles;
//
//    int initialTriangle = 0;
//    int numComponents = 0;
//    int totalVisited = 0;
//
//    while (heap.numElems > 0) {
//        // Add new empty component to the list
//        linkedComponents->push_back(LinkedComponent());
//
//        int initialTriangleIndexInHeap = heap.heapIndex[initialTriangle];
//        heap.elements[initialTriangleIndexInHeap] = 0;  // Zero distance
//        heap.bubbleUp(initialTriangleIndexInHeap);
//
//        int numVisited = 0;
//
//        while (heap.numElems > 0 && heap.root() < infinity) {
//            //... int currentDistance = heap.root();
//            int current = heap.rootIndex();
//            heapIndex[current] = (-1);  // Remove triangle from heap
//            heap.removeRoot();
//            linkedComponents->back().push_back(current);
//            ++numVisited;
//
//            if (heap.numElems == 0)
//                break;
//
//            //... int newDistance = currentDistance + 1;
//            for (int j = 0; j < 3; ++j) {
//                int neighbor =
//                    adjacentTriangles->at(current).adjacentTriangles[j];
//                if (neighbor < 0)
//                    break;
//                int neighborIndexInHeap = heapIndex[neighbor];
//                if (neighborIndexInHeap < 0)    // Already visited
//                    continue;
//
//                if (heap.elements[neighborIndexInHeap] != 0) {
//                    heap.elements[neighborIndexInHeap] = 0; // Mark as visited
//                    heap.bubbleUp(neighborIndexInHeap);
//                }
//            }
//        } // end while (heap.root() < infinity)
//
//        totalVisited += numVisited;
//        ++numComponents;
//
//        //qDebug() <<
//        //    "Linked component " << numComponents <<
//        //    " is extracted: numTriangles=" << numVisited <<
//        //    " totalVisited=" << totalVisited << "\n";
//
//        // Define the first unvisited triangle
//        if (heap.numElems > 0)
//            initialTriangle = heap.rootIndex();
//        else
//            initialTriangle = (-1);
//    } // end while
//
//    //qDebug() <<
//    //    "Total linked components: " << numComponents <<
//    //    ", triangles: " << totalVisited << "\n";
//
//    //Q_ASSERT(totalVisited == numTriangles);
//    //Q_ASSERT(numComponents == (int)linkedComponents->size());
//
//    delete[] heapIndex;
//    delete[] triangleIndex;
//
//    linkedComponentsCalculated = true;
//    return numComponents;
//}

//void Triangulation::copyMaximalComponent(Triangulation& t) const {
//    t.clear();
//    //Q_ASSERT(linkedComponentsCalculated && linkedComponents != 0);
//    if (!linkedComponentsCalculated || linkedComponents->size() == 0)
//        return;
//
//    // Find maximal component
//    int indMax = 0;
//    int numMax = (int)linkedComponents->at(0).size();
//    for (int i = 1; i < (int)linkedComponents->size(); ++i) {
//        if ((int)linkedComponents->at(i).size() > numMax) {
//            indMax = i;
//            numMax = (int)linkedComponents->at(i).size();
//        }
//    }
//    //qDebug() << "Maximal linked component: " << numMax << " triangles.\n";
//
//    t.triangles.resize(numMax);
//
//    std::vector<int> vertexIndices(vertices.size());
//    for (size_t i = 0; i < vertices.size(); ++i)
//        vertexIndices[i] = (-1);
//
//    const LinkedComponent& maxComp = linkedComponents->at(indMax);
//    LinkedComponent::const_iterator i = maxComp.begin();
//    int numTriangles = 0, numVert = 0;
//    while (i != maxComp.end()) {
//        int triangleIdx = *i;
//        const Triangle& triangle = triangles.at(triangleIdx);
//        Triangle newTriangle;
//        for (int j = 0; j < 3; ++j) {
//            int v = triangle[j];
//            if (vertexIndices[v] < 0) {
//                t.vertices.push_back(vertices[v]);
//                vertexIndices[v] = numVert;
//                ++numVert;
//            }
//            newTriangle[j] = vertexIndices[v];
//        }
//        t.triangles.at(numTriangles) = newTriangle;
//        ++numTriangles;
//        ++i;
//    }
//    //qDebug() << "Maximal linked component is copied.\n";
//}

//bool Triangulation::contains(const R3Point& p) const {
//    double s = 0.;
//    for (int i = 0; i < (int)triangles.size(); ++i) {
//        const Triangle& tr = triangles[i];
//        double a = (R3Vector ()).signedSolidAngle(
//            vertices[tr[0]].point - p,
//            vertices[tr[1]].point - p,
//            vertices[tr[2]].point - p
//        );
//        s += a;
//    }
//    s = fabs(s);        // Safety
//    // s must be 4*pi for a point inside and 0 for a point outside.
//    return (s > 3. * PI);
//}

void Triangulation::computeVertexRing(
    int vertexIdx, VertexRing& vertexRing
) const {
    assert(trianglesOfVerticesCalculated);
    assert(trianglesOfVertices.size() == vertices.size());

    std::multimap<int, int> ringEdges;

    const TrianglesOfVertex& vertexTriangles =
        trianglesOfVertices.at(vertexIdx);
    int initialVertex = (-1);
    for (size_t t = 0; t < vertexTriangles.size(); ++t) {
        int triangleIdx = vertexTriangles[t];
        const Triangle& triangle = triangles.at(triangleIdx);
        int k = 0;
        int triangleVertices[2];
        for (int j = 0; j < 3; ++j) {
            int w = triangle[j];
            if (w != vertexIdx) {
                assert(k < 2);
                triangleVertices[k] = w;
                ++k;
                if (initialVertex < 0)
                    initialVertex = w;
            }
        }
        assert(k == 2);
        ringEdges.insert(std::make_pair(
            triangleVertices[0], triangleVertices[1]
        ));
        ringEdges.insert(std::make_pair(
            triangleVertices[1], triangleVertices[0]
        ));
    }

    assert(initialVertex >= 0);
    std::deque<int> ring;
    ring.push_back(initialVertex);
    size_t ringSize = 1;
    while (true) {
        // Growing at the end of ring
        auto range = ringEdges.equal_range(ring.back());
        auto i = range.first;
        while (i != range.second) {
            int vNext = i->second;
            if (
                ring.size() <= 1 ||
                (
                    vNext != ring[ring.size() - 2] &&
                    vNext != ring.front()
                    )
                ) {
                ring.push_back(vNext);
                break;
            }
            ++i;
        }

        // Growing at the beginning of ring
        range = ringEdges.equal_range(ring.front());
        i = range.first;
        while (i != range.second) {
            int vPrev = i->second;
            if (
                ring.size() <= 1 ||
                (
                    vPrev != ring[1] &&
                    vPrev != ring.back()
                    )
                ) {
                ring.push_front(vPrev);
                break;
            }
            ++i;
        }

        if (ring.size() == ringSize)    // No change =>
            break;                      //     end of ring construction
        ringSize = ring.size(); // Continue the ring construstion
    }

    vertexRing.resize(ring.size());
    for (size_t i = 0; i < ring.size(); ++i) {
        vertexRing[i] = ring[i];
    }

    //assert(
    //    vertexRing.size() >= vertexTriangles.size()
    //);

    vertexRing.borderVertex = (
        vertexRing.size() > vertexTriangles.size()
        );
}


//void Triangulation::checkBorderEdges(
//    int& numManifoldEdges,
//    int& numBorderEdges,
//    int& numNonManifoldEdges // Edges that belong to more then 2 triangles
//) const {
//    numManifoldEdges = 0;
//    numBorderEdges = 0;
//    numNonManifoldEdges = 0;
//
//    std::map<Edge, int> edgeDegree;
//    for (size_t i = 0; i < triangles.size(); ++i) {
//        const Triangulation::Triangle& t = triangles.at(i);
//        for (int i = 0; i < 3; ++i) {
//            int v0 = t[i];
//            int v1;
//            if (i < 2)
//                v1 = t[i + 1];
//            else
//                v1 = t[0];
//            Edge te(v0, v1);
//            if (edgeDegree.find(te) == edgeDegree.end()) {
//                edgeDegree[te] = 1;
//            }
//            else {
//                ++(edgeDegree[te]);
//            }
//        }
//    }
//
//    std::map<Edge, int>::const_iterator i = edgeDegree.cbegin();
//    while (i != edgeDegree.cend()) {
//        assert(i->second > 0);
//        if (i->second == 2) {
//            ++numManifoldEdges;
//        }
//        else if (i->second == 1) {
//            ++numBorderEdges;
//        }
//        else {
//            ++numNonManifoldEdges;
//        }
//
//        ++i;
//    }
//}

//void Triangulation::annihilateSmallTriangles(
//    int& numRemovedEdges,
//    int& numRemovedCompleteTriangles,
//    int& numRemovedTriangles,
//    double eps /* = 0.01 */
//) {
//    std::vector<Vertex> newVertices;
//    std::vector<Triangle> newTriangles;
//
//    std::map<int, int> replacedVertex;
//    std::vector<Vertex> modifiedVertices;
//    std::set<int> removedTriangles;
//    std::set<Edge> removedEdges;
//
//    numRemovedEdges = 0;
//    numRemovedCompleteTriangles = 0;
//
//    for (size_t i = 0; i < triangles.size(); ++i) {
//        const Triangle& t = triangles.at(i);
//
//        int vertexIdx0 = t.indices[0];
//        int vertexIdx1 = t.indices[1];
//        int vertexIdx2 = t.indices[2];
//
//        Edge edge0(vertexIdx0, vertexIdx1);
//        Edge edge1(vertexIdx1, vertexIdx2);
//        Edge edge2(vertexIdx2, vertexIdx0);
//
//        Vertex v0 = vertices.at(vertexIdx0);
//        Vertex v1 = vertices.at(vertexIdx1);
//        Vertex v2 = vertices.at(vertexIdx2);
//
//        R3Point p0 = v0.point;
//        R3Point p1 = v1.point;
//        R3Point p2 = v2.point;
//
//        double len[3];
//        len[0] = p0.distance(p1);
//        len[1] = p1.distance(p2);
//        len[2] = p2.distance(p0);
//
//        if (
//            len[0] <= eps &&
//            len[1] <= eps &&
//            len[2] <= eps &&
//            replacedVertex.count(vertexIdx0) == 0 &&
//            replacedVertex.count(vertexIdx1) == 0 &&
//            replacedVertex.count(vertexIdx2) == 0
//
//            && numRemovedCompleteTriangles < TEST_MAX_TRIANGLES //???
//            )
//        {
//            // Replace complete triangle by its center
//            R3Point center(
//                (p0.x + p1.x + p2.x) / 3.,
//                (p0.y + p1.y + p2.y) / 3.,
//                (p0.z + p1.z + p2.z) / 3.
//            );
//            R3Vector normal = v0.normal + v1.normal + v2.normal;
//            normal.normalize();
//            modifiedVertices.push_back(
//                Vertex(center, normal)
//            );
//            int replacedVertexIdx = int(modifiedVertices.size() - 1);
//            replacedVertex[vertexIdx0] = replacedVertexIdx;
//            replacedVertex[vertexIdx1] = replacedVertexIdx;
//            replacedVertex[vertexIdx2] = replacedVertexIdx;
//            ++numRemovedCompleteTriangles;
//
//            removedEdges.insert(edge0);
//            removedEdges.insert(edge1);
//            removedEdges.insert(edge2);
//
//            //... removedTriangles.insert(i); // Remove this triangle
//            continue;
//        }
//
//        int numSmallEdges = 0;
//        int smallEdgeIdx = (-1);
//        if (
//            len[0] <= eps &&
//            replacedVertex.count(vertexIdx0) == 0 &&
//            replacedVertex.count(vertexIdx1) == 0
//            ) {
//            smallEdgeIdx = 0;
//            ++numSmallEdges;
//        }
//        if (
//            len[1] <= eps &&
//            replacedVertex.count(vertexIdx1) == 0 &&
//            replacedVertex.count(vertexIdx2) == 0 &&
//            (
//                smallEdgeIdx < 0 ||
//                len[1] < len[smallEdgeIdx]
//                )
//            ) {
//            smallEdgeIdx = 1;
//            ++numSmallEdges;
//        }
//        if (
//            len[2] <= eps &&
//            replacedVertex.count(vertexIdx2) == 0 &&
//            replacedVertex.count(vertexIdx0) == 0 &&
//            (
//                smallEdgeIdx < 0 ||
//                len[2] < len[smallEdgeIdx]
//                )
//            ) {
//            smallEdgeIdx = 2;
//            ++numSmallEdges;
//        }
//        if (numSmallEdges == 0)
//            continue;
//
//        if (numRemovedEdges >= TEST_MAX_EDGES)
//            continue;           //?????
//
//        int vertIdx0 = t.indices[smallEdgeIdx];
//        int vertIdx1 = t.indices[(smallEdgeIdx + 1) % 3];
//        Edge smallEdge(vertIdx0, vertIdx1);
//
//        assert(removedEdges.count(smallEdge) == 0);
//        assert(
//            replacedVertex.count(vertIdx0) == 0 &&
//            replacedVertex.count(vertIdx1) == 0
//        );
//
//        // Replace 2 vertices by the center of edge
//        R3Point center =
//            vertices.at(vertIdx0).point +
//            (vertices.at(vertIdx1).point - vertices.at(vertIdx0).point) * 0.5;
//        R3Vector normal = vertices.at(vertIdx0).normal +
//            vertices.at(vertIdx1).normal;
//        normal.normalize();
//        modifiedVertices.push_back(
//            Vertex(center, normal)
//        );
//        int replacedVertexIdx = int(modifiedVertices.size() - 1);
//        replacedVertex[vertIdx0] = replacedVertexIdx;
//        replacedVertex[vertIdx1] = replacedVertexIdx;
//
//        removedEdges.insert(smallEdge);
//        ++numRemovedEdges;
//
//        //... removedTriangles.insert(i); // Remove this triangle
//    } // end for
//
//    numRemovedTriangles = 0;
//    for (size_t i = 0; i < triangles.size(); ++i) {
//        const Triangle& t = triangles.at(i);
//
//        int vertexIdx0 = t.indices[0];
//        int vertexIdx1 = t.indices[1];
//        int vertexIdx2 = t.indices[2];
//
//        Edge edge0(vertexIdx0, vertexIdx1);
//        Edge edge1(vertexIdx1, vertexIdx2);
//        Edge edge2(vertexIdx2, vertexIdx0);
//
//        if (
//            removedEdges.count(edge0) != 0 ||
//            removedEdges.count(edge1) != 0 ||
//            removedEdges.count(edge2) != 0
//            ) {
//            removedTriangles.insert(i); // Remove this triangle
//            ++numRemovedTriangles;
//        }
//    }
//    assert(numRemovedTriangles == int(removedTriangles.size()));
//
//    std::vector<int> newVertexIdx(vertices.size());
//    int numModifiedVertices = int(modifiedVertices.size());
//    newVertices.resize(
//        numModifiedVertices +
//        (vertices.size() - replacedVertex.size())
//    );
//    // Copy modified vertices to the beginning of newVertices array
//    for (size_t i = 0; i < modifiedVertices.size(); ++i) {
//        newVertices.at(i) = modifiedVertices.at(i);
//    }
//    int currentVertexIdx = numModifiedVertices;
//    for (size_t i = 0; i < vertices.size(); ++i) {
//        if (replacedVertex.count(i) == 0) {
//            // The i-th vertex is not replaced
//            newVertices.at(currentVertexIdx) = vertices.at(i);
//            newVertexIdx.at(i) = currentVertexIdx;
//            ++currentVertexIdx;
//        }
//        else {
//            // The i-th vertex was replaced
//            newVertexIdx.at(i) = replacedVertex[i];
//        }
//    }
//
//    newTriangles.resize(
//        triangles.size() - numRemovedTriangles
//    );
//    int currentTriangleIdx = 0;
//    for (size_t i = 0; i < triangles.size(); ++i) {
//        if (removedTriangles.count(i) != 0)
//            continue;
//        const Triangle& t = triangles.at(i);
//        newTriangles.at(currentTriangleIdx) =
//            Triangle(
//                newVertexIdx.at(t[0]),
//                newVertexIdx.at(t[1]),
//                newVertexIdx.at(t[2])
//            );
//        ++currentTriangleIdx;
//    }
//
//    // Check edges
//    int numManifoldEdges = 0;
//    int numBorderEdges = 0;
//    int numNonManifoldEdges = 0;
//
//    //qDebug() << "New vertices: " << newVertices.size() << endl;
//    //qDebug() << "New triangles: " << newTriangles.size() << endl;
//
//    std::map<Edge, int> edgeDegree;
//    for (size_t i = 0; i < newTriangles.size(); ++i) {
//        const Triangulation::Triangle& t = newTriangles.at(i);
//        for (int i = 0; i < 3; ++i) {
//            int v0 = t[i];
//            int v1;
//            if (i < 2)
//                v1 = t[i + 1];
//            else
//                v1 = t[0];
//            Edge te(v0, v1);
//            if (edgeDegree.count(te) == 0) {
//                edgeDegree[te] = 1;
//            }
//            else {
//                ++(edgeDegree[te]);
//            }
//        }
//    }
//
//    std::map<Edge, int>::const_iterator i = edgeDegree.cbegin();
//    while (i != edgeDegree.cend()) {
//        assert(i->second > 0);
//        if (i->second == 2) {
//            ++numManifoldEdges;
//        }
//        else if (i->second == 1) {
//            ++numBorderEdges;
//            //qDebug() << "Border edge: " << i->first.vertIdx[0]
//            //    << " " << i->first.vertIdx[1] << endl;
//        }
//        else {
//            ++numNonManifoldEdges;
//            //qDebug() << "Non-manifold edge: " << i->first.vertIdx[0]
//            //    << " " << i->first.vertIdx[1] << endl;
//            //qDebug() << "Edge degree: " << i->second << endl;
//        }
//        ++i;
//    }
//
//    vertices = newVertices;
//    triangles = newTriangles;
//}

double Triangulation::simplify(
    int& numRemovedEdges,
    int& numRemovedTriangles,
    double eps // = 0.01
) {
    defineTrianglesOfVertices();

    std::vector<Vertex> newVertices;
    std::vector<Triangle> newTriangles;

    std::map<int, int> replacedVertex;
    std::vector<Vertex> modifiedVertices;
    std::set<int> removedTriangles;
    std::set<Edge> removedEdges;
    std::set<int> frozenVertices;

    int numComplicatedEdges = 0;
    numRemovedEdges = 0;

    for (size_t i = 0; i < triangles.size(); ++i) {
        const Triangle& t = triangles.at(i);

        int vertexIdx0 = t.indices[0];
        int vertexIdx1 = t.indices[1];
        int vertexIdx2 = t.indices[2];

        Edge edge0(vertexIdx0, vertexIdx1);
        Edge edge1(vertexIdx1, vertexIdx2);
        Edge edge2(vertexIdx2, vertexIdx0);

        Vertex v0 = vertices.at(vertexIdx0);
        Vertex v1 = vertices.at(vertexIdx1);
        Vertex v2 = vertices.at(vertexIdx2);

        R3Point p0 = v0.point;
        R3Point p1 = v1.point;
        R3Point p2 = v2.point;

        double len[3];
        len[0] = p0.distance(p1);
        len[1] = p1.distance(p2);
        len[2] = p2.distance(p0);

        int numSmallEdges = 0;
        int smallEdgeIdx = (-1);
        if (
            len[0] <= eps &&
            replacedVertex.count(vertexIdx0) == 0 &&
            replacedVertex.count(vertexIdx1) == 0 &&

            frozenVertices.count(vertexIdx0) == 0 &&
            frozenVertices.count(vertexIdx1) == 0
        ) {
            smallEdgeIdx = 0;
            ++numSmallEdges;
        }
        if (
            len[1] <= eps &&
            replacedVertex.count(vertexIdx1) == 0 &&
            replacedVertex.count(vertexIdx2) == 0 &&

            frozenVertices.count(vertexIdx1) == 0 &&
            frozenVertices.count(vertexIdx2) == 0 &&
            (
                smallEdgeIdx < 0 ||
                len[1] < len[smallEdgeIdx]
            )
        ) {
            smallEdgeIdx = 1;
            ++numSmallEdges;
        }
        if (
            len[2] <= eps &&
            replacedVertex.count(vertexIdx2) == 0 &&
            replacedVertex.count(vertexIdx0) == 0 &&

            frozenVertices.count(vertexIdx2) == 0 &&
            frozenVertices.count(vertexIdx0) == 0 &&
            (
                smallEdgeIdx < 0 ||
                len[2] < len[smallEdgeIdx]
            )
        ) {
            smallEdgeIdx = 2;
            ++numSmallEdges;
        }
        if (numSmallEdges == 0)
            continue;

        //... if (numRemovedEdges >= TEST_MAX_EDGES)
        //...     continue;           //?????

        // Check whether we can remove this edge
        // preserving a topology (manifold edges)
        int vertIdx0 = t.indices[smallEdgeIdx];
        int vertIdx1 = t.indices[(smallEdgeIdx + 1)%3];
        const VertexStar& star0 = starsOfVertices.at(vertIdx0);
        const VertexStar& star1 = starsOfVertices.at(vertIdx1);

        int numCommonVertices = 0;

        VertexStar::const_iterator s1 = star1.cbegin();
        while (s1 != star1.cend()) {
            if (star0.count(*s1) != 0) {
                ++numCommonVertices;

                /*...
                int commonVertex = *s1;
                qDebug() << "    Common vertex: " << commonVertex;
                // Print a star of a common vertex
                qDebug() << "    Star of vertex " << commonVertex << ":";
                VertexStar::const_iterator cvs =
                    starsOfVertices.at(commonVertex).cbegin();
                while (cvs != starsOfVertices.at(commonVertex).cend()) {
                    qDebug() << "    " << *cvs;
                    ++cvs;
                }
                qDebug() << "    ---- End of Star";
                ...*/

            }
            ++s1;
        }
        if (numCommonVertices > 2) {

            /*...
            qDebug() << "Complicated edge: "
                << vertIdx0 << vertIdx1;
            qDebug() << "  Number of common vertices: " << numCommonVertices;

            //------------------------------------------------
            qDebug() << "  Star of vertex " << vertIdx0 << ":";
            s1 = star0.cbegin();
            while (s1 != star0.cend()) {
                qDebug() << *s1 << " ";
                ++s1;
            }

            qDebug() << "  Star of vertex " << vertIdx1 << ":";
            s1 = star1.cbegin();
            while (s1 != star1.cend()) {
                qDebug() << *s1 << " ";
                ++s1;
            }
            //------------------------------------------------
            ...*/

            ++numComplicatedEdges;
            continue;   // Do not remove this edge
        }

        Edge smallEdge(vertIdx0, vertIdx1);

        assert(removedEdges.count(smallEdge) == 0);
        assert(
            replacedVertex.count(vertIdx0) == 0 &&
            replacedVertex.count(vertIdx1) == 0
        );

        // Replace 2 vertices by the center of edge
        R3Point center =
            vertices.at(vertIdx0).point +
            (vertices.at(vertIdx1).point - vertices.at(vertIdx0).point)*0.5;
        R3Vector normal = vertices.at(vertIdx0).normal +
            vertices.at(vertIdx1).normal;
        normal.normalize();
        modifiedVertices.push_back(
            Vertex(center, normal)
        );
        int replacedVertexIdx = int(modifiedVertices.size() - 1);
        replacedVertex[vertIdx0] = replacedVertexIdx;
        replacedVertex[vertIdx1] = replacedVertexIdx;

        removedEdges.insert(smallEdge);
        ++numRemovedEdges;

        VertexStar::const_iterator s0 = star0.cbegin();
        while (s0 != star0.cend()) {
            int v = *s0;
            if (v != vertIdx0 && v != vertIdx1)
                frozenVertices.insert(v);
            ++s0;
        }
        s1 = star1.cbegin();
        while (s1 != star1.cend()) {
            int v = *s1;
            if (v != vertIdx0 && v != vertIdx1)
                frozenVertices.insert(v);
            ++s1;
        }

        //... removedTriangles.insert(i); // Remove this triangle
    } // end for

    //qDebug() << "numComplicatedEdges: " << numComplicatedEdges/2;
    //qDebug() << "Removed simple edges: " << removedEdges.size();

    numRemovedTriangles = 0;
    for (size_t i = 0; i < triangles.size(); ++i) {
        const Triangle& t = triangles.at(i);

        int vertexIdx0 = t.indices[0];
        int vertexIdx1 = t.indices[1];
        int vertexIdx2 = t.indices[2];

        Edge edge0(vertexIdx0, vertexIdx1);
        Edge edge1(vertexIdx1, vertexIdx2);
        Edge edge2(vertexIdx2, vertexIdx0);

        if (
            removedEdges.count(edge0) != 0 ||
            removedEdges.count(edge1) != 0 ||
            removedEdges.count(edge2) != 0
        ) {
            removedTriangles.insert(i); // Remove this triangle
            ++numRemovedTriangles;
        }
    }
    assert(numRemovedTriangles == int(removedTriangles.size()));

    std::vector<int> newVertexIdx(vertices.size());
    int numModifiedVertices = int(modifiedVertices.size());
    newVertices.resize(
        numModifiedVertices +
        (vertices.size() - replacedVertex.size())
    );
    // Copy modified vertices to the beginning of newVertices array
    for (size_t i = 0; i < modifiedVertices.size(); ++i) {
        newVertices.at(i) = modifiedVertices.at(i);
    }
    int currentVertexIdx = numModifiedVertices;
    for (size_t i = 0; i < vertices.size(); ++i) {
        if (replacedVertex.count(i) == 0) {
            // The i-th vertex is not replaced
            newVertices.at(currentVertexIdx) = vertices.at(i);
            newVertexIdx.at(i) = currentVertexIdx;
            ++currentVertexIdx;
        } else {
            // The i-th vertex was replaced
            newVertexIdx.at(i) = replacedVertex[i];
        }
    }

    newTriangles.resize(
        triangles.size() - numRemovedTriangles
    );
    int currentTriangleIdx = 0;
    for (size_t i = 0; i < triangles.size(); ++i) {
        if (removedTriangles.count(i) != 0)
            continue;
        const Triangle& t = triangles.at(i);
        newTriangles.at(currentTriangleIdx) =
            Triangle(
                newVertexIdx.at(t[0]),
                newVertexIdx.at(t[1]),
                newVertexIdx.at(t[2])
            );
        ++currentTriangleIdx;
    }

    // Check edges
    int numManifoldEdges = 0;
    int numBorderEdges = 0;
    int numNonManifoldEdges = 0;

    //qDebug() << "New vertices: " << newVertices.size();
    //qDebug() << "New triangles: " << newTriangles.size();

    std::map<Edge, int> edgeDegree;
    for (size_t i = 0; i < newTriangles.size(); ++i) {
        const Triangulation::Triangle& t = newTriangles.at(i);
        for (int i = 0; i < 3; ++i) {
            int v0 = t[i];
            int v1;
            if (i < 2)
                v1 = t[i + 1];
            else
                v1 = t[0];
            Edge te(v0, v1);
            if (edgeDegree.count(te) == 0) {
                edgeDegree[te] = 1;
            } else {
                ++(edgeDegree[te]);
            }
        }
    }

    std::map<Edge, int>::const_iterator i = edgeDegree.cbegin();
    while (i != edgeDegree.cend()) {
        assert(i->second > 0);
        if (i->second == 2) {
            ++numManifoldEdges;
        } else if (i->second == 1) {
            ++numBorderEdges;
            //qDebug() << "Border edge: " << i->first.vertIdx[0]
            //    << " " << i->first.vertIdx[1];
        } else {
            ++numNonManifoldEdges;
            //qDebug() << "Non-manifold edge: " << i->first.vertIdx[0]
            //    << " " << i->first.vertIdx[1];
            //qDebug() << "Edge degree: " << i->second;
        }
        ++i;
    }

    vertices = newVertices;
    triangles = newTriangles;

    invalidateAdjacentTriangles();
    invalidateTrianglesOfVertices();

    double minLen = 1e+30;
    for (size_t i = 0; i < triangles.size(); ++i) {
        const Triangle& t = triangles.at(i);

        int vertexIdx0 = t.indices[0];
        int vertexIdx1 = t.indices[1];
        int vertexIdx2 = t.indices[2];
        Vertex v0 = vertices.at(vertexIdx0);
        Vertex v1 = vertices.at(vertexIdx1);
        Vertex v2 = vertices.at(vertexIdx2);

        R3Point p0 = v0.point;
        R3Point p1 = v1.point;
        R3Point p2 = v2.point;

        double len = p0.distance(p1);
        if (len < minLen)
            minLen = len;
        len = p1.distance(p2);
        if (len < minLen)
            minLen = len;
        len = p2.distance(p0);
        if (len < minLen)
            minLen = len;
    }
    return minLen;
}

void Triangulation::cotangentLaplaceSmoothing(
    double lambda /* = 0.330 */  // May be negative for inflation/Taubin smooth
) {
    //... assert(adjacentTrianglesCalculated);
    if (!adjacentTrianglesCalculated)
        defineAdjacentTriangles();
    // Also starsOfVertices are calculated
    //... assert(starsOfVertices.size() == vertices.size());
    if (starsOfVertices.size() != vertices.size())
        defineTrianglesOfVertices();
    assert(ringsOfVertices.size() == vertices.size());

    std::vector<R3Point> newVertices(vertices.size());
    for (int i = 0; i < int(vertices.size()); ++i) {
        if (!ringsOfVertices.at(i).isBorderVertex()) {
            R3Vector laplaceShift = cotangentLaplace(i);
            newVertices[i] = vertices[i].point + laplaceShift * lambda;
        }
        else {
            newVertices[i] = vertices[i].point;
        }
    }

    // Copy shifted vertices back
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].point = newVertices[i];
    }
}

void Triangulation::uniformLaplaceSmoothing(
    double lambda /* = 0.330 */  // May be negative for inflation/Taubin smooth
) {
    //... assert(adjacentTrianglesCalculated);
    if (!adjacentTrianglesCalculated)
        defineAdjacentTriangles();
    // Also starsOfVertices are calculated
    //... assert(starsOfVertices.size() == vertices.size());
    if (starsOfVertices.size() != vertices.size())
        defineTrianglesOfVertices();
    assert(ringsOfVertices.size() == vertices.size());

    std::vector<R3Point> newVertices(vertices.size());
    for (int i = 0; i < int(vertices.size()); ++i) {
        if (!ringsOfVertices.at(i).isBorderVertex()) {
            R3Vector laplaceShift = uniformLaplace(i);
            newVertices[i] = vertices[i].point + laplaceShift * lambda;
        }
        else {
            newVertices[i] = vertices[i].point;
        }
    }

    // Copy shifted vertices back
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].point = newVertices[i];
    }
}


void Triangulation::taubinSmoothing(
    int iterations /* = 1 */,
    double lambda /* = 0.330 */,
    double mu /* = 0.331 */,
    bool useCotangentLaplace /* = false */
) {
    defineTrianglesOfVertices();
    assert(starsOfVertices.size() == vertices.size());

    if (useCotangentLaplace) {
        for (int i = 0; i < iterations; ++i) {
            cotangentLaplaceSmoothing(
                lambda
            );
            cotangentLaplaceSmoothing(
                -mu
            );
        }
    }
    else {
        for (int i = 0; i < iterations; ++i) {
            uniformLaplaceSmoothing(
                lambda
            );
            uniformLaplaceSmoothing(
                -mu
            );
        }
    }
}

R3Vector Triangulation::cotangentLaplace(int vertexIdx) const {
    if (!adjacentTrianglesCalculated)
        defineAdjacentTriangles();
    if (!trianglesOfVerticesCalculated)
        defineTrianglesOfVertices();
    assert(ringsOfVertices.size() == vertices.size());

    const VertexRing& vertexRing = ringsOfVertices.at(vertexIdx);
    if (vertexRing.isBorderVertex() || vertexRing.size() == 0)
        return R3Vector(0., 0., 0.);

    int n = int(vertexRing.size());
    R3Point t = vertices.at(vertexIdx).point;
    if (n == 1) {
        int idx = vertexRing.at(0);
        R3Point p = vertices.at(idx).point;
        return p - t;
    }
    double weightsSum = 0.;
    R3Point centroid(0., 0., 0.); // Weighted centroid

    // Loop for each edge of umbrella
    for (int i = 0; i < n; ++i) {
        int i0 = i - 1;
        if (i0 < 0)
            i0 = n - 1;
        int i1 = i + 1;
        if (i1 >= n)
            i1 = 0;
        int idx = vertexRing[i];
        int idx0 = vertexRing[i0];
        int idx1 = vertexRing[i1];

        R3Point p = vertices.at(idx).point;
        R3Point p0 = vertices.at(idx0).point;
        R3Point p1 = vertices.at(idx1).point;
        double cotanAlpha = R3Vector::cotan(
            p - p0, t - p0
        );
        double cotanBeta = R3Vector::cotan(
            p - p1, t - p1
        );
        double w = (cotanAlpha + cotanBeta) / 2.;
        weightsSum += w;
        centroid += p * w;
    }
    centroid *= (1. / weightsSum);
    return centroid - t;
}

R3Vector Triangulation::uniformLaplace(int vertexIdx) const {
    if (!adjacentTrianglesCalculated)
        defineAdjacentTriangles();
    if (!trianglesOfVerticesCalculated)
        defineTrianglesOfVertices();
    assert(ringsOfVertices.size() == vertices.size());

    const VertexRing& vertexRing = ringsOfVertices.at(vertexIdx);
    if (vertexRing.isBorderVertex() || vertexRing.size() == 0)
        return R3Vector(0., 0., 0.);

    int n = int(vertexRing.size());
    R3Point t = vertices.at(vertexIdx).point;
    if (n == 1) {
        int idx = vertexRing.at(0);
        R3Point p = vertices.at(idx).point;
        return p - t;
    }
    if (n == 2) {
        int idx = vertexRing.at(0);
        R3Point p = vertices.at(idx).point;
        idx = vertexRing.at(1);
        p += vertices.at(idx).point;
        p *= 0.5;
        return p - t;
    }

    R3Point centroid(0., 0., 0.); // Centroid of wireframe
    double weight = 0.;

    int idx0 = vertexRing[0];
    R3Point p0 = vertices.at(idx0).point;
    for (int i = 1; i <= n; ++i) {
        int j = i;
        if (j >= n)
            j = 0;
        int idx1 = vertexRing[j];
        R3Point p1 = vertices.at(idx1).point;
        double w = p0.distance(p1);
        centroid += (p0 + p1) * (w / 2.);
        weight += w;
        p0 = p1;
    }
    if (weight > R3_EPSILON)
        centroid *= (1. / weight);
    return centroid - t;
}
// TODO : Implement Skala mesh at surface voxels
void Triangulation::TriangulationOfTetrahedron(R3Graph::Tetrahedron& tetrahedron)
{
    std::vector<int> TriangleIndices;
    TriangleIndices.clear();
    std::set<int> EdgesWithPoints;

    for (int i = 0; i < 6; ++i)
    {
        double thrfun1 = tetrahedron.edges[i].A.second;
        double thrfun2 = tetrahedron.edges[i].B.second;
        if (thrfun1 * thrfun2 < 0.)
        {
            if (tetrahedron.edges[i].index == 0)
            {

                R3Point v = tetrahedron.edges[i].PointOnEdge();
                vertices.push_back(v);

                tetrahedron.edges[i].index = (int)vertices.size() - 1;
                TriangleIndices.push_back(tetrahedron.edges[i].index);
            }
            else
            {
                // vertex has already computed
                TriangleIndices.push_back(tetrahedron.edges[i].index);
            }
            EdgesWithPoints.insert(i);
        }
    }
    // Add triangles
    if (TriangleIndices.size() == 3)
    {
        Triangle t1(TriangleIndices[0],
            TriangleIndices[1],
            TriangleIndices[2]);
        t1.OutwardDirected(tetrahedron.Outward, vertices[TriangleIndices[0]].point,
            vertices[TriangleIndices[1]].point, vertices[TriangleIndices[2]].point);
        triangles.push_back(t1);
    }
    else if (TriangleIndices.size() == 4)
    {
        Triangle t1(TriangleIndices[0],
            TriangleIndices[1],
            TriangleIndices[2]);
        t1.OutwardDirected(tetrahedron.Outward, vertices[TriangleIndices[0]].point,
            vertices[TriangleIndices[1]].point, vertices[TriangleIndices[2]].point);
        triangles.push_back(t1);

        Triangle t2(TriangleIndices[1],
            TriangleIndices[2],
            TriangleIndices[3]);
        t2.OutwardDirected(tetrahedron.Outward, vertices[TriangleIndices[1]].point,
            vertices[TriangleIndices[2]].point, vertices[TriangleIndices[3]].point);
        triangles.push_back(t2);
    }
}
