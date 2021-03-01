#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <clocale>
#include <cassert>
#include <set>
#include "Triangulation.h"

using namespace R3Graph;

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
}

void Triangulation::orientate() {
    for (size_t i = 0; i < triangles.size(); ++i) {
        Triangle t = triangles.at(i);
        Vertex v0 = vertices.at(t.indices[0]);
        Vertex v1 = vertices.at(t.indices[1]);
        Vertex v2 = vertices.at(t.indices[2]);
        R3Point p0 = v0.point;
        R3Point p1 = v1.point;
        R3Point p2 = v2.point;
        //R3Vector meanN = (v0.normal + v1.normal + v2.normal)*(1./3.);
        t.Normal = (p1 - p0).vectorProduct(p2 - p0);
        t.Normal.normalize();
        //R3Vector Outward = p0 - imageCenter;
        //if (Outward.scalarProduct(t.Normal) < 0.){}
        //    t.RightHand();
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

