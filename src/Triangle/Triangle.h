#ifndef TRIANGLE
#define TRIANGLE

#include "../Vector/Vector.h"
#include "../Ray/Ray.h"
#include "../Object/Object.h"

class BoundingBox
{
public:
    BoundingBox();

    Vector min, max;
};

class TriangleIndices
{
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false);
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;    // indices within the uv coordinates array
    int ni, nj, nk;       // indices within the normals array
    int group;            // face group
};

class TriangleMesh : public Object
{
public:
    TriangleMesh(const Vector &albedo, bool isMirror = false, bool isTransp = false, bool isLight = false, double refraction_index = 1.4);
    ~TriangleMesh();

    void readOBJ(const char *obj);

    virtual bool intersect(const Ray &r, Vector &P, Vector &N, double &t) const;
    bool intersect_with_triangle(const TriangleIndices &triangle, const Ray &ray, Vector &P, Vector &N, double &t) const;

    std::vector<TriangleIndices> indices; // Indices des sommets qui constituent chaque face
    std::vector<Vector> vertices;         // Coordonnées des sommets
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;

    BoundingBox bbox;
    bool intersect_bbox(const Ray &r) const;
    void get_bbox();
};

#endif