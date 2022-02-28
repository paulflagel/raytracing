#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <string.h>
#include <cmath>

#include "../Vector/Vector.h"
#include "Triangle.h"

BoundingBox::BoundingBox() {}

TriangleIndices::TriangleIndices(int vtxi, int vtxj, int vtxk, int ni, int nj, int nk, int uvi, int uvj, int uvk, int group, bool added)
{
    this->vtxi = vtxi;
    this->vtxj = vtxj;
    this->vtxk = vtxk;

    this->uvi = uvi;
    this->uvj = uvj;
    this->uvk = uvk;

    this->ni = ni;
    this->nj = nj;
    this->nk = nk;

    this->group = group;
}

TriangleMesh::TriangleMesh(const Vector &albedo, bool isMirror, bool isTransp, bool isLight, double refraction_index)
{
    this->albedo = albedo;
    this->isMirror = isMirror;
    this->isTransp = isTransp;
    this->isLight = isLight;
    this->refraction_index = refraction_index;
};

TriangleMesh::~TriangleMesh(){};

void TriangleMesh::readOBJ(const char *obj)
{

    char matfile[255];
    char grp[255];

    FILE *f;
    f = fopen(obj, "r");
    int curGroup = -1;
    while (!feof(f))
    {
        char line[255];
        if (!fgets(line, 255, f))
            break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's')
        {
            sscanf(line, "usemtl %[^\n]\n", grp);
            curGroup++;
        }

        if (line[0] == 'v' && line[1] == ' ')
        {
            Vector vec;

            Vector col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
            {
                col[0] = std::min(1., std::max(0., col[0]));
                col[1] = std::min(1., std::max(0., col[1]));
                col[2] = std::min(1., std::max(0., col[2]));

                vertices.push_back(vec);
                vertexcolors.push_back(col);
            }
            else
            {
                sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n')
        {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't')
        {
            Vector vec;
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f')
        {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;
            t.group = curGroup;

            char *consumedline = line + 1;
            int offset;

            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9)
            {
                if (i0 < 0)
                    t.vtxi = vertices.size() + i0;
                else
                    t.vtxi = i0 - 1;
                if (i1 < 0)
                    t.vtxj = vertices.size() + i1;
                else
                    t.vtxj = i1 - 1;
                if (i2 < 0)
                    t.vtxk = vertices.size() + i2;
                else
                    t.vtxk = i2 - 1;
                if (j0 < 0)
                    t.uvi = uvs.size() + j0;
                else
                    t.uvi = j0 - 1;
                if (j1 < 0)
                    t.uvj = uvs.size() + j1;
                else
                    t.uvj = j1 - 1;
                if (j2 < 0)
                    t.uvk = uvs.size() + j2;
                else
                    t.uvk = j2 - 1;
                if (k0 < 0)
                    t.ni = normals.size() + k0;
                else
                    t.ni = k0 - 1;
                if (k1 < 0)
                    t.nj = normals.size() + k1;
                else
                    t.nj = k1 - 1;
                if (k2 < 0)
                    t.nk = normals.size() + k2;
                else
                    t.nk = k2 - 1;
                indices.push_back(t);
            }
            else
            {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                if (nn == 6)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (k0 < 0)
                            t.ni = normals.size() + k0;
                        else
                            t.ni = k0 - 1;
                        if (k1 < 0)
                            t.nj = normals.size() + k1;
                        else
                            t.nj = k1 - 1;
                        if (k2 < 0)
                            t.nk = normals.size() + k2;
                        else
                            t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }

            consumedline = consumedline + offset;

            while (true)
            {
                if (consumedline[0] == '\n')
                    break;
                if (consumedline[0] == '\0')
                    break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.group = curGroup;
                if (nn == 3)
                {
                    if (i0 < 0)
                        t2.vtxi = vertices.size() + i0;
                    else
                        t2.vtxi = i0 - 1;
                    if (i2 < 0)
                        t2.vtxj = vertices.size() + i2;
                    else
                        t2.vtxj = i2 - 1;
                    if (i3 < 0)
                        t2.vtxk = vertices.size() + i3;
                    else
                        t2.vtxk = i3 - 1;
                    if (j0 < 0)
                        t2.uvi = uvs.size() + j0;
                    else
                        t2.uvi = j0 - 1;
                    if (j2 < 0)
                        t2.uvj = uvs.size() + j2;
                    else
                        t2.uvj = j2 - 1;
                    if (j3 < 0)
                        t2.uvk = uvs.size() + j3;
                    else
                        t2.uvk = j3 - 1;
                    if (k0 < 0)
                        t2.ni = normals.size() + k0;
                    else
                        t2.ni = k0 - 1;
                    if (k2 < 0)
                        t2.nj = normals.size() + k2;
                    else
                        t2.nj = k2 - 1;
                    if (k3 < 0)
                        t2.nk = normals.size() + k3;
                    else
                        t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (k0 < 0)
                                t2.ni = normals.size() + k0;
                            else
                                t2.ni = k0 - 1;
                            if (k2 < 0)
                                t2.nj = normals.size() + k2;
                            else
                                t2.nj = k2 - 1;
                            if (k3 < 0)
                                t2.nk = normals.size() + k3;
                            else
                                t2.nk = k3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(f);
};

void TriangleMesh::get_bbox()
{
    double xMin, xMax, yMin, yMax, zMin, zMax;
    for (int k = 0; k < vertices.size(); k++)
    {
        xMin = std::min(xMin, vertices[k][0]);
        xMax = std::max(xMax, vertices[k][0]);
        yMin = std::min(yMin, vertices[k][1]);
        yMax = std::max(yMax, vertices[k][1]);
        zMin = std::min(zMin, vertices[k][2]);
        zMax = std::max(zMax, vertices[k][2]);
    }
    this->bbox.min = Vector(xMin, yMin, zMin);
    this->bbox.max = Vector(xMax, yMax, zMax);
}

bool TriangleMesh::intersect(const Ray &r, Vector &P, Vector &N, double &t) const
{
    // On teste d'abord si on intersecte la bounding box du mesh
    if (!intersect_bbox(r))
    {
        return false;
    }

    t = 1E99; // distance au triangle le plus proche
    bool triangle_intersection = false;

    // On parcourt tous les triangles du mesh
    for (int k = 0; k < indices.size(); k++)
    {
        double t_intersect;
        Vector Ptriangle, Ntriangle;
        if (intersect_with_triangle(indices[k], r, Ptriangle, Ntriangle, t_intersect))
        {
            triangle_intersection = true;
            if (t_intersect < t)
            {
                t = t_intersect;
                P = Ptriangle;
                N = Ntriangle;
            }
        }
    }
    return triangle_intersection;
}

bool TriangleMesh::intersect_bbox(const Ray &r) const
{

    double tx1 = (bbox.min[0] - r.C[0]) / r.u[0];
    double tx2 = (bbox.max[0] - r.C[0]) / r.u[0];
    double txmin = std::min(tx1, tx2);
    double txmax = std::max(tx1, tx2);

    double ty1 = (bbox.min[1] - r.C[1]) / r.u[1];
    double ty2 = (bbox.max[1] - r.C[1]) / r.u[1];
    double tymin = std::min(ty1, ty2);
    double tymax = std::max(ty1, ty2);

    double tz1 = (bbox.min[2] - r.C[2]) / r.u[2];
    double tz2 = (bbox.max[2] - r.C[2]) / r.u[2];
    double tzmin = std::min(tz1, tz2);
    double tzmax = std::max(tz1, tz2);

    double tmin_max = std::min(txmax, std::min(tymax, tzmax));
    double tmax_min = std::max(txmin, std::max(tymin, tzmin));

    if (tmax_min <= tmin_max) // est ce que le max des min est inférieur au min des max
    {
        return true;
    }
    return false;
}

bool TriangleMesh::intersect_with_triangle(const TriangleIndices &triangle, const Ray &ray, Vector &P, Vector &N, double &t) const
{
    // Points du triangle
    Vector vertice_i = vertices[triangle.vtxi];
    Vector vertice_j = vertices[triangle.vtxj];
    Vector vertice_k = vertices[triangle.vtxk];

    Vector e_1 = vertice_j - vertice_i;
    Vector e_2 = vertice_k - vertice_i;

    N = cross(e_1, e_2);

    Vector i_to_origin = (ray.C - vertice_i);
    Vector i_to_origin_cross_u = cross(i_to_origin, ray.u);
    double u_dot_n = dot(ray.u, N);

    double beta = -dot(e_2, i_to_origin_cross_u) / u_dot_n;
    double gamma = +dot(e_1, i_to_origin_cross_u) / u_dot_n;
    double alpha = 1. - beta - gamma;
    t = -dot(i_to_origin, N) / u_dot_n;
    if (!(0 <= alpha && alpha <= 1 && 0 <= beta && beta <= 1 && 0 <= gamma && gamma <= 1 && t >= 0))
    {
        return false;
    }
    P = vertice_i + e_1 * beta + e_2 * gamma;
    N.normalize();
    return true;
}