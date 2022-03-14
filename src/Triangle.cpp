#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <string.h>
#include <cmath>
#include <list>

#include "Vector.h"
#include "Triangle.h"

BoundingBox::BoundingBox(Vector min, Vector max)
{
    this->min = min;
    this->max = max;
}

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
    Vector bbmin(vertices[0][0], vertices[0][1], vertices[0][2]);
    Vector bbmax(vertices[0][0], vertices[0][1], vertices[0][2]);
    for (int k = 0; k < vertices.size(); k++)
    {
        for (int dim = 0; dim < 3; dim++)
        {
            bbmin[dim] = std::min(bbmin[dim], vertices[k][dim]);
            bbmax[dim] = std::max(bbmax[dim], vertices[k][dim]);
        }
    }
    this->bbox.min = bbmin;
    this->bbox.max = bbmax;
}

BoundingBox TriangleMesh::get_bbox(int lower, int upper)
{
    Vector bbmin(vertices[indices[lower].vtxi][0], vertices[indices[lower].vtxi][1], vertices[indices[lower].vtxi][2]);
    Vector bbmax(vertices[indices[lower].vtxi][0], vertices[indices[lower].vtxi][1], vertices[indices[lower].vtxi][2]);
    TriangleIndices triangle;
    for (int k = lower; k < upper; k++)
    {
        triangle = indices[k];
        for (int dim = 0; dim < 3; dim++)
        {
            bbmin[dim] = std::min(bbmin[dim], std::min(vertices[triangle.vtxi][dim], std::min(vertices[triangle.vtxj][dim], vertices[triangle.vtxk][dim])));
            bbmax[dim] = std::max(bbmax[dim], std::max(vertices[triangle.vtxi][dim], std::max(vertices[triangle.vtxj][dim], vertices[triangle.vtxk][dim])));
        }
    }
    return BoundingBox(bbmin, bbmax);
}

bool TriangleMesh::intersect(const Ray &r, Vector &P, Vector &N, double &t) const
{
    t = 1E99; // distance au triangle le plus proche
    bool triangle_intersection = false;

    if (use_bvh)
    {
        if (!(bvh.bbox.intersect(r)))
            return false;

        std::list<const BVH *> nodes;
        nodes.push_back(&bvh);

        while (!nodes.empty())
        {
            const BVH *currentBVH = nodes.front();
            nodes.pop_front();

            if (currentBVH->left_child)
            {
                if (currentBVH->left_child->bbox.intersect(r))
                    nodes.push_front(currentBVH->left_child);
                if (currentBVH->right_child->bbox.intersect(r))
                    nodes.push_front(currentBVH->right_child);
            }

            else
            {
                for (int k = currentBVH->lower_index; k < currentBVH->upper_index; k++)
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
            }
        }
        return triangle_intersection;
    }

    else
    {
        if (!bbox.intersect(r))
        {
            // On teste d'abord si on intersecte la bounding box du mesh
            return false;
        }

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
}

bool BoundingBox::intersect(const Ray &r) const
{
    double tx1 = (min[0] - r.C[0]) / r.u[0];
    double tx2 = (max[0] - r.C[0]) / r.u[0];
    double txmin = std::min(tx1, tx2);
    double txmax = std::max(tx1, tx2);

    double ty1 = (min[1] - r.C[1]) / r.u[1];
    double ty2 = (max[1] - r.C[1]) / r.u[1];
    double tymin = std::min(ty1, ty2);
    double tymax = std::max(ty1, ty2);

    double tz1 = (min[2] - r.C[2]) / r.u[2];
    double tz2 = (max[2] - r.C[2]) / r.u[2];
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

    if (normals_interpolation)
    {
        // Interpolation des normales
        Vector ni = normals[triangle.ni];
        Vector nj = normals[triangle.nj];
        Vector nk = normals[triangle.nk];
        N = (ni * alpha) + (ni * beta) + (nk * gamma);
    }
    // N = -N;
    N.normalize();
    return true;
}

void TriangleMesh::swapAxis(int axis1, int axis2)
{
    for (int k = 0; k < normals.size(); k++)
    {
        std::swap(normals[k][axis1], normals[k][axis2]);
    }

    for (int k = 0; k < vertices.size(); k++)
    {
        std::swap(vertices[k][axis1], vertices[k][axis2]);
    }
}

void TriangleMesh::invertNormals()
{
    for (int k = 0; k < normals.size(); k++)
    {
        normals[k] = (-1) * normals[k];
    }
}

void TriangleMesh::scale(float factor, Vector offset)
{
    for (int k = 0; k < vertices.size(); k++)
    {
        vertices[k] = factor * vertices[k] + offset;
    }
};

void TriangleMesh::rotate(int axis, double angle)
{
    angle = angle * M_PI / 180;
    double cosinus = cos(angle);
    double sinus = sin(angle);
    int dim1 = (axis + 1) % 3;
    int dim2 = (axis + 2) % 3;

    for (int k = 0; k < vertices.size(); k++)
    {
        double tmp1 = vertices[k][dim1];
        double tmp2 = vertices[k][dim2];
        vertices[k][dim1] = tmp1 * cosinus + tmp2 * sinus;
        vertices[k][dim2] = -tmp1 * sinus + tmp2 * cosinus;
    }
    for (int k = 0; k < normals.size(); k++)
    {
        double tmp1 = normals[k][dim1];
        double tmp2 = normals[k][dim2];
        normals[k][dim1] = tmp1 * cosinus + tmp2 * sinus;
        normals[k][dim2] = -tmp1 * sinus + tmp2 * cosinus;
    }
}

void TriangleMesh::build_BVH(BVH *n, int lower, int upper)
{
    n->bbox = get_bbox(lower, upper);
    n->lower_index = lower;
    n->upper_index = upper;
    n->left_child = NULL;
    n->right_child = NULL;

    Vector diag = n->bbox.max - n->bbox.min; // Diagonale de la bbox
    int dim = -1;                            // Axe selon lequel on sépare les bbox récursivement

    if (diag[0] >= std::max(diag[1], diag[2]))
        dim = 0;
    else if (diag[1] >= std::max(diag[0], diag[2]))
        dim = 1;
    else if (diag[2] >= std::max(diag[0], diag[1]))
        dim = 2;

    double split_value = (n->bbox.min[dim] + n->bbox.max[dim]) / 2;

    int pivot = lower;

    for (int k = lower; k < upper; k++)
    {
        TriangleIndices triangle = indices[k];
        double centre_gravite = (vertices[triangle.vtxi][dim] + vertices[triangle.vtxj][dim] + vertices[triangle.vtxk][dim]) / 3;

        if (centre_gravite < split_value)
        {
            std::swap(indices[k], indices[pivot]); // quick sort
            pivot++;
        }
    }

    if ((pivot <= lower) || (pivot >= upper - 1) || (upper - lower <= 5))
        return;

    n->left_child = new BVH;
    n->right_child = new BVH;
    build_BVH(n->left_child, lower, pivot + 1);
    build_BVH(n->right_child, pivot + 1, upper);
}

void TriangleMesh::init_BVH()
{
    build_BVH(&bvh, 0, indices.size());
}