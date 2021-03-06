#pragma once

#include "Object.hpp"

#include <cstring>

static bool insideTriangle(float x, float y, Vector3f* _v)
{
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    Vector3f p(x, y, 1);
    Vector3f a(_v[0].x, _v[0].y, 1);
    Vector3f b(_v[1].x, _v[1].y, 1);
    Vector3f c(_v[2].x, _v[2].y, 1);

    auto pa = p - a;
    auto pb = p - b;
    auto pc = p - c;

    auto z1 = crossProduct(pa, pb).z;
    auto z2 = crossProduct(pb, pc).z;
    auto z3 = crossProduct(pc, pa).z;

    return z1 * z2 > 0 && z2 * z3 > 0;
}

bool rayTriangleIntersect(const Vector3f& v0, const Vector3f& v1, const Vector3f& v2, const Vector3f& orig,
                          const Vector3f& dir, float& tnear, float& u, float& v)
{
    // TODO: Implement this function that tests whether the triangle
    // that's specified bt v0, v1 and v2 intersects with the ray (whose
    // origin is *orig* and direction is *dir*)
    // Also don't forget to update tnear, u and v.
    
    /*
    Vector3f e1 = v1 - v0;
    Vector3f e2 = v2 - v0;
    Vector3f s0 = orig - v0;
    Vector3f s1 = crossProduct(dir, e2);
    Vector3f s2 = crossProduct(s0, e1);

    Vector3f s = Vector3f(dotProduct(s2, e2), dotProduct(s1, s0), dotProduct(s2, dir)) / dotProduct(s1, e1);
    tnear = s.x;
    u = s.y;
    v = s.z;

    if (tnear >= 0 && u >= 0 && v >= 0 && (u + v) <= 1)
    {
        return true;
    }

    return false;
    */
    
    Vector3f E1 = v1 - v0;
    Vector3f E2 = v2 - v0;
    Vector3f S = orig - v0;
    Vector3f S1 = crossProduct(dir, E2);
    Vector3f S2 = crossProduct(S, E1);

    tnear = dotProduct(S2, E2)  / dotProduct(S1, E1);
    u     = dotProduct(S1, S)   / dotProduct(S1, E1);
    v =     dotProduct(S2, dir) / dotProduct(S1, E1);
    if (tnear >= 0 && u >= 0 && v >= 0 && (u + v) <= 1)
    {
        return true;
    }
    return false;
    
    //Vector3f interPoint = orig + tnear * dir;
    //u = interPoint.x; v = interPoint.y;
    //Vector3f tri[3]; tri[0] = v0; tri[1] = v1; tri[2] = v2;

    //return insideTriangle(u, v, tri);
}

class MeshTriangle : public Object
{
public:
    MeshTriangle(const Vector3f* verts, const uint32_t* vertsIndex, const uint32_t& numTris, const Vector2f* st)
    {
        uint32_t maxIndex = 0;
        for (uint32_t i = 0; i < numTris * 3; ++i)
            if (vertsIndex[i] > maxIndex)
                maxIndex = vertsIndex[i];
        maxIndex += 1;
        vertices = std::unique_ptr<Vector3f[]>(new Vector3f[maxIndex]);
        memcpy(vertices.get(), verts, sizeof(Vector3f) * maxIndex);
        vertexIndex = std::unique_ptr<uint32_t[]>(new uint32_t[numTris * 3]);
        memcpy(vertexIndex.get(), vertsIndex, sizeof(uint32_t) * numTris * 3);
        numTriangles = numTris;
        stCoordinates = std::unique_ptr<Vector2f[]>(new Vector2f[maxIndex]);
        memcpy(stCoordinates.get(), st, sizeof(Vector2f) * maxIndex);
    }

    bool intersect(const Vector3f& orig, const Vector3f& dir, float& tnear, uint32_t& index,
                   Vector2f& uv) const override
    {
        bool intersect = false;
        for (uint32_t k = 0; k < numTriangles; ++k)
        {
            const Vector3f& v0 = vertices[vertexIndex[k * 3]];
            const Vector3f& v1 = vertices[vertexIndex[k * 3 + 1]];
            const Vector3f& v2 = vertices[vertexIndex[k * 3 + 2]];
            float t, u, v;
            if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tnear)
            {
                tnear = t;
                uv.x = u;
                uv.y = v;
                index = k;
                intersect |= true;
            }
        }

        return intersect;
    }

    void getSurfaceProperties(const Vector3f&, const Vector3f&, const uint32_t& index, const Vector2f& uv, Vector3f& N,
                              Vector2f& st) const override
    {
        const Vector3f& v0 = vertices[vertexIndex[index * 3]];
        const Vector3f& v1 = vertices[vertexIndex[index * 3 + 1]];
        const Vector3f& v2 = vertices[vertexIndex[index * 3 + 2]];
        Vector3f e0 = normalize(v1 - v0);
        Vector3f e1 = normalize(v2 - v1);
        N = normalize(crossProduct(e0, e1));
        const Vector2f& st0 = stCoordinates[vertexIndex[index * 3]];
        const Vector2f& st1 = stCoordinates[vertexIndex[index * 3 + 1]];
        const Vector2f& st2 = stCoordinates[vertexIndex[index * 3 + 2]];
        st = st0 * (1 - uv.x - uv.y) + st1 * uv.x + st2 * uv.y;
    }

    Vector3f evalDiffuseColor(const Vector2f& st) const override
    {
        float scale = 5;
        float pattern = (fmodf(st.x * scale, 1) > 0.5) ^ (fmodf(st.y * scale, 1) > 0.5);
        return lerp(Vector3f(0.815, 0.235, 0.031), Vector3f(0.937, 0.937, 0.231), pattern);
    }

    std::unique_ptr<Vector3f[]> vertices;
    uint32_t numTriangles;
    std::unique_ptr<uint32_t[]> vertexIndex;
    std::unique_ptr<Vector2f[]> stCoordinates;
};
