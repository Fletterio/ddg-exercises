// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

using namespace std::complex_literals;

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {
    size_t v1 = he.tipVertex().getIndex();
    size_t v2 = he.tailVertex().getIndex();

    size_t v3 = he.next().tipVertex().getIndex();

    auto u = vertexPositions[v1] - vertexPositions[v3];
    auto v = vertexPositions[v2] - vertexPositions[v3];

    double denominator = cross(u, v).norm();

    return dot(u, v) / denominator;

    //Default method
    //return 2.0 * halfedgeCotanWeight(he);      
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */

//In case we have not yet covered it in class, 
//the barycentric dual area associated with a vertex is equal to one-third the area of all triangles ijk touching i.
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {
    double area = 0.0;
    for (Face f : v.adjacentFaces()) {
        if (!f.isBoundaryLoop()) area += faceArea(f);
    }
    return area / 3;
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {

    auto o = c.vertex();
    auto p = c.halfedge().tipVertex();
    auto q = c.halfedge().next().tipVertex();

    auto v = vertexPositions[p] - vertexPositions[o];
    auto w = vertexPositions[q] - vertexPositions[o];


    return clamp(acos(dot(v, w) / (v.norm() * w.norm())), 0.0, PI);
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {

    auto n1 = faceNormals[he.face()];
    auto n2 = faceNormals[he.twin().face()];
    auto v = vertexPositions[he.tipVertex()] - vertexPositions[he.tailVertex()];

    return atan2(dot(v.normalize(), cross(n1, n2)), dot(n1, n2));
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {

    auto n = Vector3::zero();

    for (const auto& f : v.adjacentFaces()) {
        n += faceNormals[f.getIndex()];
    }

    return n.normalize();
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {

    auto n = Vector3::zero();
    for (const auto& c : v.adjacentCorners()) {
        n += angle(c) * faceNormals[c.face()];
    }
    return n.normalize();

}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    auto first = v.halfedge();
    auto left = first;
    auto right = left.twin().next();
    auto n = Vector3::zero();
    do {
        auto v = vertexPositions[right.tipVertex()] - vertexPositions[right.tailVertex()];
        auto w = vertexPositions[left.tipVertex()] - vertexPositions[left.tailVertex()];
        n += cross(v, w) / (v.norm2() * w.norm2());
        left = right;
        right = right.twin().next();
    } while (left != first);
    return n.normalize();
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    auto n = Vector3::zero();
    for (auto const& f : v.adjacentFaces()) {
        n += faceAreas[f] * faceNormals[f];
    }
    return n.normalize();
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    auto n = Vector3::zero();
    for (auto const& he : v.outgoingHalfedges()) {
        auto u = vertexPositions[he.tipVertex()] - vertexPositions[he.tailVertex()];
        n += dihedralAngle(he) * u.normalize();
    }
    return n.normalize();
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    auto n = Vector3::zero();
    for (auto const& he : v.outgoingHalfedges()) {
        auto u = vertexPositions[he.tipVertex()] - vertexPositions[he.tailVertex()];
        n += edgeCotanWeights[he.edge()] * u;
    }
    return n.normalize();
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    double angleSum = 0.0;
    for (auto const& c : v.adjacentCorners()) {
        angleSum += angle(c);
    }
    return 2 * PI - angleSum;
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {

    double defectSum = 0.0;
    for (auto const& v : mesh.vertices()) {
        defectSum += angleDefect(v);
    }
    return defectSum;
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {
    
    double sum = 0.0;
    for (auto const& he : v.outgoingHalfedges()) {
        auto u = vertexPositions[he.tipVertex()] - vertexPositions[he.tailVertex()]; 
        sum += dihedralAngle(he) * u.norm();
    }
    return sum / 2;
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {

    auto first = v.halfedge();
    auto left = first;
    auto right = left.twin().next();
    double area = 0.0;
    do {
        auto v = vertexPositions[right.tipVertex()] - vertexPositions[right.tailVertex()];
        auto w = vertexPositions[left.tipVertex()] - vertexPositions[left.tailVertex()];
        area += halfedgeCotanWeight(left.twin()) * w.norm2() + halfedgeCotanWeight(right) * v.norm2();
        left = right;
        right = right.twin().next();
    } while (left != first);
    return area / 4;
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    auto K = angleDefect(v) / circumcentricDualArea(v);
    auto H = scalarMeanCurvature(v) / circumcentricDualArea(v);
    auto k2 = H + sqrt(H * H - K);
    auto k1 = H - sqrt(H * H - K);
    if (k1 > k2) std::swap(k1, k2);
    return std::make_pair(k1, k2); // placeholder
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    SparseMatrix<double> laplace(mesh.nVertices(), mesh.nVertices());
    std::vector<Eigen::Triplet<double>> entries;

    // row i represents the image of vertex i
    for (Vertex v : mesh.vertices()) {
        double cum = 0.0;
        for (Halfedge he : v.outgoingHalfedges()) {
            cum += edgeCotanWeight(he.edge());
            entries.push_back(Eigen::Triplet<double>(static_cast<double>(vertexIndices[v]),
                                                     static_cast<double>(vertexIndices[he.tipVertex()]),
                                                     - edgeCotanWeight(he.edge())));
        }
        entries.push_back(
            Eigen::Triplet<double>(static_cast<double>(vertexIndices[v]), static_cast<double>(vertexIndices[v]), cum + 1e-8));
    }
    laplace.setFromTriplets(entries.cbegin(), entries.cend());

    return laplace;
}
/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    SparseMatrix<double> mass(mesh.nVertices(), mesh.nVertices());
    std::vector<Eigen::Triplet<double>> entries;
    for (Vertex v : mesh.vertices()) {
        entries.push_back(
            Eigen::Triplet<double>(static_cast<double>(vertexIndices[v]), static_cast<double>(vertexIndices[v]), barycentricDualArea(v)));
    }
    mass.setFromTriplets(entries.cbegin(), entries.cend());
    return mass;
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    SparseMatrix<std::complex<double>> laplace(mesh.nVertices(), mesh.nVertices());
    std::vector<Eigen::Triplet<std::complex<double>>> entries;

    // row i represents the image of vertex i
    for (Vertex v : mesh.vertices()) {
        std::complex<double> cum(0.0,0.0);
        for (Halfedge he : v.outgoingHalfedges()) {
            cum += edgeCotanWeight(he.edge());
            entries.push_back(Eigen::Triplet<std::complex<double>>(
                static_cast<double>(vertexIndices[v]),
                static_cast<double>(vertexIndices[he.tipVertex()]), -edgeCotanWeight(he.edge())));
        }
        entries.push_back(Eigen::Triplet<std::complex<double>>(static_cast<double>(vertexIndices[v]),
                                                               static_cast<double>(vertexIndices[v]),
                                                               cum + 1e-8));
    }
    laplace.setFromTriplets(entries.cbegin(), entries.cend());

    return laplace;
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral