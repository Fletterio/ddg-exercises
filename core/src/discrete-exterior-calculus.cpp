// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    SparseMatrix<double> h0(mesh.nVertices(), mesh.nVertices());
    std::vector<Eigen::Triplet<double>> entries;
    for (Vertex v : mesh.vertices()) {
        size_t idx = vertexIndices[v];
        entries.push_back(
            Eigen::Triplet<double>(static_cast<double>(idx), static_cast<double>(idx), barycentricDualArea(v)));
    }
    h0.setFromTriplets(entries.cbegin(), entries.cend());
    return h0;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    SparseMatrix<double> h1(mesh.nEdges(), mesh.nEdges());
    std::vector<Eigen::Triplet<double>> entries;
    for (Edge e : mesh.edges()) {
        double idx = static_cast<double>(edgeIndices[e]);
        entries.push_back(
            Eigen::Triplet<double>(idx, idx, edgeCotanWeight(e)));
    }
    h1.setFromTriplets(entries.cbegin(), entries.cend());
    return h1;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    SparseMatrix<double> h2(mesh.nFaces(), mesh.nFaces());
    std::vector<Eigen::Triplet<double>> entries;
    for (Face f : mesh.faces()) {
        double idx = static_cast<double>(faceIndices[f]);
        entries.push_back(
            Eigen::Triplet<double>(idx, idx, 1.0 / faceArea(f)));
    }
    h2.setFromTriplets(entries.cbegin(), entries.cend());
    return h2;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    SparseMatrix<double> signedAdj(mesh.nEdges(), mesh.nVertices());
    std::vector<Eigen::Triplet<double>> entries;

    for (Edge e : mesh.edges()) {
        double edgeIdx = static_cast<double>(edgeIndices[e]);
        entries.push_back(Eigen::Triplet<double>(edgeIdx, static_cast<double>(vertexIndices[e.firstVertex()]), -1.0));
        entries.push_back(Eigen::Triplet<double>(edgeIdx, static_cast<double>(vertexIndices[e.secondVertex()]), 1.0));
    }

    signedAdj.setFromTriplets(entries.cbegin(), entries.cend());
    return signedAdj;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    SparseMatrix<double> signedAdj(mesh.nFaces(), mesh.nEdges());
    std::vector<Eigen::Triplet<double>> entries;

    for (Face f : mesh.faces()) {
        double faceIdx = static_cast<double>(faceIndices[f]);
        Halfedge he = f.halfedge();
        for (auto i = 0u; i < 3; i++) {
            double orientation = he.tailVertex() == he.edge().firstVertex() ? 1.0 : -1.0;
            entries.push_back(Eigen::Triplet<double>(faceIdx, static_cast<double>(edgeIndices[he.edge()]), orientation));
            he = he.next();
        }
    }

    signedAdj.setFromTriplets(entries.cbegin(), entries.cend());
    return signedAdj;
}

} // namespace surface
} // namespace geometrycentral