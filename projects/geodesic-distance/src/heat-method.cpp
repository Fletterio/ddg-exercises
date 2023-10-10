// Implement member functions HeatMethod class.
#include "heat-method.h"
#include "geometrycentral/numerical/linear_solvers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;
    geometry->requireVertexIndices();
    geometry->requireVertexPositions();
    geometry->requireFaceAreas();
    geometry->requireFaceNormals();
    geometry->requireHalfedgeCotanWeights();

    // TODO: Build Laplace and flow matrices.
    // Note: core/geometry.cpp has meanEdgeLength() function
    this->A = geometry->laplaceMatrix(); 
    double h = geo->meanEdgeLength();
    double t = h * h;
    this->F = geometry->massMatrix() + t * this->A;
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {

    auto vectorField = FaceData<Vector3>(*mesh, {0, 0, 0});

    for (const Face& f : mesh->faces()) {
        Vector3 nabla_u = Vector3::zero();
        for (const Halfedge& he : f.adjacentHalfedges()) {
            auto orientedEdge = - geometry->vertexPositions[he.tipVertex()] + geometry->vertexPositions[he.tailVertex()];
            nabla_u += u[geometry->vertexIndices[he.next().tipVertex()]] * cross(geometry->faceNormals[f], orientedEdge);
        }
        vectorField[f] = nabla_u.normalize();

    }
    return vectorField; 
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {

    Vector<double> integratedDivergence = Vector<double>::Zero(mesh->nVertices());
    for (const Vertex& v : mesh->vertices()) {
        double intDiv = 0;
        for (const Face& f : v.adjacentFaces()) {
            auto e1 = f.halfedge();
            Halfedge e2;
            if (e1.tipVertex() != v && e1.tailVertex() != v) {
                e1 = e1.next();
            }
            if (e1.tipVertex() == v) {
                e1 = e1.next();
            }
            e2 = e1.next().next();
            //e1 is the right he in cc order, e2 is inverted in cc order
            auto e1Or = geometry->vertexPositions[e1.tipVertex()] - geometry->vertexPositions[v];
            auto e2Or = geometry->vertexPositions[e2.tailVertex()] - geometry->vertexPositions[v];
            intDiv += geometry->halfedgeCotanWeights[e1] * dot(e1Or, X[f]) +
                      geometry->halfedgeCotanWeights[e2] * dot(e2Or, X[f]);
        }
        integratedDivergence[geometry->vertexIndices[v]] = intDiv;
    }
    return integratedDivergence; 
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {

    SparseMatrix<double> flow = this->F;
    auto u = solvePositiveDefinite(flow, delta);
    auto grad_u = computeVectorField(u);
    Vector<double> divergence = -1.0 * computeDivergence(grad_u);

    SparseMatrix<double> laplace = this->A;
    Vector<double> phi = solvePositiveDefinite(laplace, divergence);

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}