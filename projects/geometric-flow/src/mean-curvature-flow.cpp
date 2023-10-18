// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"
#include "geometrycentral/numerical/linear_solvers.h"


/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {
    geometry->requireVertexIndices();
    return M + h * geometry->laplaceMatrix();
}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {

    geometry->requireVertexIndices();
    //solve system (M - hC)x_{k+1} = Mx_k
    //first we generate a |V| x 3 matrix, each row representing a vertex in space
    DenseMatrix<double> immersion(mesh->nVertices(), 3);
    auto i = 0u;
    for (Vertex v : mesh->vertices()) {
        auto coords = geometry->inputVertexPositions[v];
        immersion.row(i) << coords.x , coords.y , coords.z;
        i++;
    }
    auto M = geometry->massMatrix();
    //generate the right hand side
    Vector<double> x = M * immersion.col(0);
    Vector<double> y = M * immersion.col(1);
    Vector<double> z = M * immersion.col(2);

    //get the flow operator
    auto op = buildFlowOperator(M, h);
    
    //solve the system to get x_{k+1}
    x = solvePositiveDefinite(op, x);
    y = solvePositiveDefinite(op, y);
    z = solvePositiveDefinite(op, z);

    i = 0u;
    for (Vertex v : mesh->vertices()) {
        geometry->inputVertexPositions[v] = Vector3{x[i], y[i], z[i]}; 
        i++;
    }
}