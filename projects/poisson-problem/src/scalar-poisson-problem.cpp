// Implement member functions for ScalarPoissonProblem class.
#include "scalar-poisson-problem.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ScalarPoissonProblem::ScalarPoissonProblem(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
    geometry->requireVertexIndices();
    geometry->requireFaceAreas();

    // TODO: Build member variables A (Laplace matrix), M (mass matrix), total area
    this->A = geometry->laplaceMatrix(); 
    this->M = geometry->massMatrix(); 
    this->totalArea = geometry->totalArea(); 

}

/*
 * Computes the solution of the poisson problem Ax = -M(rho - rhoBar), where A is the POSITIVE DEFINITE Laplace matrix
 * and M is the mass matrix.
 *
 * Input: <rho>, the density of vertices in the mesh.
 * Returns: The solution vector.
 */
Vector<double> ScalarPoissonProblem::solve(const Vector<double>& rho) const {
    SparseMatrix<double> matrix = this->A;
    double rhobar = (M * rho).sum() / totalArea;
    Vector<double> rhs = -M * (rho - Vector<double>::Constant(rho.size(), rhobar));
    return solvePositiveDefinite(matrix, rhs);
}