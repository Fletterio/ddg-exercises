// Implement member functions for HodgeDecomposition class.
#include "hodge-decomposition.h"
#include "./solvers.h"
#include "geometrycentral/numerical/linear_solvers.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HodgeDecomposition::HodgeDecomposition(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: build DEC operators
    this->hodge1 = geometry->buildHodgeStar1Form();
    this->hodge2 = geometry->buildHodgeStar2Form();
    this->d0 = geometry->buildExteriorDerivative0Form();
    this->d1 = geometry->buildExteriorDerivative1Form();

    // TODO: Build operator inverses.
    // Hint: Use the sparseInverseDiagonal() in utils/src/solvers.cpp to invert sparse diagonal matrices.
    this->hodge1Inv = sparseInverseDiagonal(this->hodge1);
    this->hodge2Inv = sparseInverseDiagonal(this->hodge2);
    this->d0T = this->d0.transpose();
    this->d1T = this->d1.transpose();

    // TODO: Construct 0-form Laplace matrix.
    // Shift matrix by a small constant (1e-8) to make it positive definite.
    this->A = geometry->laplaceMatrix();

    // TODO: Construct 2-form matrix.
    this->B = this->d1 * this->hodge1Inv * this->d1T;
}

/*
 * Compute the 0-form potential α by solving the system 𝛿dα = 𝛿ω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The exact component dα of ω.
 */
Vector<double> HodgeDecomposition::computeExactComponent(const Vector<double>& omega) const {

    Vector<double> rhs = this->d0T * this->hodge1 * omega;
    auto laplaceMatrix = this->A;
    Vector<double> alpha = solvePositiveDefinite(laplaceMatrix, rhs);
    return this->d0 * alpha;
}

/*
 * Compute the 2-form potential β by solving the system d𝛿β = dω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The coexact component 𝛿β of ω.
 */
Vector<double> HodgeDecomposition::computeCoExactComponent(const Vector<double>& omega) const {

    SparseMatrix<double> twoFormLaplace = this->B;
    Vector<double> rhs = this->d1 * omega;
    Vector<double> betaTilde = solveSquare(twoFormLaplace, rhs);
    return this->hodge1Inv * this->d1T * betaTilde;
}

/*
 * Compute the harmonic component γ = ω - dα - 𝛿β of ω.
 *
 * Input: A primal 1-form <omega> on the edges of the input mesh, the exact component <dAlpha> of ω, and the coexact
 * component <deltaBeta> of ω.
 * Returns: The coexact component 𝛿β of ω.
 */
Vector<double> HodgeDecomposition::computeHarmonicComponent(const Vector<double>& omega, const Vector<double>& dAlpha,
                                                            const Vector<double>& deltaBeta) const {

    return omega - dAlpha - deltaBeta;
}