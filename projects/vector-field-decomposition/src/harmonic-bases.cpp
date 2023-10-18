// Implement member functions for HarmonicBases class.
#include "harmonic-bases.h"
#include "solvers.h"

#include "geometrycentral/numerical/linear_solvers.h"
//-------------------------------------------- Here the actual file begins ----------------------------------------------------

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HarmonicBases::HarmonicBases(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build a closed, but not exact, primal 1-form ω.
 *
 * Input: A std::vector of Halfedges representing a homology generator of the mesh.
 * Returns: A vector representing a closed primal 1-form.
 */
Vector<double> HarmonicBases::buildClosedPrimalOneForm(const std::vector<Halfedge>& generator) const {

    Vector<double> oneForm = Vector<double>::Zero(mesh->nEdges());
    int sign = 1;
    Halfedge previous, current;
    previous = generator[generator.size() - 1];
    for (auto& he : generator) {
        current = he;
        if (previous.tipVertex() == current.tailVertex()) {
            sign = -sign;
        }
        if (he.edge().halfedge() == he) {
            oneForm[he.edge().getIndex()] = sign;
        } else {
            oneForm[he.edge().getIndex()] = -sign;
        }
        previous = he;
    }
    return oneForm;
}

/*
 * Compute the harmonic bases [γ1, γ2 ... γn] of the input mesh.
 *
 * Input: A std::vector of homology generators of the mesh (which are in turn represented as std::vectors of halfedges),
 * and a HodgeDecomposition object. Returns:
 */
std::vector<Vector<double>> HarmonicBases::compute(const std::vector<std::vector<Halfedge>>& generators,
                                                   const HodgeDecomposition& hodgeDecomposition) const {
    std::vector<Vector<double>> gammas;
    SparseMatrix<double> hodge0, hodge0Inverse(mesh->nVertices(), mesh->nVertices()), diff0, hodge1, codifferential;
    hodge0 = geometry->buildHodgeStar0Form();
    hodge0Inverse = sparseInverseDiagonal(hodge0);
    diff0 = geometry->buildExteriorDerivative0Form();
    hodge1 = geometry->buildHodgeStar1Form();
    codifferential = hodge0Inverse * diff0.transpose() * hodge1;

    for (auto& generator : generators) {
        Vector<double> omega = buildClosedPrimalOneForm(generator);
        Vector<double> rhs = diff0.transpose() * hodge1 * omega;
        auto laplaceMatrix = geometry->laplaceMatrix();
        Vector<double> alpha = solvePositiveDefinite(laplaceMatrix, rhs);
        Vector<double> gamma = omega - diff0 * alpha;
        gammas.push_back(gamma);
    }
    return gammas; // placeholder
}