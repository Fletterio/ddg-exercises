// Implement member functions for TrivialConnections class.
#include "trivial-connections.h"
#include "geometrycentral/numerical/linear_solvers.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
TrivialConnections::TrivialConnections(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
    geometry->requireVertexPositions();

    TreeCotree t = TreeCotree(mesh, geometry);
    t.buildGenerators();
    HarmonicBases hb = HarmonicBases(mesh, geometry);
    HodgeDecomposition hodge = HodgeDecomposition(mesh, geometry);

    this->generators = t.generators;
    this->bases = hb.compute(t.generators, hodge);

    // Build period matrix.
    this->P = this->buildPeriodMatrix();

    // TODO: Store DEC operators
    this->A = geometry->laplaceMatrix();
    this->hodge1 = geometry->buildHodgeStar1Form();
    this->d0 = geometry->buildExteriorDerivative0Form();
}

/*
 * Builds the period matrix Pij = ∑_{ek ∈ li} (ξj)k, where li is the ith homology generator, ek is a dual edge in li and
 * ξj is the jth harmonic 1-form basis.
 *
 * Input:
 * Returns: A sparse matrix represending the period matrix.
 */
SparseMatrix<double> TrivialConnections::buildPeriodMatrix() const {
    SparseMatrix<double> pm(this->generators.size(), this->generators.size());
    std::vector<Eigen::Triplet<double>> entries;
    for (int i = 0; i < this->generators.size(); i++) {
        for (int j = 0; j < this->bases.size(); j++) {
            const auto& generator = this->generators[i];
            const auto& base = this->bases[j];
            double integral = 0.0;
            for (const auto& he : generator) {
                double integratedFormOverEdge = base[he.edge().getIndex()];
                integral += integratedFormOverEdge * (he.edge().halfedge() == he ? 1.0 : -1.0);
            }
            entries.push_back({i, j, integral});
        }
    }
    pm.setFromTriplets(entries.cbegin(), entries.cend());
    return pm;
}

/*
 * Determine if a mesh satisfies Gauss-Bonnet.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: True if mesh satisfies Gauss-Bonnet, false otherwise.
 */
bool TrivialConnections::satsifyGaussBonnet(const Vector<double>& singularity) const {

    return (abs(singularity.sum() - geometry->eulerCharacteristic()) < 1e-8);
}

/*
 * Compute the dual 0-form potential β by solving the system d𝛿β = -K + 2π * singularity.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: The coexact component 𝛿β.
 */
Vector<double> TrivialConnections::computeCoExactComponent(const Vector<double>& singularity) const {

    //construct K
    Vector<double> K(mesh->nVertices());
    for (const auto& v : mesh->vertices()) {
        K[v.getIndex()] = geometry->angleDefect(v);
    }
    Vector<double> rhs = - K + 2 * PI * singularity;
    SparseMatrix<double> laplace = this->A;
    Vector<double> betaTilde = solvePositiveDefinite(laplace, rhs);
    return this->hodge1 * this->d0 * betaTilde;
}


/*
 * Given an initial angle αi in face i, this function computes the new angle αj in the neighboring face j as
 * αj = αi - θij + θji, where θij and θji are the angles between the shared edge e and an arbitrary but fixed reference
 * direction in faces i and j. Repeating this procedure for n consecutive dual edges in a generator gives a sequence of
 * angles α0, . . . , αn with a resulting total angle defect equal to αn - α0. This corresponds to transporting a vector
 * around a generator by unfolding, sliding and refolding it across neighboring faces without any extra in plane
 * rotation.
 *
 * Input: A halfedge lying on the shared edge between face i and j, and the initial angle αi.
 * Returns: The new angle αj.
 */
double TrivialConnections::transportNoRotation(Halfedge he, double alphaI) const {

    Vector3 u = geometry->halfedgeVector(he);

    // Compute two orthonormal tangent vectors for each face.
    Face fi = he.face();
    Face fj = he.twin().face();
    Vector3 e1 = geometry->halfedgeVector(fi.halfedge()).normalize();
    Vector3 e2 = cross(geometry->faceNormal(fi), e1);
    Vector3 f1 = geometry->halfedgeVector(fj.halfedge()).normalize();
    Vector3 f2 = cross(geometry->faceNormal(fj), f1);
    double thetaIJ = atan2(dot(u, e2), dot(u, e1));
    double thetaJI = atan2(dot(u, f2), dot(u, f1));

    return alphaI - thetaIJ + thetaJI;
}

/*
 * Compute the harmonic component γ = ∑_{i = 1, ..., 2g} zi ξi by solving the system Pz = v - ∑𝛿β.
 * v - ∑𝛿β should be normalized to lie between -π and π.
 *
 * Input: The coexact component 𝛿β.
 * Returns: The harmonic component γ.
 */
Vector<double> TrivialConnections::computeHarmonicComponent(const Vector<double>& deltaBeta) const {

    if (this->generators.size() == 0) {
        return Vector<double>::Zero(mesh->nEdges());
    }
    //construct rhs
    Vector<double> rhs(this->generators.size());
    for (int i = 0; i < generators.size(); i++) {
        const auto& generator = generators[i];
        double holonomy = 0.0;
        double integralOverGenerator = 0.0;
        for (const auto& he : generator) {
            holonomy = transportNoRotation(he, holonomy);
            integralOverGenerator += deltaBeta[he.edge().getIndex()] * (he.edge().halfedge() == he ? 1.0 : -1.0);
        }
        rhs[i] = - holonomy - integralOverGenerator;
        rhs[i] = rhs[i] - 2 * PI * floor(rhs[i] / (2 * PI)) - PI;
    }
    SparseMatrix<double> pm = this->P;
    Vector<double> z = solveSquare(pm, rhs);
    Vector<double> harmonicComponent = Vector<double>::Zero(mesh->nEdges());
    for (int i = 0; i < this->bases.size(); i++) {
        harmonicComponent += z[i] * bases[i];
    }
    return harmonicComponent;
}

/*
 * Compute the dual 1-form connections φ = 𝛿β + γ.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: A vector representing the connections.
 */
Vector<double> TrivialConnections::computeConnections(const Vector<double>& singularity) const {

    if (!this->satsifyGaussBonnet(singularity)) {
        std::cerr << "Singularities do not add up to the Euler characteristic of the mesh, which is: "
                  << geometry->eulerCharacteristic() << std::endl;
        return Vector<double>::Zero(mesh->nEdges());
    }
    Vector<double> dBeta = computeCoExactComponent(singularity);
    Vector<double> harmonic = computeHarmonicComponent(dBeta);
    return dBeta + harmonic;
}