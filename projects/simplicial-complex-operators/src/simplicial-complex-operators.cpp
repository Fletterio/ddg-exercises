// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

//-------------------------AUX FUNCTIONS------------------------------------
bool empty_intersection(const std::set<size_t>& x, const std::set<size_t>& y) {
    std::set<size_t>::const_iterator i = x.begin();
    std::set<size_t>::const_iterator j = y.begin();
    while (i != x.end() && j != y.end()) {
        if (*i == *j)
            return false;
        else if (*i < *j)
            ++i;
        else
            ++j;
    }
    return true;
}

bool unique_intersection(const std::set<size_t>& x, const std::set<size_t>& y) {
    int count = 0;
    std::set<size_t>::const_iterator i = x.begin();
    std::set<size_t>::const_iterator j = y.begin();
    while (i != x.end() && j != y.end()) {
        if (*i == *j) {
            if (count++) return false;
            ++i;
        }  
        else if (*i < *j)
            ++i;
        else
            ++j;
    }
    return true;
}

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();

    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.

    SparseMatrix<size_t> adj(mesh->nEdges(), mesh->nVertices());
    std::vector<Eigen::Triplet<size_t>> entries;

    for (Edge e : mesh->edges()) {
        size_t edgeIdx = geometry->edgeIndices[e];
        entries.push_back(Eigen::Triplet<size_t>(edgeIdx, geometry->vertexIndices[e.firstVertex()], 1));
        entries.push_back(Eigen::Triplet<size_t>(edgeIdx, geometry->vertexIndices[e.secondVertex()], 1));
    }

    adj.setFromTriplets(entries.cbegin(), entries.cend());
    return adj;

}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    SparseMatrix<size_t> adj(mesh->nFaces(), mesh->nEdges());
    std::vector<Eigen::Triplet<size_t>> entries;

    for (Face f : mesh->faces()) {
        size_t faceIdx = geometry->faceIndices[f];
        auto firstHalfEdge = f.halfedge();
        auto currHalfEdge = firstHalfEdge;
        do {
            entries.push_back(Eigen::Triplet<size_t>(faceIdx, geometry->edgeIndices[currHalfEdge.edge()], 1));
            currHalfEdge = currHalfEdge.next();
        } while (currHalfEdge != firstHalfEdge);
    }
    adj.setFromTriplets(entries.cbegin(), entries.cend());
    return adj;

}

SparseMatrix<size_t> SimplicialComplexOperators::buildVertexFaceAdjacencyMatrix() const {
    return this->A1 * this->A0;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector<size_t> vertexVector = Vector<size_t>::Zero(mesh->nVertices());
    for (size_t vIdx : subset.vertices) {
        vertexVector[vIdx] = 1;
   }
    return vertexVector;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector<size_t> edgeVector = Vector<size_t>::Zero(mesh->nEdges());
    for (size_t eIdx : subset.edges) {
        edgeVector[eIdx] = 1;
    }
    return edgeVector;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    Vector<size_t> faceVector = Vector<size_t>::Zero(mesh->nFaces());
    for (size_t fIdx : subset.faces) {
        faceVector[fIdx] = 1;
    }
    return faceVector;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;

    for (size_t vIdx : subset.vertices) {
        vertices.insert(vIdx);
        for (Eigen::SparseMatrix<size_t, Eigen::ColMajor>::InnerIterator it(A0, vIdx); it; ++it) {
            edges.insert(it.row());
        }
        for (Eigen::SparseMatrix<size_t, Eigen::ColMajor>::InnerIterator it(A1A0, vIdx); it; ++it) {
            faces.insert(it.row());
        }
    }
    for (size_t eIdx : subset.edges) {
        edges.insert(eIdx);
        for (Eigen::SparseMatrix<size_t, Eigen::ColMajor>::InnerIterator it(A1, eIdx); it; ++it) {
            faces.insert(it.row());
        }
    }

    for (size_t fIdx : subset.faces) {
        faces.insert(fIdx);
    }

    return MeshSubset(vertices, edges, faces); 
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;

    SparseMatrix<size_t> A0T = A0.transpose();
    SparseMatrix<size_t> A1T = A1.transpose();
    SparseMatrix<size_t> A1A0T = A1A0.transpose();


    for (size_t fIdx : subset.faces) {
        faces.insert(fIdx);
        for (Eigen::SparseMatrix<size_t>::InnerIterator it(A1T, fIdx); it; ++it) {
            edges.insert(it.row());
        }
        for (Eigen::SparseMatrix<size_t>::InnerIterator it(A1A0T, fIdx); it; ++it) {
            vertices.insert(it.row());
        }
    }

    for (size_t eIdx : subset.edges) {
        edges.insert(eIdx);
        for (Eigen::SparseMatrix<size_t>::InnerIterator it(A0T, eIdx); it; ++it) {
            vertices.insert(it.row());
        }
    }

    for (size_t vIdx : subset.vertices) {
        vertices.insert(vIdx);
    }

    return MeshSubset(vertices, edges, faces);
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    auto ClSt = closure(star(subset));
    auto StCl = star(closure(subset));

    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;

    std::set_difference(ClSt.vertices.begin(), ClSt.vertices.end(), StCl.vertices.begin(), StCl.vertices.end(), std::inserter(vertices, vertices.begin()));
    std::set_difference(ClSt.edges.begin(), ClSt.edges.end(), StCl.edges.begin(), StCl.edges.end(), std::inserter(edges, edges.begin()));
    std::set_difference(ClSt.faces.begin(), ClSt.faces.end(), StCl.faces.begin(), StCl.faces.end(), std::inserter(faces, faces.begin()));

    return MeshSubset(vertices, edges, faces);
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    return subset.equals(closure(subset)); // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    if (! isComplex(subset)) {
        return -1;
    }

    bool canBe0Complex = subset.edges.empty() && subset.faces.empty();
    bool canBe1Complex = subset.faces.empty() && ! canBe0Complex;
    bool canBe2Complex = ! canBe1Complex;

    if (canBe2Complex) {
        for (size_t vIdx : subset.vertices) {
            std::set<size_t> adjacentFaces;
            for (SparseMatrix<size_t>::InnerIterator it(A1A0, vIdx); it; ++it) {
                adjacentFaces.insert(it.row());
            }
            if (empty_intersection(adjacentFaces, subset.faces)) {
                canBe2Complex = false;
                break;
            }
        }
        if (canBe2Complex) {
            for (size_t eIdx : subset.edges) {
                std::set<size_t> adjacentFaces;
                for (SparseMatrix<size_t>::InnerIterator it(A1, eIdx); it; ++it) {
                    adjacentFaces.insert(it.row());
                }
                if (empty_intersection(adjacentFaces, subset.faces)) {
                    canBe2Complex = false;
                    break;
                }
            }
        }
    }
    if (canBe2Complex) return 2;
    if (canBe1Complex) {
        for (size_t vIdx : subset.vertices) {
            std::set<size_t> adjacentEdges;
            for (SparseMatrix<size_t>::InnerIterator it(A0, vIdx); it; ++it) {
                adjacentEdges.insert(it.row());
            }
            if (empty_intersection(adjacentEdges, subset.edges)) {
                canBe1Complex = false;
                break;
            }
        }
    }
    if (canBe1Complex) return 1;

    return canBe0Complex ? 0 : -1;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    std::set<size_t> vertices, edges, faces;

    std::set<size_t> aux;

    int k = isPureComplex(subset);
    switch (k){ 
    case -1:
        std::cerr << "Only pure k-complexes have boundary!\n";
        return MeshSubset();

    case 0:
        return MeshSubset();

    case 1:
        for (size_t vIdx : subset.vertices) {
            std::set<size_t> adjacentEdges;
            for (SparseMatrix<size_t>::InnerIterator it(A0, vIdx); it; ++it) {
                adjacentEdges.insert(it.row());
            }
            if (unique_intersection(adjacentEdges, subset.edges)) vertices.insert(vIdx);
        }
        return MeshSubset(vertices, edges, faces);

    case 2:
        for (size_t eIdx : subset.edges) {
            std::set<size_t> adjacentFaces;
            for (SparseMatrix<size_t>::InnerIterator it(A1, eIdx); it; ++it) {
                adjacentFaces.insert(it.row());
            }
            if (unique_intersection(adjacentFaces, subset.faces)) edges.insert(eIdx);
        }
        return closure(MeshSubset(vertices, edges, faces));
    }


}