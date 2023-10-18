// Implement member functions for TreeCotree class.
#include "tree-cotree.h"
#include <queue>
#include <set>

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
TreeCotree::TreeCotree(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
    vertexParent = std::map<Vertex, Vertex>();
    geometry->requireVertexIndices();
}

/*
 * Build a primal spanning tree on a mesh without boundary. More specifically, populate the member variable
 * <vertexParent>, which is a std::map that maps each vertex of the input mesh to its parent in the primal spanning
 * tree.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildPrimalSpanningTree() {

    Vertex root = *(mesh->vertices().begin());
    std::set<Vertex> visited = {root};
    vertexParent[root] = root;
    std::queue<Vertex> vertexQueue;
    vertexQueue.push(root);
    while (!vertexQueue.empty()) {
        Vertex parent = vertexQueue.front();
        for (const Vertex& v : parent.adjacentVertices()) {
            if (0 == visited.count(v)) {
                visited.insert(v);
                vertexParent[v] = parent;
                vertexQueue.push(v);
            }
        }
        vertexQueue.pop();
    }
}

/*
 * Check whether a halfedge is in the primal spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the primal spanning tree, false otherwise.
 */
bool TreeCotree::inPrimalSpanningTree(Halfedge he) {

    return he.tailVertex() == vertexParent[he.tipVertex()];
}

/*
 * Build a dual spanning tree on a mesh without boundary. More specificially, populate the member variable <faceParent>,
 * which is a std::map that maps each face of the input mesh to its parent in the dual spanning tree.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildDualSpanningCoTree() {

    Face root = *(mesh->faces().begin());
    std::set<Face> visited = {root};
    faceParent[root] = root;
    std::queue<Face> faceQueue;
    faceQueue.push(root);
    while (!faceQueue.empty()) {
        Face parent = faceQueue.front();
        for (const Halfedge& he : parent.adjacentHalfedges()) {
            if ((!inPrimalSpanningTree(he)) && (!inPrimalSpanningTree(he.twin())) && 0 == visited.count(he.twin().face())) {
                Face f = he.twin().face();
                visited.insert(f);
                faceParent[f] = parent;
                faceQueue.push(f);
            }
        }
        faceQueue.pop();
    }
}

/*
 * Check whether a halfedge is in the dual spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the dual spanning tree, false otherwise.
 */
bool TreeCotree::inDualSpanningCotree(Halfedge he) {

    return he.face() == faceParent[he.twin().face()];
}

/*
 * Returns a halfedge lying on the shared edge between face f and g.
 *
 * Input: Two adjacent faces <f> and <g>.
 * Returns: A halfedge lying on the shared edge between face f and g.
 */
Halfedge TreeCotree::sharedHalfedge(Face f, Face g) const {

    for (Halfedge he : f.adjacentHalfedges()) {
        if (he.twin().face() == g) {
            return he;
        }
    }
    // Should never get here!
    std::cerr << "Oops, TreeCotree::sharedHalfedge() received bad input." << std::endl;
    return f.halfedge();
}

/*
 * Compute the homology generators of the mesh.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildGenerators() {

    // order doesn't matter in a mesh without boundary
    buildPrimalSpanningTree();
    buildDualSpanningCoTree();

    //Build generators and populate this->generators
    //First we find all edges not in either tree
    std::vector<Edge> loopEdges;
    for (const auto& e : mesh->edges()) {
        if ((!inPrimalSpanningTree(e.halfedge())) && (!inDualSpanningCotree(e.halfedge())) &&
            (!inPrimalSpanningTree(e.halfedge().twin())) && (!inDualSpanningCotree(e.halfedge().twin()))) {
            loopEdges.push_back(e);
        }
    }
    for (const Edge& e : loopEdges) {
        auto FaceIt = e.adjacentFaces().begin();
        Face faceForward = *FaceIt;
        ++FaceIt;
        Face faceBackward = *FaceIt;
        std::vector<Face> forwardFaces, backwardFaces;
        do {
            forwardFaces.push_back(faceForward);
            faceForward = faceParent[faceForward];
        } while (faceForward != faceParent[faceForward]);
        do {
            backwardFaces.push_back(faceBackward);
            faceBackward = faceParent[faceBackward];
        } while (faceBackward != faceParent[faceBackward]);
        while (backwardFaces.back() == forwardFaces.back() &&
               backwardFaces[backwardFaces.size() - 2] == forwardFaces[forwardFaces.size() - 2]) {
            backwardFaces.pop_back();
            forwardFaces.pop_back();
        }
        backwardFaces.pop_back();
        auto FaceVecIt = backwardFaces.crbegin();
        while (FaceVecIt != backwardFaces.crend()) {
            forwardFaces.push_back(*FaceVecIt);
            FaceVecIt++;
        }
        //Here I have the whole loop of faces stores in forwardFaces
        std::vector<Halfedge> generator;
        auto FaceDualEdgeIt = forwardFaces.cbegin();
        Face previous = *FaceDualEdgeIt;
        FaceDualEdgeIt++;
        Face current = *FaceDualEdgeIt;
        FaceDualEdgeIt++;
        Face next;
        Edge fromPrevToCurr;
        Edge fromCurrToNext;
        while (FaceDualEdgeIt != forwardFaces.cend()) {
            next = *FaceDualEdgeIt;
            for (const Edge e : previous.adjacentEdges()) {
                for (auto it = current.adjacentEdges().begin(); it != current.adjacentEdges().end(); ++it) {
                    if (*it == e) fromPrevToCurr = e;
                }
            }
            for (const Edge e : current.adjacentEdges()) {
                for (auto it = next.adjacentEdges().begin(); it != next.adjacentEdges().end(); ++it) {
                    if (*it == e) fromCurrToNext = e;
                }
            }
            Halfedge genEdge = fromPrevToCurr.halfedge();
            if ((genEdge.tipVertex() != fromCurrToNext.firstVertex()) &&
                (genEdge.tipVertex() != fromCurrToNext.secondVertex())) {
                genEdge = genEdge.twin();
            }
            generator.push_back(genEdge);
            previous = current;
            current = next;
            FaceDualEdgeIt++;
        }
        //last two missing
        next = forwardFaces[0];
        for (const Edge e : previous.adjacentEdges()) {
            for (auto it = current.adjacentEdges().begin(); it != current.adjacentEdges().end(); ++it) {
                if (*it == e) fromPrevToCurr = e;
            }
        }
        for (const Edge e : current.adjacentEdges()) {
            for (auto it = next.adjacentEdges().begin(); it != next.adjacentEdges().end(); ++it) {
                if (*it == e) fromCurrToNext = e;
            }
        }
        Halfedge genEdge = fromPrevToCurr.halfedge();
        if ((genEdge.tipVertex() != fromCurrToNext.firstVertex()) &&
            (genEdge.tipVertex() != fromCurrToNext.secondVertex())) {
            genEdge = genEdge.twin();
        }
        generator.push_back(genEdge);
        previous = current;
        current = next;
        next = forwardFaces[1];
        for (const Edge e : previous.adjacentEdges()) {
            for (auto it = current.adjacentEdges().begin(); it != current.adjacentEdges().end(); ++it) {
                if (*it == e) fromPrevToCurr = e;
            }
        }
        for (const Edge e : current.adjacentEdges()) {
            for (auto it = next.adjacentEdges().begin(); it != next.adjacentEdges().end(); ++it) {
                if (*it == e) fromCurrToNext = e;
            }
        }
        genEdge = fromPrevToCurr.halfedge();
        if ((genEdge.tipVertex() != fromCurrToNext.firstVertex()) &&
            (genEdge.tipVertex() != fromCurrToNext.secondVertex())) {
            genEdge = genEdge.twin();
        }
        generator.push_back(genEdge);
        generators.push_back(generator);
    }
}