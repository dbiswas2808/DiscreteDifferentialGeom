"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                let vertices = this.mesh.vertices
                for (let i = 0; i < vertices.length; i++) {
                        vertices[i].index = i;
                }
                
                let edges = this.mesh.edges
                for (let i = 0; i < edges.length; ++i) {
                        edges[i].index = i;
                }

                let faces = this.mesh.faces
                for (let i = 0; i < faces.length; ++i) {
                        faces[i].index = i;
                }

                // TODO
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                let tr = new Triplet(mesh.edges.length, mesh.vertices.length);
                let halfedges = mesh.halfedges;
                for (let he of halfedges) {
                        tr.addEntry(1, he.edge.index, he.vertex.index);
                }
                return SparseMatrix.fromTriplet(tr);
                // TODO
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                let tr = new Triplet(mesh.faces.length, mesh.edges.length);
                let faces = mesh.faces;
                for (let f of faces) {
                        for (let e of f.adjacentEdges()) {
                                tr.addEntry(1, f.index, e.index);
                        }
                }
                
                return SparseMatrix.fromTriplet(tr);
                // TODO
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                let vec = DenseMatrix.zeros(this.mesh.vertices.length);

                for (let v of subset.vertices) {
                        vec.set(1, v);
                }

                return vec;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                let vec = DenseMatrix.zeros(this.mesh.edges.length);

                for (let e of subset.edges) {
                        vec.set(1, e);
                }

                return vec;
                // TODO
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                let vec = DenseMatrix.zeros(this.mesh.faces.length);

                for (let f of subset.faces) {
                        vec.set(1, f);
                }

                return vec;
                // TODO
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                // TODO
                let newSubset = new MeshSubset();
                for (let v of subset.vertices) {
                        newSubset.addVertex(v);
                        let vertex = this.mesh.vertices[v];
                        for (let e of vertex.adjacentEdges()) {
                                newSubset.addEdge(e.index);
                        }

                        for (let f of vertex.adjacentFaces()) {
                                newSubset.addFace(f.index);
                        }
                }

                for (let e of subset.edges) {
                        newSubset.addEdge(e);
                        let edge = this.mesh.edges[e];
                        
                        for (let f of edge.adjacentFaces()) {
                                newSubset.addFace(f.index);
                        }
                }

                newSubset.addFaces(subset.faces);
                return newSubset; // placeholder
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                let newSubset = new MeshSubset();
                newSubset.addVertices(subset.vertices);

                for (let fIdx of subset.faces) {
                        let f = this.mesh.faces[fIdx];
                        newSubset.addFace(f.index);
                        for (let e of f.adjacentEdges()) {
                                newSubset.addEdge(e.index);
                        }

                        for (let v of f.adjacentVertices()) {
                                newSubset.addVertex(v.index);
                        }
                }

                for (let eIdx of subset.edges) {
                        let e = this.mesh.edges[eIdx];
                        newSubset.addEdge(e.index);

                        for (let v of e.vertices()) {
                                newSubset.addVertex(v.index);
                        }
                }

                // closure of vertex is the vertex itself
                return newSubset; // placeholder
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                let closure = this.closure(this.star(subset));
                closure.deleteSubset(this.star(this.closure(subset))); // placeholder
                return closure;
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                return subset.equals(this.closure(subset));
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                if (this.isComplex(subset)) {
                        let edgeVector = this.buildEdgeVector(subset);
                        let faceVector = this.buildFaceVector(subset);

                        if (subset.faces.size > 0){
                                let faceEdgeAdjacencyMatrix = this.A1;

                                for (let e of subset.edges) {
                                        let numEnts = faceEdgeAdjacencyMatrix.subMatrix(0, this.mesh.faces.length, e, e + 1).transpose().timesDense(faceVector);
                                        if (numEnts.get(0) == 0) {
                                                return -1;
                                        }
                                }
                        }
                        
                        if (subset.edges.size > 0) {
                                let edgeVertexAdjacencyMatrix = this.A0;
                                for (let v of subset.vertices) {
                                        let numEnts = edgeVertexAdjacencyMatrix.subMatrix(0, this.mesh.edges.length, v, v + 1).transpose().timesDense(edgeVector);
                                        if (numEnts.get(0) == 0) {
                                                return -1;
                                        }
                                }
                        }

                        if (subset.faces.size > 0) {
                                return 2;
                        } else if (subset.edges.size > 0) {
                                return 1;
                        } else if (subset.vertices.size > 0) {
                                return 0;
                        }
                }

                return -1;
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                let boundary = new MeshSubset();
                if (this.isPureComplex(subset)) {
                        let edgeVector = this.buildEdgeVector(subset);
                        let faceVector = this.buildFaceVector(subset);
                        let faceEdgeAdjacency = this.A1;
                        let edgeVertexAjacency = this.A0;

                        if (subset.faces.size > 0) {
                                for (let e of subset.edges) {
                                        let numEntity = faceEdgeAdjacency.subMatrix(0, this.mesh.faces.length, e, e + 1).transpose().timesDense(faceVector).get(0);
                                        if ( numEntity == 1) {
                                                boundary.addEdge(e);
                                        }
                                }
                                boundary = this.closure(boundary);
                        } else if (subset.edges.size > 0) {
                                for (let v of subset.vertices) {
                                        let numEntity = edgeVertexAjacency.subMatrix(0, this.mesh.edges.length, v, v + 1).transpose().timesDense(edgeVector).get(0);
                                        if ( numEntity == 1) {
                                                boundary.addVertex(v);
                                        }
                                }
                                boundary = this.closure(boundary);
                        }
                }

                return boundary; // placeholder
        }
}
