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
                let tr = new Triplet(this.mesh.edges.length, this.mesh.vertices.length);
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
                let tr = new Triplet(this.mesh.faces.length, this.mesh.edges.length);
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
                let newSubset = subset;
                for (let v of subset.vertices) {
                        let vertex = this.mesh.vertices[v];
                        for (let e of vertex.adjacentEdges()) {
                                newSubset.addEdge(e.index);
                        }

                        // vertex = this.mesh.vertices[v];
                        // for (let f of vertex.adjacentFaces()) {
                        //         newSubset.addFace(f.index);
                        // }
                }

                for (let e of newSubset.edges) {
                        let edge = this.mesh.edges[e];
                        let size = edge.vertices().size;
                        for (let i = 0; i < size; i++) {
                                newSubset.addVertex(edge.vertices[i].index);
                        }
                        
                        size = edge.adjacentFaces().size;
                        for (let i = 0; i < size; i++) {
                                newSubset.addFace(edge.adjacentFaces()[i].index);
                        }
                }

                for (let f of subset.faces) {
                        let face = this.mesh.faces[f];
                        for (let e of face.adjacentEdges()) {
                                newSubset.addEdge(e.index);
                        }
                        for (let v of face.adjacentCorners()) {
                                newSubset.addVertex(v.index);
                        }
                }

                return newSubset; // placeholder
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                let newSubset = new MeshSubset();
                // newSubset = subset;
                // for (let f of subset.faces) {
                //         for (let e of f.adjacentEdges()) {
                //                 newSubset.addEdge(e);
                //         }

                //         for (let v of f.adjacentCorners()) {
                //                 newSubset.addVertex(v);
                //         }
                // }

                // for (let e of subset.edges) {
                //         newSubset.addVertex(e.halfedge.vertex);
                //         newSubset.addVertex(e.halfedge.twin.vertex);  
                // }

                // closure of vertex is the vertex itself

                return newSubset; // placeholder
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                return this.closure(this.star(subset)).deleteSubset(this.star(this.closure(subset))); // placeholder
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
                let newSubset = new MeshSubset();
                // if (subset.faces.size != 0) {
                //         for (let f of subset.faces) {
                //                 newSubset.addFace(f);
                //                 newSubset.addVertices(f.adjacentCorners());
                //                 newSubset.addEdges(f.adjacentEdges());
                //         }
                //         if (subset.equals(newSubset)) {
                //                 return 2;
                //         }
                // } else if (subset.edges.size != 0) {
                //         for (let e of subset.edges) {
                //                 newSubset.addEdge(e);
                //                 newSubset.addVertices(e.vertices());
                //         }
                //         if (subset.equals(newSubset)) {
                //                 return 1;
                //         }
                // } else if (subset.vertices.size != 0) { 
                //         newSubset.addVertices(e.vertices());
                        
                //         if (subset.equals(newSubset)) {
                //                 return 0;
                //         }
                // }

                return -1;
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                let newSubset = new MeshSubset();
                // for (let v of subset.vertices) {
                //         let count = 0;
                //         for (let e of v.adjacentEdges) {
                //                 subset.edges.has(e);
                //                 count++;
                //                 if (count == 2) {
                //                         break;
                //                 }
                //         }

                //         if (count == 1) {
                //                 newSubset.addVertex(v);
                //         }
                // }

                // for (let e of subset.edges) {
                //         let count = 0;
                //         for (let f of e.adjacentFaces) {
                //                 subset.faces.has(f);
                //                 count++;
                //                 if (count == 2) {
                //                         break;
                //                 }

                //         }

                //         if (count == 1) {
                //                 newSubset.addEdge(e);
                //         }
                // }

                // TODO
                return this.closure(newSubset); // placeholder
        }
}
