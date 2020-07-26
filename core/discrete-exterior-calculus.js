"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		let diagonal = DenseMatrix.zeros(geometry.mesh.vertices.length);
		for (let vertex of geometry.mesh.vertices) {
			diagonal.set(geometry.barycentricDualArea(vertex), vertexIndex[vertex]);
		}

		return SparseMatrix.diag(diagonal);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		let diagonal = DenseMatrix.zeros(geometry.mesh.edges.length);
		for (let edge of geometry.mesh.edges) {
			let he = edge.halfedge;
			let dual = geometry.cotan(he) + geometry.cotan(he.twin);
			diagonal.set(dual / 2, edgeIndex[edge]);
		}

		return SparseMatrix.diag(diagonal);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		let diagonal = DenseMatrix.zeros(geometry.mesh.faces.length);
		for (let f of geometry.mesh.faces) {
			let dual = 1 / geometry.area(f);
			diagonal.set(dual, faceIndex[f]);
		}

		return SparseMatrix.diag(diagonal);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		let tr = new Triplet(geometry.mesh.edges.length, geometry.mesh.vertices.length);
		for (let e of geometry.mesh.edges) {
			let he = e.halfedge;
			tr.addEntry(1, edgeIndex[e], vertexIndex[he.vertex]);
			tr.addEntry(-1, edgeIndex[e], vertexIndex[he.twin.vertex]);
		}
		
		return SparseMatrix.fromTriplet(tr); // placeholder
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		let tr = new Triplet(geometry.mesh.faces.length, geometry.mesh.edges.length);
		for (let f of geometry.mesh.faces) {
			let he = f.halfedge;
			for (let i = 0; i < 3; i++) {
				tr.addEntry((he.index < he.twin.index ? 1 : -1), faceIndex[f], edgeIndex[he.edge]);
				he = he.next;
			}
		}

		return SparseMatrix.fromTriplet(tr); // placeholder
	}
}
