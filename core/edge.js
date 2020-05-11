"use strict";

class Edge {
	/**
	 * This class represents an edge in a {@link module:Core.Mesh Mesh}.
	 * @constructor module:Core.Edge
	 * @property {module:Core.Halfedge} halfedge One of the halfedges associated with this edge.
	 */
	constructor() {
		this.halfedge = undefined;
		this.index = -1; // an ID between 0 and |E| - 1, where |E| is the number of edges in a mesh
	}

	/**
	 * Checks whether this edge lies on a boundary.
	 * @method module:Core.Edge#onBoundary
	 * @returns {boolean}
	 */
	onBoundary() {
		return (this.halfedge.onBoundary || this.halfedge.twin.onBoundary);
	}

	/**
	 * Gets the adjacent faces of the edge
	 * @method model:Core.Edge#adjacentFaces
	 * @returns {array of Face}
	 */
	adjacentFaces() {
		let adajcentFaces = [];
		if (!this.halfedge.onBoundary) {
			adajcentFaces.push(this.halfedge.face);
		}

		if (!this.halfedge.twin.onBoundary) {
			adajcentFaces.push(this.halfedge.twin.face);
		}

		return adajcentFaces;
	}

	/**
	 * Gets the vertices of the edge
	 * @method model:Core.Edge#vertices
	 * @returns {array of vertices}
	 */
	vertices() {
		let vertices = [];
		vertices.push(this.halfedge.vertex);
		vertices.push(this.halfedge.twin.vertex);
		return vertices;
	}

	/**
	 * Defines a string representation for this edge as its index.
	 * @ignore
	 * @method module:Core.Edge#toString
	 * @returns {string}
	 */
	toString() {
		return this.index;
	}
}
