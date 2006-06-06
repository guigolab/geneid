package qalign.algo.msa;

/**

 * Die Beschreibung des Typs hier eingeben.

 * Erstellungsdatum: (17.01.2002 21:46:07)

 * @author: Hager-Dummy

 */

public class Edge {

	public Vertex tail, head; // edge tail and head vertices

	public int dist; // distance to head from source along edge

	public int refer;
	// how many backtrack edges point to me; needed for detecting no longer used edges/vertices (S.16)

	public Edge next, prev; // edge adjacency list links

	public Edge heap_succ, heap_pred; // heap bucket links

	public Edge nonextracted_next, nonextracted_prev;
	// nonextracted inedges links

	public Edge backtrack;
	// edge to previous edge in path; preceeding edge in some optimal path (S.16)

	/**
	 * Insert the method's description here.
	 * Creation date: (20.01.2002 12:42:29)
	 */
	public Edge() {
	}
	/**
	 * Insert the method's description here.
	 * Creation date: (20.01.2002 12:25:31)
	 * @param newHead align.Vertex
	 */
	public Edge(Vertex newHead) {
	}


	public String toString() {

			// funny thing, edges are in backtracking orientation..
		return (tail.toString() + " --> " + head.toString());
	}

}
