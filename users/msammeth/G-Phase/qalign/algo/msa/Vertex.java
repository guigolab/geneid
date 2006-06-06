package qalign.algo.msa;

/**

 * lattice vertex

 * Erstellungsdatum: (17.01.2002 21:43:42)

 * @author: Hager-Dummy

 */

public class Vertex extends Coordinate {

	public Edge out; // outgoing edge adjacency list

	public CoordinateValues prev_coord_val; // father in array

	public Edge nonextracted_inedges;
	// incoming edges still not extracted from heap

	/**
	 * Insert the method's description here.
	 * Creation date: (20.01.2002 03:12:16)
	 * @param prev_coord_val align.CoordinateValues
	 */
	public Vertex() {
	}


	public String toString() {

		String result= "{";
		CoordinateValues cv= prev_coord_val;
		while (cv!= null) {

			result+= cv.value;
			cv= cv.curr_coord.prev_coord_val;
			if (cv!= null)
				result+= ", ";
		}
		result+= "}";

		return result;
	}
}
