package qalign.algo.msa;

/**

 * Array for accessing vertices by lattice coordinate.

 * Erstellungsdatum: (17.01.2002 21:24:42)

 * @author: Hager-Dummy

 */

public class Coordinate {

	public int lo, hi;							// lower and upper limits on array indices

	public CoordinateValues[] coord_vals;		// next coordinate array indices, <create_vertex> IS ARRAY!!!

	public CoordinateValues  prev_coord_val;	// previous coordinate array index

	public Coordinate next_on_free_list;		// maintains available records

	public int refer;							// how many valid coordinate values do I have




}
