package qalign.algo.msa;

/**

 * Index in array.

 * Erstellungsdatum: (17.01.2002 21:30:34)

 * @author: Hager-Dummy

 */

public class CoordinateValues {

	public Coordinate next_coord; // next coordinate array

	public Coordinate curr_coord; // current coordinate array

	// value of this coordinate in absolute terms, used by coord(), safe_coord(), create_coord(), short-> int
	public int value;

	/**
	 * Null-Constructor.
	 * <main.create_coordinate()>
	 * Erstellungsdatum: (21.01.2002 15:41:15)
	 */
	public CoordinateValues() {

		curr_coord = null;
		next_coord = null;
	}
}
