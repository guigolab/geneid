package qalign.algo.msa;

/**

 * Die Beschreibung des Typs hier eingeben.

 * Erstellungsdatum: (17.01.2002 21:57:50)

 * @author: Hager-Dummy

 */

public class Heap {

	public int min, max;			// minimum and maximum buckets of heap

	public Edge[] bucket;			// buckets of edges

/**
 * Die Beschreibung der Methode hier eingeben.
 * Erstellungsdatum: (21.01.2002 18:37:12)
 * @param newMax int
 */
public Heap(int newMax) {

	bucket= new Edge[newMax+2];		// +1 (newMax mal upcount), +1 (length)
}
}
