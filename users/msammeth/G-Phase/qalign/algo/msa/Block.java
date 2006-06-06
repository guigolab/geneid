package qalign.algo.msa;

/**
 * data structure for fixed points of each seq. pair --
 * (fixer.c)
 * Creation date: (19.01.2002 20:25:25)
 * @author:
 */
public class Block {

	public int i,j,len;		// i,j are startingpoints of fixed BLOCK, len is length
	public Block next;			// next points to next BLOCK for seq. pair

}
