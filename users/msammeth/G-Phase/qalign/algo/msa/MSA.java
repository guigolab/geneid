package qalign.algo.msa;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import java.util.StringTokenizer;

import qalign.algo.CancelException;
import qalign.gui.GUIProxy;


/**
 * This is a 1:1 porting of the

import qalign.gui.MSAProxy;


/**
 * This is a 1:1 porting of the C-code.
 * (meaning it runs, scores and crashes consistent with the original version).
 * <br>
 * <b>May the force be with you.</b>
 * <br>
 * Bugs:
 * - in dist(): division result cast to (int) and then re-cast to (float).
 * 				(consistent with c-code).
 * - in newal(): 'BIG not big enough'. Nice bug fix - the real bug is within the
 * 				  nested loops, cause increasing BIG does not fix the problem.
 * 				  (checked, completely consistent with the original).
 *
 * @see #dist()
 * @see #newal()
 *
 * @author: NCBI, (porting by Michael Sammeth)
 */

public class MSA {

	protected boolean already= false;

		// Constants
	private final static int TRUE= 1;					// boolean truth
	private final static int FALSE= 0;					// boolean falsity
	private final static short SIGMA = 128;				// alphabet size
	private final static byte DIAG = 0;					// code for traceback
	private final static byte VERT = 1;					// code for traceback
	private final static byte HORZ = 2;					// code for traceback
	private final static String USAGE=
		"Usage:	 msa [-m] [-g] [-b]	[-o] [-d<delta>] [-e <epsilon file>] [-c <cost file>] [-f <fixer file>]	<input filename>";


	/**
	 * The GUI interface.
	 */
	protected GUIProxy proxy= null;

	/**
	 * For user cancels.<br>
	 * Checked in <code>sub()<code>.
	 */
	protected boolean cancel= false;

	/**
	 * Big value to init utopic upper bounds.<br>
	 * Is 999999 in original c-code, was set to 99999 for DCA tests
	 * (personal communication with Stoye), in order to reproduce the
	 * published results.<br>
	 * Does influence the results found (consistent with the original C-code !).
	 * E.g., BIG= 99999  for kinases12 -b set crashes with ArrayIndexOutOfBounds(),
	 * whilst BIG= 999999 runs ok. This is also explaining why it sometimes is
	 * not a good idea to set it to Integer.MAX_VALUE (doesn't fix the 'BIG not
	 * big enough' bug in newal() anyway up to my experience- try kinases 12 with
	 * -b -g option and pam250).<br>
	 * Made not <code>final</code> to be changeable by subclasses.
	 */
	protected int BIG = 999999;

	/**
	 * The character used for gaps.
	 */
	//protected char DASH = '-';
	
	protected char getDASH(){
		return '-';
	}

	/**
	 * Line length of output.<br>
	 * Made not <code>final</code> to be changeable by subclasses.
	 */
	protected byte LINE = 75;

	/**
	 * Max length of input sequences.<br>
	 * Made not <code>final</code> to be changeable by subclasses.
	 * Must be <code>static</code> for use in main().
	 */
	public static short LENGTH = Short.MAX_VALUE;

	/**
	 * Maximum number of input sequences.<br>
	 * Made not <code>final</code> to be changeable by subclasses.
	 * Must be <code>static</code> for use in main().
	 */
	public static int NUMBER = 40;

	/**
	 * Maximum size of input filename.<br>
	 * Made not <code>final</code> to be changeable by subclasses.
	 */
	protected byte FILE_NAME_SIZE = 20;

	/**
	 * Maximum number of command line arguments.<br>
	 * Made not <code>final</code> to be changeable by subclasses.
	 * Must be <code>static</code> for use in main().
	 */
	protected static byte MAX_ARGS = 13;

	/**
	 * The input sequences.
	 */
	protected char[][] S = new char[NUMBER	+1][];

	/**
	 * The Number of input sequences.
	 */
	protected int K;

	/**
	 * The lengthts of input sequences.
	 */
	protected int[] N = new int[NUMBER +1];

	/**
	 * Freedom of heuristical 'tube'. The difference between upper and lower bound.
	 *
	 * @see Upper
	 * @see Lower
	 */
	protected int delta= (-1);

	/**
	 * Differences between projected and pairwise costs
	 * of each combination of the input sequences.
	 */
	protected int[][] epsi = new int[NUMBER +1][NUMBER +1];

	/**
	 * Pairwise cost weight scales.
	 */
	protected int[][] scale= new int[NUMBER +1][NUMBER +1];


	/**
	 * Consistency check (sort of matching array).
	 */
	protected int[][][] Con=
			new int[NUMBER +1][NUMBER +1][LENGTH +1];

	/**
	 * Projected costs (pairwise projections of the hyperspace).
	 */
	protected int[][] proj= new int[NUMBER +1][NUMBER +1];

	/**
	 * Name of cost file, inited in main() (max: char[20]).
	 */
	protected String cname;

	/**
	 * Upper bound on alignment distance.
	 *
	 * @see #delta
	 * @see #Lower
	 */
	protected int Upper;


	/**
	 * Lower bound on alignment distance.
	 *
	 * @see #delta
	 * @see #Upper
	 */
	protected int Lower;

//
	/**
	 * Command line flag:<br>
	 *  0= unweighted summing (DO NOT calculates evolutionary tree)<br>
	 *  1= do calculate guiding tree<br>
	 *  (default= 0)
	 */
	protected int bflag= 0;
	/**
	 * Command line flag:<br>
	 *  0= end gaps PENALIZED (charged the same as internal ones)<br>
	 *  1= do not penalize terminal gaps<br>
	 *  (default= 1)
	 */
	protected int gflag= 1;
	/**
	 * Command line flag:<br>
	 *  0= no forced positions provided<br>
	 *  1= file with forced positions provided<br>
	 *  (default= 0)
	 */
	protected int fflag= 0;
	/**
	 * Command line flag:<br>
	 *  0= suppress output<br>
	 *  1= enable output<br>
	 *  (default= 1)
	 */
	protected int oflag= 1;

	/**
	 * Forward diagonal distance.
	 */
	protected int[][] dd= new int[LENGTH+1][];

	/**
	 * Forward horizontal distance
	 */
	protected int[][] hh= new int[LENGTH+1][];

	/**
	 * Forward vertical distance
	 */
	protected int[][] vv= new int[LENGTH+1][];

	/**
	 * Virtual vertex before source of hyperspace
	 * to assign to tail of first edge.
	 */
	protected Vertex presource;

	/**
	 * Symbol distance (i.e., cost matrix).
	 */
	protected short[][] D= new short[SIGMA][SIGMA];

	/**
	 * Altschul gap counts. Array to determine pairwise gap openings or extensions.
	 */
	protected int[][][][] T= new int[3][3][3][3];

	/**
	 * Cuts off first two dimensions of T, so simple gap detection ?!
	 */
	protected int[][][][] Tpointer=
			new int[NUMBER+1][NUMBER+1][3][];

	/**
	 * Gap opening costs.
	 *
	 * @ see #GG
	 */
	protected int G;

	/**
	 * Terminal Gap costs ('0' for not penalised).
	 *
	 * @ see #G
	 * @ see #gflag
	 */
	protected int GG;

	/**
	 * Priority queue to make the way to the hyperspace.
	 */
	protected Heap h;

	/**
	 * Fixed positions.
	 *
	 * @see #fix(int, int[], int[], int)
	 */
	protected Block[][] pFIX=
					new Block[NUMBER + 1][NUMBER + 1]; 	// BLOCKS, indexed by seq. pair <fixer.c>

	/**
	 * Marking Root node in guiding tree.
	 *
	 * @see #INOD
	 * @see #bias()
	 */
	protected final int ROOT= -9;						// <bias> (internal) root node's seq# (to be distinctable from 'normal' internals)

	/**
	 * Marking an Internal node's sequence number in guiding tree.
	 *
	 * @see #ROOT
	 * @see #bias()
	 */
	protected final int INOD= -1;

	/**
	 * @see #bias()
	 * @see #init_data()
	 */
	protected Node[] node_index;

	/**
	 * Pairwise distances in guiding tree.
	 * @see #bias()
	 * @see #init_data()
	 */
	protected int[][] Dij= new int[NUMBER][];

	/**
	 * Nodes constructed for guiding tree.
	 *
	 * @see #bias()
	 */
	protected int indexlen;

	/**
	 * All nodes in guiding tree.
	 *
	 * @see #bias()
	 * @see Node.#Node()
	 */
	protected int count;

	/**
	 * @see #br_len(int, int)
	 */
	protected int count2;

	/**
	 * @see #bias()
	 * @see #whts()
	 * @see #trace()
	 */
	protected float[][] B2= new float[NUMBER][NUMBER];		// super-local renamed due to already defined

	/**
	 * @see #bias()
	 * @see #whts()
	 */
	protected Node[] pN;									// super-local

	/**
	 * Minimal value for epsilon (heuristics).
	 *
	 * @see #ecalc()
	 */
	protected int MINE= 5;

	/**
	 * Maximal value for epsilon (heuristics).
	 *
	 * @see #ecalc()
	 */
	protected final int MAXE= 99999;

	/**
	 * Associative array to efficiently encode <code>char</code>s
	 * to </code>int<code>s.
	 *
	 * @see #ecalc()
	 */
	protected final char[] Num= 							// Num converted from int[] to char[]
		{'0','1','2','3','4','5','6','7','8','9',
		'A','B','C','D','E','F','G','H','I','J','K','L','M',
		'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

	/**
	 * Re-sorted input sequences (sorted according to progressive aligning order).
	 *
	 * @see #ecalc()
	 */
	protected char[][] SS= new char[NUMBER+1][];			// 2D 'cause related to S

	/**
	 * @see #ecalc()
	 */
	protected int[] NN= new int[NUMBER+1];
//
	/**
	 * @see #ecalc()
	 */
	protected byte[][] AL= new byte[NUMBER+1][2*LENGTH];	// 'alignment': layout of CHARS, computed from NAL[][] in ecalc()

	/**
	 * @see #ecalc()
	 */
	protected int[][] NAL= new int[NUMBER+1][2*LENGTH];	// 'new alignment': contains alignment positions of set of segments, updated in change with OAL[][] by align(), newal()

	/**
	 * @see #ecalc()
	 */
	protected int[][] OAL= new int[NUMBER+1][2*LENGTH];	// 'old alignment': contains aligning positions of segments, progressively updated in change with NAL[][] by align(), newal()
//
	/**
	 * @see #ecalc()
	 */
	protected int[] Ord= new int[NUMBER+1];

	/**
	 * Faces of lattice.<br>
	 * <br>
	 * First index of face is first sequence, second index of face is second
	 * sequence. Third index of face is position in first sequence.
  	 * Value of face item is set of positions in second sequence that are
  	 * consistent with the position in the first sequence based on the
	 * pairwise alignment of the two sequences. "Consistent" means that
 	 * this pair of positions was of use in the dynamic programming graph
 	 * for the pairwise alignment.
 	 *
 	 * @see #faces()
 	 */
	protected Row[][][] face=
				new Row[NUMBER+1][NUMBER+1][LENGTH+1];

	/**
	 * Atomary pairwise costs (i.e., cost table, matrix).
	 *
	 * @see #primer()
	 * @see #faces()
	 * @see #ecalc()
	 */
	protected int[][] costs= new int[NUMBER +1][NUMBER +1];
	
	protected int alicost= -1;

	/**
	 * @see #source()
	 */
	protected Coordinate[] A2= new Coordinate[NUMBER+2];	// super-global, 2nd definition of A !!!

	/**
	 * Reincarnation of Vertices.
	 *
	 * @see #create_vertex(CoordinateValues)
	 */
	protected Vertex avail_vertex= null;

	/**
	 * Reincarnation of Edges.
	 *
	 * @see #create_edge(Vertex,Vertex)
	 */
	protected Edge avail_edge= null;

	/**
	 * Reincarnation of Coordinates.
	 *
	 * @see #create_coordinate(int[],int,CoordinateValues)
	 */
	protected Coordinate avail_coordinate= null;

	/**
	 * @see #column(int[], int[])
	 */
	protected char[][] M= new char[NUMBER+1][LINE];

	/**
	 * @see #column(int[], int[])
	 */
	protected int C= 0;

	/**
	 * @see #ecalc()
	 * @see #align(int,int,int)
	 * @see #newal(int,int,int,int)
	 */
	protected char[] A3;									// super-global

	/**
	 * @see fixread(String)
	 */
	protected final int MAXLINE= 250;



/**
 * Generate adjacent vertices within region.
 * <main>
 * Creation date: (20.01.2002 17:21:36)
 * @param e align.Edge incoming edge u(p) --e--> v(q) (outgoing edges of v to be computed)
 * @param q int[] Point (Coordinates) of vertex u (e.tail)
 */
protected void adjacent(Edge e, int[] q) {

	// uses Coordinate[] A2
	int I, i;
	int[] index= new int[LENGTH + 1];
	int n;
	Coordinate a; // *a (!)
	CoordinateValues f; // *f (!)
	int[] P= new int[NUMBER + 1]; // coordinates of generated vertex
	int[] B= new int[NUMBER + 1]; // bound on coordinates of generated vertex

				// GET COORDINATES of e.head= v(q) and GENERATE adjacent points to q
	for (i= K, f= e.head.prev_coord_val;	// vertex is tail of coordinate-ll
		i >= 1;								// iterate all Dimensions for all Coordinates
		i--, f= f.curr_coord.prev_coord_val) {
		A2[i]= f.curr_coord;  			// get all Coordinates of the Vertex e.head
		B[i]= (P[i]= q[i]) + 1; 		// P: q-point; B:point max (+1) in each dimension
	}

				//
	for (I= K; I <= K;) {

		for (; I >= 1; I--)		// find first dimension (top-down)
			if (P[I] < B[I]) 	// which is still below bound value (can only be equal, not greater!?!)
				break;
		if (I < 1)				// if we reached bound values for all dimensions
			return;				// we are done (all adjacent edges are generated)



		if ((i= P[I]= B[I]) >= (a= A2[I]).lo							// first test is to see whether i is in the proper range,
				&& i <= a.hi
				&& (f= a.coord_vals[i - a.lo]).curr_coord != null) {	// second test is to see whether i is a valid coordinate

			if (f.next_coord == null) 									// if(child in trie does not exist)
				if (I < K)													// if(not yet generated last dimension of this branch)
					if ((n= intersect(P, I + 1, index)) > 0)					// if(there is a valid intersection)
						f.next_coord= create_coordinate(index, n, f);				// create new coordinate for this dimension
					else {														// else(no valid intersection)
						f.curr_coord= null;											// remove this dimension
						if (free_coordinate(a) != 0)								// if(there still are ungenerated adjacencies)
							I--;													// backstep one dimension
						continue;													// and try next loop
					}
				else													// else(trie-branch was completely generated spanning all dimensions)
					f.next_coord= (Coordinate) create_vertex(f); 			// chain Vertex @ trie-coordinates to end of ll

			for (a= A2[++I]= f.next_coord;								// if and only if first generated coordinate was valid
				I <= K;													// for all other dimensions coordinates
				a= A2[++I]= f.next_coord) 								// do the same...
				if ((i= P[I]= B[I] - 1) >= a.lo
					&& i <= a.hi
					&& (f= a.coord_vals[i - a.lo]).curr_coord != null) {
					if (f.next_coord == null)
						if (I < K)
							if ((n= intersect(P, I + 1, index)) > 0)
								f.next_coord= create_coordinate(index, n, f);
							else {
								f.curr_coord= null;
								if (free_coordinate(a) != 0)
									I--;
								break;
							}
						else
							f.next_coord= (Coordinate) create_vertex(f); // !!!??!!!
				} else
					break;

		} // if(i is ok)

				// CREATE the edge to newly generated, adjacent vertex
		if (I == (K + 1)) {						// last nextCoordinate* of an Coordinate-Set
			create_edge(e.head, (Vertex) a); 	// points to the vertex @ these coordinates
			I= K;								// continue loop
		}

	} // while(within dimensions)

} // adjacent()
/**
 * Subroutine to calculate an optimal progressive
 * alignment of the segments in the order chosen.<br>
 * <ecalc>
 * FIXED A3= SS + offset (low) for start; delegated to newal()
 * since there's the only reference of it (read-only)
 * CHANGED A3= SS[..+low] adding low offset delegated
 * to subroutine newal()
 * @globalReturn NAL[][] new alignment layout
 * (calculated one after another with OAL[][] with align())
 * Creation date: (20.01.2002 01:14:10)
 * @return int length of layout comprising all aligned segments
 * @param oldp1 int start of segment
 * @param p1 int end of segment
 * @param len int maximal allowed length of aligned segments' layout
 * @throws Exception (passing on from newal)
 * @see #ecalc()
 * @see #newal(int,int,int,int)
 */
protected int align(int oldp1, int p1, int len) throws Exception {

	int I,J,	// seqs loop counters
		i,		// position loop counter
		low,	// starting position of segment in seq Ord[I]
		m;		// actual length of layout

	for (I=1;I<=K;++I)
		OAL[I][0]=0;

	A3=SS[1];						// load Ord[1] seq; +low (start-offset) delegated to newal()
	low=Con[1][Ord[1]][oldp1];		// start of segment in seq Ord[1] (@ which what pos the start-pos in seq 1 was aligned)
	m=Con[1][Ord[1]][p1]-low-1;		// length of segment in seq Ord[1]

	for (i=1;i<=m;++i)
		NAL[1][m-i]=low+i;			// init backward new alignment with segments aligning pos's with seq Ord[1]

	for (I=2;I<=K;++I) {

		for (J=1;J<I;++J) 							// copy last determined new alignment
			for (i=1;i<=m;++i) 						// to old alignment
				OAL[J][i]=NAL[J][m-i]; 				// (mirrored!!!)

		A3=SS[I];									// new seq of which the segment is to be aligned
		low=Con[1][Ord[I]][oldp1];					// starting pos of segment in seq Ord[I]

		m=newal(I,low,Con[1][Ord[I]][p1]-low-1,m); 	// returns new length of layout

		if (m>len && I<K) 							// layout exceeds max. length (= |longest seq|*2 -2)
			fatal("Heuristic alignment segment is longer than "+len);
	}

	return m; 	// length of (complete) layout of aligned segments

}
/**
 * MSA subroutines.
 * Version 1.0   June 27, 1989
 * <br>
 * Routines for computing pairwise weights for SP alignment.
 * One first estimates an unrooted tree and then uses the topology
 * (Rationale-1) or the topology and branch lengths (Rationale-2)
 * to determine the weights.<br>
 * <br>
 * ---------------------------------------------------------------
 * <br>
 * Tree Program written by R Miner, with assistance by DJ Lipman<br>
 * <br>
 * Program to make unrooted evolutionary tree from pairwise distance
 * matrix using the Neighbor Joining method.<br>
 * <br>
 * See:	Saitou & Nei, Mol. Biol. Evol. 4 (1987) 406-425<br>
 * <br>
 * Start with array (node_index) of NODES, where each NODE is an operational
 * taxonomic unit (OTU).  Array indicates starting assumption of tree
 * topology.  Find best neighbors and coalesce -> remove from array, and
 * replace with NODE which is ancestor of the neighbors.  Continue this
 * process until there are only 3 NODES in node_index (only one topology
 * possible with 3 OTU's, so problem is solved).  Kludge slightly to
 * form into binary tree with appropriate root.  NOTE: Root has no biological
 * significance, but is for convenience in dealing with structure.
 * Nodes are labelled with 1) positive numbers corresponding to original
 * distance matrix 2) Negative one's  corresponding to internal nodes
 * 3) Negative 9 corresponding to ancestral node.<br>
 * <br>
 * Output is preorder traversal of tree (depth first).<br>
 * See associated program "read.c" for input routine of tree structure<br>
 * <bias.c>
 * Creation date: (19.01.2002 20:39:33)
 */
protected void bias() {

	Node ancestor;
	float len;

						// read in dist. matrix, initialize starting data structure
	init_data();

						// determine optimum node pairs, and make them neighbors
	minimize_Sij();

						// finish up turning initial data structure into binary tree
						// by forming root node
	ancestor = makenode(ROOT,(float)0.0,node_index[0],node_index[1]);
	count2 = 1;
	len = dist(0,1);
	len -= subdist(node_index[0],(float)0.0) / count2;
	count2 = 1;
	len -= subdist(node_index[1],(float)0.0) / count2;
	node_index[0].ltop = len;

						//	Print out rooted tree
	rpt(ancestor);
//	printf("\n");

						//	Calculate weights
	whts();

}
/**
 * computes the least squares estimate for the length of the
 * edge between node_index[i] and node_index[j].
 * <bias.c>
 * Creation date: (19.01.2002 21:59:42)
 * @return float
 * @param i int indices to node_index
 * @param j int
 */
protected float br_len(int i, int j) {

	float diz=0, djz=0;		// arithmetic middle distance from i,j to the others
	int t;					// seq loop counter
	count2= 1;				// init of global variable

	for (t=0; t<indexlen; ++t)
		if (t!=i && t!=j) {
			diz += dist(i,t);	// sum up distances from i to all other seqs
			djz += dist(j,t);	// and j to others
		}
	diz = diz / (indexlen - 2); // arithmetic middle (all - i - j)
	djz = djz / (indexlen - 2);

	float fff= ((dist(i,j) + diz - djz)/2 					// half distance of i-j middled by diz-diz
		- subdist(node_index[i], (float) 0.0)/count2);	// distance from Node i to leaves (count2 is altered recursively!)
//	System.out.println("br_len("+i+","+j+")= "+dist(i,j)+"+"+diz+"-"+djz+"- "+subdist(node_index[i], 0f)+" => "+fff);
	return fff;
}
/**
 * creates an internal node with node_index[i] and
 * node_indes[j] as children, and replaces node_index[i]
 * and node_index[j] with a single pointer to the new node
 * <bias.c>
 * Creation date: (19.01.2002 22:21:25)
 * @param i int indices to node_index
 * @param j int
 */
protected void coalesce(int i, int j) {

	Node par;	// new parent node

	node_index[i].ltop= br_len(i,j);	// get the distance from children
	node_index[j].ltop= br_len(j,i);	// to their new father

	par= new Node(INOD, (float) 0.0, node_index[i], node_index[j]);	// seq#= INOD= (-1)= internal Node
	node_index[i] = par;				// create new father and save it @ index i

	node_index[j] = node_index[indexlen-1];	// delete Node @ index j
	--indexlen;

}
/**
 * Compute	multiple sequence alignment	column for edge	<p,q>.
 * <main>
 * Erstellungsdatum: (21.01.2002 10:51:41)
 * @param p int[]
 * @param q int[]
 */
protected void column(int[] p, int[] q) {

	int	k;

	for	(k=1;k<=K;k++)
		  M[k][C] =	S[k][q[k] *	(q[k]-p[k])];
	if (++C	>= LINE)
		  output();
}
/**
 * Computes the sum of path lengths resulting from making
 * nodes node_index[i] and node_index[j] neighbors
 * <bias.c>
 * Creation date: (19.01.2002 22:11:08)
 * @return float
 * @param i int indices to node_index
 * @param j int
 */
protected float compute_S(int i, int j) {

	int t, tt;			// loop counters
	float s1=0, s2=0;	// sums of distances

	for (t=0;t<indexlen;t++)
		if (t!=i && t!=j) 			// sum up the distances for i,j to all others
			s1+= dist(i,t)+dist(j,t);
	s1= s1 / (2 * (indexlen- 2));	// calculate half of (for adding up both) arithmetic middle

	for (t=0; t<indexlen; t++)
		for (tt=t+1; tt<indexlen; tt++)	// sum up the distances of all (but i,j!) to each other
			if (t!=i && t!=j && tt!=i && tt!=j)
	     		s2+= dist(t,tt);
	s2= s2 / (indexlen- 2);			// kind of arithmetic middle ??!

	float fff= (s1 + s2 + dist(i,j) / 2);
//	System.err.println("computeS("+i+","+j+") "+fff);
	return fff; // add up dist(i,j) and half it for both seqs
}
/**
 * Calculate percent mismatch between sequences and convert it to
 * a distance measure for use in calculating an evolutionary tree.<br>
 * Traceback of 2D-Alignment (input from vv[], hh[], dd[]),
 * init of consistency matrix (Con[][][]).<br>
 * Creation date: (19.01.2002 19:46:15)
 * @return int (evolutionary) distance measure
 * @param I,J int seq# of pairwise alignment
 * @param n,m int length of the pairwise compared seqs
 */
protected int convert(int I, int J, int n, int m) {
/*
				// Evolutionary distances based on PAM model of molecular evolution
	ifdef PAMDIST		// PAMDIST commented out
	int[] ED= {...} 	// new int[101] filled with corresponding values
*/
	int i,j,V,H,M;
	int dir= DIAG;		// Direction of previous edge in traceback
	int match= 0;		// Number of matches in an optimal alignment
	float f,g;

	for(i=n,j=m;((i!=0)||(j!=0));) {	// init traceback in corner of aMatrix, count down ends if both are counted down
		V=vv[i][j]- (dir==VERT ? (j==m ? GG:G) : 0); 	// determine V.ertical cost without gap opening cost (added from beginning)
		H=hh[i][j]- (dir==HORZ ? (i==n ? GG:G) : 0); 	// determine H.orizontal real cost
		M= min3(V,H,dd[i][j]);							// get minimum
		if (j==0||M==V) {							// case 1: 2nd seq is empty OR vertical edge has lowest value
			dir=VERT;									// traceback vertical
			--i;										// count down aligned pos
		} else
			if (i==0||M==H) {						// case 2: 1st seq is empty OR horiz. edge has lowest value
				dir= HORZ;								// traceback horizontal
				--j;									// count down...
			} else {								// case 3: both seqs have chars left AND diagonal has lowest value
				dir= DIAG;								// traceback diagonal
				match+= S[I][i]==S[J][j] ? 1:0;			// boolean to int adaption of "S[][]==S[][]"
				Con[I][J][i]= j;						// mark in each consitency combination in the aligned pos
				Con[J][I][j]=i;							// the pos of the other seq, which what it was aligned
				--i;
				--j;									// count down both seq (diag!)
			}
	} // end of for(traceback)

/* 	ifdef PAMDIST			// PAMDIST commented out
		i=f= (100.0*(n-match+m-match))/(n+m);
		g= f-i;				// decimal fragment of the term above
		return ((int) (0.5+g*ED[i+1]+(1-g*ED[i]));
	else
*/
	return ((int) (0.5+1000.0*(n-match+ m-match)/(n+m)));

}
/**
 * compute lattice coordinate of vertex.
 * <main>
 * Creation date: (20.01.2002 17:34:38)
 * @param v align.Vertex vertex of which the coordinates are to be saved in p-Array
 * @param p int[] Array which stores the (tail) coordinates of edges
 */
protected void coord(Vertex v, int[] p) {

	int	i;
	CoordinateValues a;	// *a (!)
	Coordinate b;		// *b (!)

	if (v!=presource)	// presource (global) vertex before source, tail of first edge
		for (i=K, a=v.prev_coord_val; i>=1; i--, a=b.prev_coord_val) {
			b = a.curr_coord;
			p[i] = a.value;						// KICKED: p[i] = a - b.coord_vals + b.lo;
	  	}
	else				// if presource (vertex @ tail of first edge)
		for (i=K;i>=1;i--)
			p[i] = -1;

/* cost : compute cost of edge <q,r> preceded by edge <p,q> */
/*int cost (p, q, r)
int	p[], q[], r[];
{
	register	char	C [NUMBER+1];
	register	int	w, I, J, t [2] [NUMBER+1];

	for (I=1;I<=K;I++) {
		C[I] = (r[I]>q[I] ? S[I][r[I]] : DASH);
		t[0][I] = q[I] - p[I];
		t[1][I] = r[I] - q[I];
	}
	if(gflag)
          for (I=1;I<=K;I++)
            if (q[I]==0||q[I]==N[I])
              t[0][I]=t[1][I]=2;
	for (I=2,w=0;I<=K;I++)
          for (J=1;J<I;J++)
	    w+=scale[J][I]*(D[C[J]][C[I]]+T[t[0][I]][t[0][J]][t[1][I]][t[1][J]]);
	return w;
}*/

}
/**
 * create an coordinate.
 * variable counter deleted (does nothing -> program speedup 1/7)
 * Creation date: (20.01.2002 12:51:45)
 * @return align.Coordinate
 * @param index int[]
 * @param n int
 * @param prev_coord_val align.CoordinateValues
 */
protected Coordinate create_coordinate(int[] index, int n, CoordinateValues prev_coord_val) {
	int i, counter;					// loop counter for Pointer-Array
	Coordinate a;

	if (avail_coordinate!= null) {							// look for available COORDINATE on free list (RE-INCARNATION)
		a = avail_coordinate;
		avail_coordinate = a.next_on_free_list;
	} else
		a = new Coordinate();

	a.lo = index[0];
	a.hi = index[n-1];
	a.prev_coord_val = prev_coord_val;
	a.refer = n;

	a.coord_vals = new CoordinateValues[a.hi- a.lo +1];		// init Array
	for (counter= a.hi, i= a.hi-a.lo; i>=0; i--, counter--) {
		a.coord_vals[i]= new CoordinateValues();			// init Positions in Array: CoordinateValues(next_coord= null, curr_coord= null)
		a.coord_vals[i].value= counter;						// value of coordinate in absolute terms, needed for safe_coord(), coord()
	}
	for (n--;n>=0;n--) 										// link each new COORDINATE_VALUE to its parent COORDINATE and adjust
		a.coord_vals[index[n]- a.lo].curr_coord = a;    	// 	bounds; note that this loop works only on the possible values which
															// 	are stored in the array index
	return a;

}
/**
 * create edge v -> w.
 * Creation date: (20.01.2002 12:35:38)
 * <main>
 * @return align.Edge v --e--> w
 * @param v align.Vertex e.tail
 * @param w align.Vertex e.head
 */
protected Edge create_edge(Vertex v, Vertex w) {

	Edge e;		// the edge to be created

				// GET an Edge-Object
	if (avail_edge!=null) {		// re-incarnation
		e= avail_edge;
		avail_edge= e.next;		// abused as ll for avaiable edges
	} else
		e= new Edge();			// or new creation

				// CROSS-CHAIN e with its tail v
	e.tail = v;					// end-point of e
	e.next = v.out;				// add to double chained linked-list of v's outgoing edges
	v.out = e;					// e= new start of ll
	if (e.next!= null)			// if it has a successor, let him know
		e.next.prev = e;
	e.prev = null;				// e is new start of ll, so theres no precursor


    			// CROSS-CHAIN e with its head w
	e.head = w;
	if (e.head.nonextracted_inedges!= null)				// if theres already a linked list, chain it behind e
          e.head.nonextracted_inedges.nonextracted_prev = e;
	e.nonextracted_next = e.head.nonextracted_inedges;	// cross-chain in linked list, so that e knows its successor (or null)
	e.nonextracted_prev = null;							// e inserted @ start -> no precursor
	e.head.nonextracted_inedges = e;					// e= new start of linked list

				// PRE-INIT rest of attributes
    e.refer = 0;				// e isnt yet in an backtrack pathway; no backtracking edges are pointing to it
	e.backtrack = null;			// ...and itself isnt pointing to any
	e.heap_succ = 				// will be initialized by insert() in heap
	e.heap_pred = e;
	e.dist = Upper+1;   		// set e's distance to upper bound for lack of more information

/*	if (!already) {
		String tst= e.toString();
		if ((tst.indexOf(" 0")== (-1))&& (tst.indexOf("{0")== (-1))) {
			System.out.println("created edge: "+tst);
			already= true;
		}
	}
	if (already){
		String tst= e.toString();
		if ((tst.indexOf(" 0")!= (-1))|| (tst.indexOf("{0")!= (-1))) {
			System.out.println("falling back: "+tst);
			already= false;
		}
	}
*/
	return e;

}
/**
 * Insert the method's description here.
 * Creation date: (20.01.2002 12:15:05)
 * @return align.Vertex
 */
protected Vertex create_vertex(CoordinateValues prev_coord_val) {

	Vertex v;
	if (avail_vertex!= null) {		// reincarnation
		v= avail_vertex;
		avail_vertex = v.out.head;	// (Edge leading to) next avaiable Vertex
	} else
        v= new Vertex();

	v.prev_coord_val= prev_coord_val;	// init attributes
	v.out= null;
    v.nonextracted_inedges= null;

	return v;
}
/**
 * Die Beschreibung der Methode hier eingeben.
 * Erstellungsdatum: (18.01.2002 15:23:13)
 * @param a char
 */
protected short dag(char a) throws CancelException
{

//	byte a_= (byte) (((byte) a)- (byte) '-');
	return sub(a, a);
//	return D[(byte)a][(byte)a];
}
/**
 * Die Beschreibung der Methode hier eingeben.
 * Erstellungsdatum: (18.01.2002 15:23:13)
 * @param a char
 */
protected void dag(char a, short newValue) {

//	byte a_= (byte) (((byte) a)- (byte) '-');	// 'A'= 65; '-'= 45; (!)
	sub(a, a, newValue);
//	D[(byte)a][(byte)a]= newValue;
}
/**
 * read input data.
 * (comprises getseq())
 * Erstellungsdatum: (18.01.2002 13:55:53)
 * @param stream java.io.FileInputStream
 * @param Fname char[]
 */
protected void data(BufferedReader stream, String fName) {

	int i, j;
	int n, nb;
	char 	a, b;
	int 	c;
	String p;
	String s;
	BufferedReader fp= null;
	char[][] buffer= new char[LENGTH+1][]; // protected in C is super-local !!

										// read sequences
	String[] seqsStr= new String[NUMBER];

							/* read sequences */
	try {
		fp = new BufferedReader(new FileReader(fName));
	} catch (FileNotFoundException e) {
		System.err.println(e);
	}
	if (fp == null)
		fatal("Cannot open input file");

	p= "";
	nb= -1;
	try{
	while (fp.ready()) { /*  read in next line */
		s = fgets(MAXLINE, fp).trim();
		if ((s.length()>= 1) && (s.charAt(0)== '>')) {
			if (nb>= 0) {
				seqsStr[nb]= p;		// save sequence
				p= "";
			}
			nb++;					// increase number
		} else
			p+= s;

	}
	seqsStr[nb]= p;
	fp.close();
	}catch (IOException e) {
		System.err.println(e);
	}

	if ((nb+1)> NUMBER) {
		System.out.println((nb+1)+" is too many seqs");
		return;					// System.exit(-1);
	}
	for (K= 1; K<= nb+1; K++){
		buffer[K]= new char[seqsStr[K-1].length()+1];
		for (n= 0, i= 0; i< seqsStr[K-1].length(); i++) {
			a= seqsStr[K-1].charAt(i);
			if (n> LENGTH+1) {
				System.out.println("seq "+K+" too long (MSA)");
				return;			// System.exit(-1);
			}
			if ((a>= 'A')&& (a<= 'Z')) {	// cut wrong characters
				buffer[K][i]= a;			// convert to char[]
				n++;
			}
		}
		N[K]= n;							// save in length array
		buffer[K][i]= '\0';					// append '\0' to end of buffer
		S[K]= new char[n+2];				// initial DASH and terminal '\0' ??
		S[K][0]= getDASH();						// append '-' to start of seq
		for (i= 0; i< n+1; i++)
			S[K][i+1]= buffer[K][i];
	}
	K--;
/*
	for (i= 0; i< seqsStr.length; i++) {	// test init
		for (j= 0; j< S[i].length; j++)
			System.out.print(S[i][j]);
		System.out.print(" : "+N[i]);
		System.out.println();
	}
*/

					// *** 	INIT GAP COST,
					//		POSITIVE SYMMETRIC INTEGER DISTANCE TABLES,
					//		GAP COUNT TABLE								***

	if (stream!= null) {			// if(a cost input file was provided)

		s= "";
		StringTokenizer st;
		try{
			while (s.equals("")&& stream.ready())				// skip leading empty lines
				s= stream.readLine().trim();
			G= Integer.parseInt(s);	// gap opening costs
			s= "";
			while (stream.ready()) {						// read & parse linewise throug file
				s = stream.readLine().trim();
				while (s.equals("")&& stream.ready())			// skip empty lines
					s = stream.readLine().trim();
						// bug fix
				if (!s.trim().equals("")) {
					st= new StringTokenizer(s, " ", false);

					a= st.nextToken().charAt(0);				// triples of the form
					b= st.nextToken().charAt(0);				// <a	b 	d(a,b)>
					c= Integer.parseInt(st.nextToken());

					D[a][b]	= D[b][a] =	(short) c;				// distance table
				}
			}
			stream.close();
		} catch(IOException e) {
			System.err.println(e);
		}

	} else {						// else(load default dayhoff matrix - blosum62)

		G= 8;							// gap (opening) cost
		dag('-', (short) 0);			// tuples: distance of every two equal characters
		dag('A', (short) 0);
		dag('G', (short) 0);
		dag('T', (short) 0);
		dag('C', (short) 0);
		sub('-','A', (short) 2);		// triples: distance of each pair of different chars
		sub('-','G', (short) 2);
		sub('-','T', (short) 2);
		sub('-','C', (short) 2);
		sub('A','G', (short) 1);
		sub('A','T', (short) 2);
		sub('A','C', (short) 2);
		sub('G','T', (short) 2);
		sub('G','C', (short) 2);
		sub('T','C', (short) 1);
	}

	GG= (gflag!= 0) ? 0 : G;			// end gaps not penalized or penalized

	// T[w][x][y][z]            wy xz

	// T[i][x][j][y]	| (i,j ? seq1) ^ (x,y ? seq2)



	T[0][0][0][0]= 0;					// -- : --
	T[0][0][0][1]= G;					// -- : -x
	T[0][0][1][0]= G;					// -x : --
	T[0][0][1][1]= 0;					// -x : -x
	T[0][1][0][0]= 0;					// -- : x-
	T[0][1][0][1]= 0;					// -- : xx
	T[0][1][1][0]= G;					// -x : x-
	T[0][1][1][1]= 0;					// -x : xx
	T[1][0][0][0]= 0;					// x- : --
	T[1][0][0][1]= G;					// x- : -x
	T[1][0][1][0]= 0;					// xx : --	?
	T[1][0][1][1]= 0;					// xx : -x	?
	T[1][1][0][0]= 0;					// x- : x-	?
	T[1][1][0][1]= G;					// x- : xx	?
	T[1][1][1][0]= G;					// xx : x-
	T[1][1][1][1]= 0;					// xx : xx
	T[2][0][2][0]= 0;					// -- : --
	T[2][0][2][1]= 0;					// -- : --
	T[2][1][2][0]= 0;					// -- : --
	T[2][1][2][1]= 0;					// -- : --
	T[0][2][0][2]= 0;					// -- : --
	T[0][2][1][2]= 0;					// -- : --
	T[1][2][0][2]= 0;					// -- : --
	T[1][2][1][2]= 0;					// -- : --
	T[2][2][2][2]= 0;					// -- : --






}
/**
 * heap deletion
 * Erstellungsdatum: (18.01.2002 17:23:00)
 * @param e align.Edge
 * @param h align.Heap
 */
protected void delete(Edge e, Heap h) {

	Edge b= h.bucket[e.dist];					// get bucket of edges from heap for this distance

	if (e.heap_pred!= null)						// if node to delete (=e) has a 'previous'
		e.heap_pred.heap_succ= e.heap_succ;			// chain 'previous' with 'next'
	else {										// else it is 'root of bucket'
		b= e.heap_succ;								// set 'next'/null as 'root of bucket' (=b)
		h.bucket[e.dist]= b;
	}
	if (e.heap_succ!= null)						// if e has a 'next Edge'
		e.heap_succ.heap_pred= e.heap_pred;			// set 'previous Edge' of this to 'previous' of e/null

}
/**
 * display multiple sequence alignment.
 * <main>
 * Creation date: (20.01.2002 17:40:17)
 * @param e align.Edge
 */
protected void display(Edge e) throws CancelException {

	int[] p= new int[NUMBER+1];
	int[] q= new int[NUMBER+1];
	int[] r= new int[NUMBER+1];
	int d, I, J;
	Edge f;		// *f (!)
	d= 0;		// dummy for compiler

			// shortest sink to source path in lattice within bound?
	if (e==null || (d=e.dist)>Upper) {	// heap ran empty || lastMan not good enough
		fatal("Multiple alignment within bound does not exist.");
		return;
	}

			// recover shortest path to source tracing backward from sink
	for (e.next=null;e.tail!=presource;e=f) {
        f=e.backtrack;
		coord(f.tail,p);
		safe_coord(f.head,q);
		safe_coord(e.head,r);
		f.next = e;
		project(p,q,r);
	}

			// display alignment corresponding to path
	System.out.println("\n                  ***  Optimal Multiple Alignment  ***\n");
	for (e=e.next; e!=null; e=e.next) {
		coord(e.tail,p);
		safe_coord(e.head,q);
		column(p,q);
	}
	output();

			// output statistics
	if (gflag!= 0)
          System.out.print( "End gaps not penalized.");
	System.out.println("\nCostfile:                   "+cname);
	System.out.println("Alignment cost: "+d+"     Lower bound: "+Lower);
	System.out.println("Delta:          "+(d-Lower)+"     Max. Delta:  "+delta);
	System.out.println();
	System.out.print("Sequences \tProj. Cost \tPair. Cost \tEpsilon ");
	alicost= d;
	if (bflag!= 0) {
		System.out.println("\tMax. Epsi. \tWeight \tWeight*Cost");
		for (I=1;I<K;I++)
                  for (J=I+1;J<=K;J++)
		    System.out.println(
			    I+"\t"+J+"\t\t"+proj[I][J]+"\t"+costs[I][J]+"\t\t"+(proj[I][J]-costs[I][J])+"\t\t"+
			    epsi[I][J]+"\t\t"+scale[I][J]+"\t"+scale[I][J]*proj[I][J]);
	}
	else {
		System.out.println("Max. Epsi.");
		for (I=1;I<K;I++)
                  for (J=I+1;J<=K;J++)
		    System.out.println(I+" \t"+J+" \t\t"+proj[I][J]+" \t\t"+
			    costs[I][J]+" \t\t"+(proj[I][J]-costs[I][J])+" \t\t"+epsi[I][J]);
	}

}
/**
 * computes the average distance between the nodes pointed
 * to by node_index[i] and node_index[j]
 * <bias.c>
 * Creation date: (19.01.2002 21:44:45)
 * @return float
 * @param i int indices to node_index
 * @param j int
 */
protected float dist(int i, int j) {

	count=1;
	int rdist= rdist(node_index[i], node_index[j]);
			// this is definitely false !!!
			// the result is cast to int and then back to float !
			// (but they did it in the c-code this way, I checked by console output)
	float fff= (float) (rdist/count);
//	System.out.println("dist "+i+","+j+"="+rdist+"/"+count+"="+fff);
	return fff;
}
/**
 * ***************************************************************************** <br>
 * MSA subroutines.<br>
 * Version 1.0   June 27, 1989<br>
 * <br>
 * Program by SF Altschul<br>
 * <br>
 * This program estimates epsilons for all pairs of input sequences.
 * It does so by constructing an heuristic multiple alignment, and
 * calculating the difference between the imposed pairwise costs and
 * the optimal pairwise costs.  The heuristic involves finding all
 * positions consistent among all optimal pairwise alignments, and
 * forcing them into alignment.  For intervening segments, a progessive
 * alignment approach is taken.<br>
 * <br>
 * ****************************************************************************** <br>
 * <ecalc>
 * @globalReturn epsi[][] Epsilon values
 * Creation date: (20.01.2002 00:15:24)
 * @param len int
 * @throws Exception (passing from newal() via align())
 */
protected int ecalc(int len) throws Exception {

	int	I,J,i,j,c,pos,p1,oldp1,test;

	for (I=1;I<=K;++I)
		for (J=1;J<=K;++J)
			Con[I][J][N[I]+1]=N[J]+1;	// init match after last chars in two seqs
	for (I=1;I<=K;++I)
		AL[I][0]= (byte) getDASH();					// init AL in area of seqs with '-'

	pos=oldp1=0;
	for (p1=1;p1<=N[1]+1;++p1) {		// for all positions + LastMatch of the first seq
			// Search for positions consistent among all pairwise alignments and
			// align segments of all sequences from last consistent position to
			// current one.
		for (I=test=2;test!=0 && I<=K;++I) {	// test 1 against rest I seqs on matching; (test= 2) to start loop with (test!= 0);
			test=((i=Con[1][I][p1])>=0)?1:0;	// i= pos in seq-I aligned with pos p1@seq-1; (bool->int)
			for (J=I+1;test!=0 && J<=K;++J)		// if seq-I did match with seq-1 check all other matches seq-J with seq-1 && seq-I with seq-J
				test=((j=Con[1][J][p1])>=0 && Con[I][J][i]==j)?1:0;	// (bool -> int)
		}
		if (test!=0) {
			//	Pick an order in which to construct a progressive multiple alignment
			order(oldp1,p1);

			//	Align all segments in the order chosen
			j=align(oldp1,p1,len); // aligns segments and returns length of layout

			//	Add the aligned segments to the heuristic alignment (after pos)
			for (I=1;I<=K;++I) {
				J=Ord[I];
				for (i=1;i<=j;++i)
					AL[J][pos+i]= (byte) S[J][NAL[I][j-i]]; // translate pos's from NAL[][] to chars
			}

			//	Record whether a given alignment position is consistent among all
			//	optimal alignments, or has arisen by a progressive alignment
			for (i=1;i<=j;++i) {
				 						// write in track 0 of AL
				AL[0][pos+i]= (byte) (i>K?' ':
					// bug fix by micha:
						// Num can only display 36 digits (0-9, A-Z)
						// catch ArrayOutOfBounds for Num with the next line
					((Ord[i]> (Num.length-1))? 'x':
					Num[Ord[i]]));	// print seq order for this segment and fill with spaces
			}
			AL[0][pos+=i]='*';						// border has to be consistent among all (prints over seq-order in segments< K)
			for (I=1;I<=K;++I)
				AL[I][pos]= (byte) S[I][Con[1][I][p1]];		// insert border consistent characters

			oldp1=p1;		// next loop
		} // end: if (consistent position found)
	} // end for (all chars of seq1 + e)


			// Print out the heuristic multiple alignment
	if (oflag!= 0) {
		System.out.println("\n                 ***  Heuristic Multiple Alignment  ***\n\n");
		for (i=(--pos);i>0;i-=LINE) {
			for (I=0;I<=K;++I) {
				for (j=1;j<=LINE&&j<=i;++j)
					System.out.print((char)AL[I][pos-i+j]);

				System.out.println();
			}
			System.out.println();
		}
	}

			// For each pair, calculate difference between imposed alignment cost
			// and optimal alignment cost
	for (I=1;I<K;++I)
		for (J=I+1;J<=K;++J) {
			c= pcost(I,J,pos)- costs[I][J];				// real difference between real and optimal aligning path in projection
			epsi[I][J]= (c<MINE)?MINE:(c>MAXE?MAXE:c);	// Epsilon cut to given borders
		}

	return pos;
}
/**
 * extract	minimum	distance edge from heap.
 * <main>
 * Erstellungsdatum: (21.01.2002 11:02:18)
 * @return align.Edge
 */
protected Edge extract() {

	Edge b,e;		// **b, *e;
	int i;

	b= null;		// dummy for compiler

	for	(i= h.min; i<=h.max; i++)		// find first bucket in heap
		if ((b=h.bucket[i])!=	null)
			break;

	if ((h.min= i) > h.max)				// end: no buckets in heap
		  return null;

	e =	b;								// *b

	delete(e,h);
													// remove e from the list of non-extracted	in-edges incident to
													// the	head of	e
	if (e.nonextracted_prev !=	null)
		  e.nonextracted_prev.nonextracted_next =	e.nonextracted_next;
	else
		  e.head.nonextracted_inedges	= e.nonextracted_next;
	if (e.nonextracted_next !=	null)
		  e.nonextracted_next.nonextracted_prev =	e.nonextracted_prev;
	return e;

}
/**
 * ***************************************************************************** <br>
 * MSA subroutines.<br>
 * Version 1.0   June 27, 1989<br>
 * <br>
 * Program by JD Kececioglu & SF Altschul<br>
 * <br>
 * This program calculates pairwise costs using the current cost file.
 * Certain alignment positions may be forced by calling fix().
 * Terminal gaps may be counted differently than internal gaps.<br>
 * <br>
 * ****************************************************************************** <br>
 * <faces>
 * Creation date: (20.01.2002 01:25:19)
 */
protected void faces() throws CancelException {

	char[]	A,B;					// sequences (char[]) to be compared
	//	could be super-global, but only 1 func in <faces> !!!
	int I, J;						// seq# of seqs to be compared
	int	n, m;						// length of seqs to be compared
	int	i, j;						// loop counter for pos within seq I/J
	int Gi, Gj;						// act. gapcost for pos i/j (terminal or normal)
	int q;							// forced position flag for pos. i (C[i])
	int[]	C= new int[LENGTH+1];	// Forced alignment positions

	int U;							// upper bound for a given seq-pair
										// (by adding epsilon for that seq pair to the cost
										// of optimal pw alignment)
	// int Lower 			(global): 	sum of all weighted optimal alignment costs

	int[]	col= new int[LENGTH+1];		// column over all values which fill Carillo-Lipman
	int w;								// pointer to max pos. in col

	// int[][] dd, hh, vv	(global):	re-used from primer for calculating forward matrix
	int[]	d_= new int[LENGTH+1];      // reverse diagonal   distance
    int[]   h_= new int[LENGTH+1];	    // reverse horizontal distance
	int[]	v_= new int[LENGTH+1];	    // reverse vertical   distance
	int h_l, h_r, 						// save pos's from last calculated reverse matrices column/seq
		d_l, d_r,
		v_l, v_r;
	Row	row;							// a Row saved in faces, holding all Carillo-Lipman tested coordinates




	if (oflag!=0) {						// output
		for (i=(K*K-K)/2;i!=0;--i) 			// for (all combis k(k-1)/2)
			System.err.print(".");
		System.err.println();
	}


	if (fflag==0) 						// no fixed positions
		for (i=1;i<=LENGTH;++i)
			C[i]=0;							// init fixed pos array



	for (Lower=0, I=1;I<K;I++) 			// for all combinations
		for (n=N[I],A=S[I], J=I+1;J<=K;J++) {

			if (fflag!=0)
				fix(I,J,C,n);				// force fixed positions

			m = N[J];						// load data of 2nd seq
			B = S[J];


				// FORWARD SIM-TABLE: compute distance from <0,0> to <i,j>

			dd[0][0] = 0; hh[0][0] = vv[0][0] = GG;		// init corner
			for (j=1;j<=m;j++) {						// init border m (J,B)
				vv[0][j] = dd[0][j] = BIG;				// coming from nirvana
				hh[0][j] = hh[0][j-1] + sub(getDASH(), B[j]);// sum up '-' costs; + D[DASH][B[j]];
			}
			for (i=1;i<=n;i++) {						// init border n (I,A)
				hh[i][0] = dd[i][0] = BIG;
				vv[i][0] = vv[i-1][0] + 				// add up '-' costs OR
						(C[i]!=0 ? BIG : sub(A[i],getDASH()));// forced positions: are infinitly expensive; D[A[i]][DASH]
			}

			for (i=1;i<=n;i++) 							// compute rest of matrix
				for (q=C[i], Gi= (i==n) ? GG:G,
						j=1;j<=m;j++) {
					Gj= (j==m) ? GG : G;

					dd[i][j] = min3(dd[i-1][j-1],hh[i-1][j-1],vv[i-1][j-1])
									+ (q!=0 && q!=j ? BIG : sub(A[i],B[j]));	// D[A[i]][B[j]]
					hh[i][j] = min3(dd[i][j-1]+Gi,hh[i][j-1],vv[i][j-1]+Gi)
									+ sub(getDASH(),B[j]);							// D[DASH][B[j]]
					vv[i][j] = min3(dd[i-1][j]+Gj,hh[i-1][j]+Gj,vv[i-1][j])
									+ (q!=0 ? BIG : sub(A[i],getDASH()));			// D[A[i]][DASH]
    			}

			U= (costs[I][J]= min3(dd[n][m],hh[n][m],vv[n][m])) 	// store optimal pairwise alignment score in costs[I][J]
				+ epsi[I][J];									// add epsilon for seq pair to determine upper bound for seq pair
			Lower += scale[I][J] * costs[I][J];


				// BACKWARD SIM-TABLE: compute distance from <n,m> to <i,j>

			d_[m] = 0; h_[m] = v_[m] = GG;				// init corner
			for (j = m-1; j >= 0; j--)					// init border 1
        		v_[j]= (d_[j]= h_[j]= h_[j+1] + sub(getDASH(),B[j+1])) + G;	// D[DASH][B[j+1]]
  	  		for (w= j=0;j<m;j++)										// add up border columns of matrices
				if (min3(hh[n][j]-GG,dd[n][j],vv[n][j])+h_[j]<= U) {	// Carillo-Lipman Heuristic
					col[w++] = j;
//					dist[w++] = h_[j] - GG;
				}
    		col[w++] = m;												// m-th pos always completes CarLip ??!
//			dist[w++] = 0;

			row= face[I][J][n]= new Row();				// save first Row in face[seq1][seq2][pos@seq1]: pos's in seq 2
			row.width = w;									// row.length
	    	row.column = vector(col,w);						// copy possible values of seq2 from col-container

            for (i=n-1;i>=0;i--) {						// calculate rest of backward matrix (col-wise) and add up with forward matrix
				q=C[i+1];									// flag for forced pos
				Gi= (i==0) ? GG : G;						// normal / terminal gap costs

				h_r = h_[m]; d_r = d_[m]; v_r = v_[m];		// save last calculated row ??!
				h_[m] = (d_[m] = v_[m] = v_r + sub(A[i+1],getDASH())) + G;	// D[A[i+1]][DASH]
				for ( j = m-1; j >= 0; j-- ) {
	    			Gj= (j==0) ? GG : G;

	    			h_l = h_[j]; d_l = d_[j]; v_l = v_[j];	// save last calculated row ??!
	    			d_[j]= min3(d_r, h_r, v_r) +
									(q!=0 && q!=j+1 ? BIG : sub(A[i+1],B[j+1]));	// D[A[i+1]][B[j+1]]
					h_[j]= min3(d_[j+1]+Gi,h_[j+1],v_[j+1]+Gi) +
									sub(getDASH(),B[j+1]);								// D[DASH][B[j+1]]
	    			v_[j]= min3(d_l + Gj, h_l + Gj, v_l) +
									(q!=0 ? BIG : sub(A[i+1],getDASH()));				// D[A[i+1]][DASH]
	    			h_r = h_l; d_r = d_l; v_r = v_l;
				}

	       		for (w= j=0;j<=m;j++) {				// add up FORWARD and reverse MATRIX
					Gj= ((j==0) || (j==m)) ? GG : G;	// terminal and end gap
					if (min3(hh[i][j]+min3(h_[j]-Gi,d_[j],v_[j]),
							dd[i][j]+min3(h_[j],d_[j],v_[j]),				// add up matrices and determine min
							vv[i][j]+min3(h_[j],d_[j],v_[j]-Gj) ) <=U) {	// Carillo-Lipman
						col[w++]=j;
//						dist[w++]=min3(h_[j]-Gi,d_[j],v_[j]-Gj);
					}
				}

				row = face[I][J][i]= new Row(); 	// save Row in face[seq1][seq2][pos@seq1]: Carillo-Lipman-tested positions in seq2
				row.width = w;							// column.length
				row.column= vector(col,w);				// copy positions from col (container)
			} // for (all pos. in seq I)

			if (oflag!=0)
				System.err.print("*");				// fill '.' with '*' (one for each tested combination)

		}	// for(all combinations of seqs)


	if (oflag!=0)
		System.err.println();
}


/**
 * ***************************************************************************** <br>
 * MSA subroutines.<br>
 * Version 1.0   June 27, 1989<br>
 * <br>
 * Program by JD Kececioglu & SF Altschul<br>
 * <br>
 * This program calculates pairwise costs using the current cost file.
 * Certain alignment positions may be forced by calling fix().
 * Terminal gaps may be counted differently than internal gaps.<br>
 * <br>
 * ****************************************************************************** <br>
 * <faces>
 * Creation date: (20.01.2002 01:25:19)
 */
protected void faces_check() throws CancelException {

	int lastColMin= -1;

	char[]	A,B;					// sequences (char[]) to be compared
	//	could be super-global, but only 1 func in <faces> !!!
	int I, J;						// seq# of seqs to be compared
	int	n, m;						// length of seqs to be compared
	int	i, j;						// loop counter for pos within seq I/J
	int Gi, Gj;						// act. gapcost for pos i/j (terminal or normal)
	int q;							// forced position flag for pos. i (C[i])
	int[]	C= new int[LENGTH+1];	// Forced alignment positions

	int U;							// upper bound for a given seq-pair
										// (by adding epsilon for that seq pair to the cost
										// of optimal pw alignment)
	// int Lower 			(global): 	sum of all weighted optimal alignment costs

	int[]	col= new int[LENGTH+1];		// column over all values which fill Carillo-Lipman
	int w;								// pointer to max pos. in col

	// int[][] dd, hh, vv	(global):	re-used from primer for calculating forward matrix
	int[]	d_= new int[LENGTH+1];      // reverse diagonal   distance
    int[]   h_= new int[LENGTH+1];	    // reverse horizontal distance
	int[]	v_= new int[LENGTH+1];	    // reverse vertical   distance
	int h_l, h_r, 						// save pos's from last calculated reverse matrices column/seq
		d_l, d_r,
		v_l, v_r;
	Row	row;							// a Row saved in faces, holding all Carillo-Lipman tested coordinates




	if (oflag!=0) {						// output
		for (i=(K*K-K)/2;i!=0;--i) 			// for (all combis k(k-1)/2)
			System.err.print(".");
		System.err.println();
	}


	if (fflag==0) 						// no fixed positions
		for (i=1;i<=LENGTH;++i)
			C[i]=0;							// init fixed pos array



	for (Lower=0, I=1;I<K;I++) 			// for all combinations
		for (n=N[I],A=S[I], J=I+1;J<=K;J++) {

			if (fflag!=0)
				fix(I,J,C,n);				// force fixed positions

			m = N[J];						// load data of 2nd seq
			B = S[J];


				// FORWARD SIM-TABLE: compute distance from <0,0> to <i,j>

			dd[0][0] = 0; hh[0][0] = vv[0][0] = GG;		// init corner
			for (j=1;j<=m;j++) {						// init border m (J,B)
				vv[0][j] = dd[0][j] = BIG;				// coming from nirvana
				hh[0][j] = hh[0][j-1] + sub(getDASH(), B[j]);// sum up '-' costs; + D[DASH][B[j]];
			}
			for (i=1;i<=n;i++) {						// init border n (I,A)
				hh[i][0] = dd[i][0] = BIG;
				vv[i][0] = vv[i-1][0] + 				// add up '-' costs OR
						(C[i]!=0 ? BIG : sub(A[i],getDASH()));// forced positions: are infinitly expensive; D[A[i]][DASH]
			}

			for (i=1;i<=n;i++) 							// compute rest of matrix
				for (q=C[i], Gi= (i==n) ? GG:G,
						j=1;j<=m;j++) {
					Gj= (j==m) ? GG : G;

					dd[i][j] = min3(dd[i-1][j-1],hh[i-1][j-1],vv[i-1][j-1])
									+ (q!=0 && q!=j ? BIG : sub(A[i],B[j]));	// D[A[i]][B[j]]
					hh[i][j] = min3(dd[i][j-1]+Gi,hh[i][j-1],vv[i][j-1]+Gi)
									+ sub(getDASH(),B[j]);							// D[DASH][B[j]]
					vv[i][j] = min3(dd[i-1][j]+Gj,hh[i-1][j]+Gj,vv[i-1][j])
									+ (q!=0 ? BIG : sub(A[i],getDASH()));			// D[A[i]][DASH]
    			}

			U= (costs[I][J]= min3(dd[n][m],hh[n][m],vv[n][m])) 	// store optimal pairwise alignment score in costs[I][J]
				+ epsi[I][J];									// add epsilon for seq pair to determine upper bound for seq pair
			Lower += scale[I][J] * costs[I][J];


				// BACKWARD SIM-TABLE: compute distance from <n,m> to <i,j>

			d_[m] = 0; h_[m] = v_[m] = GG;				// init corner
			for (j = m-1; j >= 0; j--)					// init border 1
        		v_[j]= (d_[j]= h_[j]= h_[j+1] + sub(getDASH(),B[j+1])) + G;	// D[DASH][B[j+1]]
  	  		for (w= j=0;j<m;j++)										// add up border columns of matrices
				if (min3(hh[n][j]-GG,dd[n][j],vv[n][j])+h_[j]<= U) {	// Carillo-Lipman Heuristic
					col[w++] = j;
//					dist[w++] = h_[j] - GG;
				}
    		col[w++] = m;												// m-th pos always completes CarLip ??!
//			dist[w++] = 0;

			row= face[I][J][n]= new Row();				// save first Row in face[seq1][seq2][pos@seq1]: pos's in seq 2
			row.width = w;									// row.length
	    	row.column = vector(col,w);						// copy possible values of seq2 from col-container

            for (i=n-1;i>=0;i--) {						// calculate rest of backward matrix (col-wise) and add up with forward matrix
				q=C[i+1];									// flag for forced pos
				Gi= (i==0) ? GG : G;						// normal / terminal gap costs

				h_r = h_[m]; d_r = d_[m]; v_r = v_[m];		// save last calculated row ??!
				h_[m] = (d_[m] = v_[m] = v_r + sub(A[i+1],getDASH())) + G;	// D[A[i+1]][DASH]
				for ( j = m-1; j >= 0; j-- ) {
	    			Gj= (j==0) ? GG : G;

	    			h_l = h_[j]; d_l = d_[j]; v_l = v_[j];	// save last calculated row ??!
	    			d_[j]= min3(d_r, h_r, v_r) +
									(q!=0 && q!=j+1 ? BIG : sub(A[i+1],B[j+1]));	// D[A[i+1]][B[j+1]]
					h_[j]= min3(d_[j+1]+Gi,h_[j+1],v_[j+1]+Gi) +
									sub(getDASH(),B[j+1]);								// D[DASH][B[j+1]]
	    			v_[j]= min3(d_l + Gj, h_l + Gj, v_l) +
									(q!=0 ? BIG : sub(A[i+1],getDASH()));				// D[A[i+1]][DASH]
	    			h_r = h_l; d_r = d_l; v_r = v_l;
				}

				int colMin= Integer.MAX_VALUE;
				int[] total= new int[m+1];
	       		for (w= j=0;j<=m;j++) {				// add up FORWARD and reverse MATRIX
					Gj= ((j==0) || (j==m)) ? GG : G;	// terminal and end gap

					total[j]=	min3(hh[i][j]+min3(h_[j]-Gi,d_[j],v_[j]),
							dd[i][j]+min3(h_[j],d_[j],v_[j]),				// add up matrices and determine min
							vv[i][j]+min3(h_[j],d_[j],v_[j]-Gj) );

						// check
					if (total[j]< colMin)
						colMin= total[j];

					if (total[j] <=U) {	// Carillo-Lipman
						col[w++]=j;
//						dist[w++]=min3(h_[j]-Gi,d_[j],v_[j]-Gj);
					}
				}
				boolean Ofound= false;
				for (j= 0; j<= m; ++j)
		       		if ((total[j]- colMin)== 0)
		       			Ofound= true;
		       	if (!Ofound)
		       		System.err.println("no 0-coordinate in row in MSA.faces: "+colMin);
				if ((lastColMin!= -1)&& (colMin!= lastColMin))
					System.err.println("faces() col error "+lastColMin+"!="+colMin+" in ["+I+","+J+"] at "+i+"!");
				else
					;// System.err.println("["+I+","+J+"] at "+i+" ok!");
				lastColMin= colMin;


				row = face[I][J][i]= new Row(); 	// save Row in face[seq1][seq2][pos@seq1]: Carillo-Lipman-tested positions in seq2
				row.width = w;							// column.length
				row.column= vector(col,w);				// copy positions from col (container)
			} // for (all pos. in seq I)

			lastColMin= -1;

			if (oflag!=0)
				System.err.print("*");				// fill '.' with '*' (one for each tested combination)

		}	// for(all combinations of seqs)


	if (oflag!=0)
		System.err.println();
}
/**
 * Die Beschreibung der Methode hier eingeben.
 * Erstellungsdatum: (22.01.2002 15:47:52)
 * @param format char[]
 * @param message1 char[]
 * @param message2 char[]
 */
protected void fatal(char[] format, char[] message1, char[] message2) {}
/**
 * fatal error.
 * <main.c>
 * Erstellungsdatum: (22.01.2002 15:47:52)
 * @param format char[]
 * @param message1 char[]
 * @param message2 char[]
 */
protected void fatal(java.lang.String message) {
				// void fatal (format,	message1, message2)
				// char	*format, *message1,	*message2;

				// fprintf(stderr,format,message1,message2);
	if (oflag!= 0)
		System.err.println(message);
	if (proxy!= null)
		proxy.setOutput(message+"\n");
//	gettimeofday(&endtime, NULL);  /* Get end time */
//	deltatime =	(endtime.tv_sec	- starttime.tv_sec)
//		 + (endtime.tv_usec	- starttime.tv_usec)/1000000.0;
//	printf("Elapsed	time = %7.3f\n",deltatime);
	return;			// System.exit(1);
}
/**
 * Die Beschreibung der Methode hier eingeben.
 * Erstellungsdatum: (29.01.2002 11:00:10)
 * @return char[]
 * @param s char[]
 * @param size int
 * @param stream java.io.FileReader
 */
protected String fgets(int size, BufferedReader stream) {

	String s= null;
	// size unused yet
	try {
		s= stream.readLine();	// +"\n"+"\0"
	} catch (IOException e) {
		System.err.println(e);
	}

	return s;

}
/**
 * *** fill takes information from input file and derives pFIX structure ***
 * Erstellungsdatum: (31.01.2002 14:08:35)
 * @param sqnum int
 * @param sqnames int[]
 * @param sqstarts int[]
 * @param Len int
 */
protected void fill(int sqnum, int[] sqnames, int[] sqstarts, int Len) {

	int I, J, ni, nj;
	/*	char * malloc();  NEW */
	Block bp1, bp2; // BLOCK*

	for (I = sqnum - 1; I!= 0; --I)
		for (J = I - 1; J >= 0; --J) {
			bp1 = new Block();
			bp2 = new Block();
			bp1.i = bp2.j = sqstarts[I];
			bp1.j = bp2.i = sqstarts[J];
			bp1.len = bp2.len = Len;
			ni = sqnames[I];
			nj = sqnames[J];
			bp1.next = pFIX[ni][nj];
			bp2.next = pFIX[nj][ni];
			pFIX[ni][nj] = bp1;
			pFIX[nj][ni] = bp2;
		}
}
/**
 * fix uses pFIX structure to "fix" points of alignment
 * (from fixer.c)
 * fix is called from primer() and faces()
 * Creation date: (19.01.2002 20:19:04)
 * @param I int 1st seq: pos1
 * @param J int	2nd seq: pos2
 * @param C int[] array with saved fixed positions (C[pos1]= pos2)
 * @param cl int length of 1st seq
 */
protected void fix(int I, int J, int[] C, int cl) {

	int i,j,k;
	Block p;

	for (i=1;i<=cl;i++) 		// init array in the interesting length (of 1st seq)
		C[i]=0;
	for (p=pFIX[I][J];p!=null;p=p.next) { // get fixed pos list (blocks) between seq I and J
		i=p.i;					// fixed start pos of block in 1st seq
		j=p.j;					// fixed start pos of block in 2nd seq
		for (k=0;k<p.len;k++) 	// extrapolate block length
			C[i++]=j++;			// by increasing pos's in both seqs over len
	}

}
/**
 * Die Beschreibung der Methode hier eingeben.
 * Erstellungsdatum: (31.01.2002 14:01:19)
 * @param fname java.lang.String
 */
protected void fixread(String fname) {

	BufferedReader fp= null;
	String s; //  = new char[MAXLINE];
	StringTokenizer st;
	String p; // char*
	String fw; // = new char[10];
	int i, j, len, sqnum;
	int[] sqnames = new int[NUMBER];
	int[] sqstarts = new int[NUMBER];

	try {
		fp = new BufferedReader(new FileReader(fname));
	} catch (FileNotFoundException e) {
		System.err.println(e);
	}
	if (fp == null)
		fatal("Cannot open " + fname);
	for (i = 1; i <= NUMBER; ++i)
		for (j = 1; j <= NUMBER; ++j)
			pFIX[i][j] = null;

	while ((s = fgets(MAXLINE, fp)) != null) { /*  read in entire line */
		p = s;
		st = new java.util.StringTokenizer(s, " ", false);
		sqnum = 0;
		p = st.nextToken();
		fw = p.trim();
		while (!fw.equals("S")) { /*get Sq names*/
			p = st.nextToken();
			sqnames[sqnum++] = Integer.parseInt(fw);
			fw = p.trim();
		}

		do { /* parse string through to end */
			i = 0;
			p = st.nextToken();
			fw = p.trim();
			while (!fw.equals("L")) { /*get Sq starts*/
				p = st.nextToken();
				sqstarts[i++] = Integer.parseInt(fw);
				fw = p.trim();
			}

			/* check consistency of block file */

			if (i != sqnum)
				fatal("Block starts do not equal sequence number.");
			p = st.nextToken();
			len = Integer.parseInt(p.trim()); /* get Block length */
			fill(sqnum, sqnames, sqstarts, len);

		}
		while (st.hasMoreTokens()); /* parse string to end */

	} /*  end loop for each line of input  */

	try {
		fp.close();
	} catch (IOException e) {
		System.err.println(e);
	}

}
/**
 * free an coordinate.
 * <main>
 * Creation date: (20.01.2002 13:17:43)
 * @param a align.Coordinate which is indicated to be freed
 * @return success status of liberation
 * (1 for success, 0 for failure due to other dependencies)
 */
protected int free_coordinate(Coordinate a) {

	Coordinate b;

	if (-- a.refer <= 0) {
		b = a.prev_coord_val.curr_coord;
		a.prev_coord_val.curr_coord = null;
//		free((char *)a.coord_vals);			// Garbage collector ?!?
		a.next_on_free_list =  avail_coordinate;
		avail_coordinate =  a;
		if (b!= null)
		  free_coordinate(b);
		return TRUE;
	}
	else return FALSE;
}
/**
 * free an edge.
 * Creation date: (20.01.2002 12:44:08)
 * @param e align.Edge
 */
protected void free_edge(Edge e) {

	if (e.prev!=null)				// remove from Edge list
          e.prev.next = e.next;
	else
          e.tail.out = e.next;

	if (e.next!=null)
		e.next.prev = e.prev;
	if ((e.backtrack!= null) && (-- e.backtrack.refer == 0))
		free_edge(e.backtrack);
	if ((e.head.out== null) && (e.refer != -1))
		free_vertex(e.head);
	e.next = avail_edge;
	avail_edge = e;

}
/**
 * free a vertex.
 * Creation date: (20.01.2002 12:27:52)
 * @param v align.Vertex
 */
protected void free_vertex(Vertex v) {

	Coordinate a;
	Edge e;

	a= v.prev_coord_val.curr_coord;
	v.prev_coord_val.curr_coord= null;
	free_coordinate(a);
    			// remove all the remaining edges coming into v from the heap
	for(e=v.nonextracted_inedges; e!= null; e= e.nonextracted_next) {
		delete(e,h);
		e.refer = -1;   // kludge telling free_edge() not to call free_vertex()
		free_edge(e);
	}
	v.out= new Edge(avail_vertex);	// (EDGE *) avail_vertex;
	avail_vertex = v;

}
/**
 * Die Beschreibung der Methode hier eingeben.
 * Erstellungsdatum: (21.01.2002 10:57:03)
 * @return align.Heap
 * @param max int
 */
protected Heap heap(int max) {

	Heap h;		// *h
	int i;

	h =	new Heap(max);		//'(HEAP *) alloc(sizeof( HEAP	) +	max*sizeof(	EDGE * ))'
	h.min = 0;				// init current indices
	h.max = max+1;			// max* upcount
	for	(i=0; i<=max; i++)
		  h.bucket[i]= null;			// init heap array of edges with 'null'
	return h;

}
/**
 * initialize adjacent vertex generator.
 * <main>
 * Creation date: (20.01.2002 17:18:54)
 * @param e align.Edge
 * @param p int[]
 * @param q int[]
 * @deprecated commented out
 */
protected void init_adjacent(Edge e, int[] p, int[] q) {
/*
	int i;
	CoordinateValues a;		// *a (!)

	coord(e.tail,p);
	safe_coord(e.head,q);
	for (i=K, a=e.head.prev_coord_val; i>=1; i--, a=a.curr_coord.prev_coord_val) {
		A[i] = a.curr_coord;
		B[i] = (P[i] = q[i]) + 1;
	}
*/
}
/**
 * reads in a pairwise distance matrix, creates leaf nodes,
 * and initializes node_index to point at them.
 * <bias.c>
 * Creation date: (19.01.2002 22:42:13)
 */
protected void init_data() {

	int i,j,t;
	indexlen = K;	// init global variable

	for (t=0;t<K;t++)
		Dij[t]= new int[K];		// '(int *) calloc(K, sizeof(int))' = allocate mem for K objects of size int
	node_index= new Node[K]; 	// '(NODE **) malloc(K * sizeof(NODE *))' = same as above ?!

	for (i=0;i<K-1;++i)
		for (j=i+1;j<K;++j) {	// double loop: all combinations of seqs
			Dij[j][i]= Dij[i][j]= scale[j+1][i+1];	// load scale factors symmetrical
//			System.out.println("init_scale ("+(j+1)+","+(i+1)+"= "+scale[j+1][i+1]+" Dji "+Dij[j][i]+" Dij "+Dij[i][j]);
		}

	for (i= Node.vcount= 0;i<K;i++)  // init vcount; for all seqs
		node_index[i]= new Node(i,(float)0.0,null,null);	// init Node array with dummy null-Nodes

}
/**
 * heap insertion
 * Erstellungsdatum: (18.01.2002 16:40:45)
 * @param e align.Edge
 * @param h align.Heap
 */
protected void insert(Edge e, Heap h) {

	Edge b= h.bucket[e.dist];		// b= root of bucket of edges containing e

	if (b!= null)					// if there is a root/bucket
		b.heap_pred= e;					// insert e before actual root
	e.heap_succ= b;					// append b/null to e (as 'next')
	e.heap_pred= null;				// e is new 'root', so previous edge in bucket= null

//	b= e;							// change new root of chain
	h.bucket[e.dist]= e;			// write back new root (pointing to rest of former bucket)

}
/**
 * intersect regions on rows of faces.<br>
 * Finds the possible values for the "seqnum"th coordinate of point,
 * given point[1],point[2],...,point[seqnum-1].  The corresponding
 * bound values are copied to bound[].<br>
 * <main>
 * Creation date: (20.01.2002 17:03:00)
 * @return int
 * @param point int[]
 * @param seqnum int
 * @param possible_values int[]
 */
protected int intersect(int[] point, int seqnum, int[] possible_values) {

	int	i, j, k;
	int[] c;				// *c (!)
	int	J, m, n;
	Row	r;					// *r (!)

			// retrieve values that are consistent with the pairwise alignment of
			// sequenecs 1 and seqnum
	r = face[1][seqnum][point[1]];			// WAS: (Row[])+ point[1]	***
	c = r.column;
	m = r.width;
	for (i=0;i<m;i++) {
		possible_values[i] = c[i];
	}
	for (J=2;J<seqnum;J++) {
			// Get values that are consistent with the pairwise
			// alignement of sequences J and seqnum
		r = face[J][seqnum][point[J]];		// WAS: (Row[])+ point[J]	***
		c = r.column;
		n = r.width;
			// compute intersection of possible values array and
			// c array for sequence j, keeping result in possible
			// values
		for (i=j=k=0; i<m && j<n;)
			if (possible_values[i]<c[j])
				i++;
			else if (possible_values[i] > c[j])
				j++;
			else {
				possible_values[k++] = possible_values[i++];
				j++;
			}
		if ((m=k)== 0) break;				// WAS: ((m=k)<= 0)	***
	}
	return m;

}
/**
 * FLAGS:
 * -m 	suppresses computation of optimal multiple alignment
 * -g 	penalizes terminal gaps
 * -b 	sets all pairwise weights to 1 (does not compute tree)
 * -o 	suppresses status reports to stderr
 * -d<delta> 		user speified delta (uppre bound for total alignment cost)
 * -e <epsilon file> user specified epsilons for each sequence pair
 *					(does not compute heuristic alignment)
 * -c <cost file>	user specified cost file (default is pam 250)
 * -f <fixer file>	allows user to force residues in alignment
 *					(see comment in fixer.c for file format)
 * <input filename>	sequences to be aligned
 * Erstellungsdatum: (17.01.2002 23:20:37)
 * @param args java.lang.String[]
 */
public static void main(String[] argv) throws CancelException {

	long starttime, endtime;
	double deltatime;

	int i, j, len, size;
	String s;
	byte eflag= 1, mflag= 1; // better than boolean

	// file names for:	(max: new char[FILE_NAME_SIZE])
	String ename= null; // epsilons,		dummy
	String fname= null; // fix positions,	dummy
	String fName= null; // sequences, 		'must open'
	java.io.BufferedReader efile= null; //	dummy
	java.io.BufferedReader stream;

	int argc= argv.length + 1; // (+1) the cmd itself (ported C)

	starttime= System.currentTimeMillis();

	// *** PROCESS ARGUMENTS ***

	MSA msa= new MSA();
	msa.oflag= 1; // (def 1);0= suppress output
	msa.bflag= 1;
	// (def 1);0= unweighted summing (DO NOT calculates evolutionary tree)
	msa.gflag= 1; // (def 1);0= end gaps PENALIZED (charged the same as internal ones)
	msa.fflag= 0; // 1= file with forced positions provided
	msa.delta= -1;
	stream= null;

	if ((argc < 1) || (argc > MAX_ARGS))
		msa.fatal(msa.USAGE);

	for (int ppos= 0;
		--argc != 0;
		ppos++) // start @ first parameter with parsing,
		// leave last parameter as seq input filename

		if (argv[ppos].charAt(0) == '-') // if (starts with switch signal)

			if (argv[ppos].charAt(1) == 'm')
				mflag= 0; //
	else
		if (argv[ppos].charAt(1) == 'g')
			msa.gflag= 0; // PENALIZE end gaps (free shift off)
	else
		if (argv[ppos].charAt(1) == 'b')
			msa.bflag= 0; // DO NOT weight seqs (no evol. tree)
	else
		if (argv[ppos].charAt(1) == 'o')
			msa.oflag= 0; // SUPPRESS output
	else
		if (argv[ppos].charAt(1) == 'd')
				msa.delta= Integer.parseInt(// delta is provided
	argv[ppos].substring(2, argv[ppos].length()));
		else
			if (argv[ppos].charAt(1) == 'e') {
				++ppos; // next arg is filename
				--argc; // one loop less
				ename= argv[ppos]; // epsilon file name
				eflag= 0; // epsilons are provided
			} else
				if (argv[ppos].charAt(1) == 'c') {
					++ppos; // next arg is filename
					--argc; // one loop less
					msa.cname= argv[ppos]; // cost file name
					try { // open cost file
						stream= new BufferedReader(new FileReader(msa.cname));
					} catch (FileNotFoundException e) {
						msa.fatal(e + "\nCannot open " + fName + "."); // 'must_open()'
					}
				} else
					if (argv[ppos].charAt(1) == 'f') {
						++ppos; // next arg is filename
						--argc; // one loop less
						fname= argv[ppos]; // force fix file name
						msa.fflag= 1; // fixed pos's are forced
					} else
						msa.fatal(msa.USAGE); // '-' w/o any argument char

	else // else (no switch)

		fName= argv[ppos]; // name of sequences input file

	// for (check all parameter)

	// *** GET INPUT ***

	msa.data(stream, fName); // read in sequences

	if (msa.fflag != 0)
		msa.fixread(fname); // open & read forced Fix Positions

	if (eflag == 0) { // get provided Epsilon Values

		try { // open file
			efile= new BufferedReader(new FileReader(ename));
		} catch (FileNotFoundException e) {
			msa.fatal(e + "\nCannot open " + fName + "."); // 'must_open()'
		}

		s= "";
		try {
			while (efile.ready())
				s += efile.readLine() + " ";
			efile.close();
		} catch (IOException e) {
			msa.fatal(e + "\nCannot read " + ename + ".");
		}

		StringTokenizer st= new StringTokenizer(s, " ", false);
		for (i= 1; i < msa.K; ++i)
			for (j= i + 1; j <= msa.K; ++j) {
				msa.epsi[i][j]= Integer.parseInt(st.nextToken());
				if (msa.epsi[i][j] < 0)
					msa.fatal("Epsilon must be positive.");
			}

	} // if (epsilons provided)

	// *** ALLOCATE MEMORY ***

	for (len= i= 1; i <= msa.K; i++) // find longest sequence
		if (msa.N[i] > len)
			len= msa.N[i];
	++len; // take one more... (1- to 0-based)
	size= len * (1 + eflag);
	for (i= 0; i < LENGTH + 1; i++) { // init rows of 2D-Arrays
		msa.dd[i]= new int[size];
		msa.hh[i]= new int[size];
		msa.vv[i]= new int[size];
	}

	// *** COMPUTE MULTIPLE SEQUENCE ALIGNMENT ***

	if (msa.K == 2) // pairwise alignment
		msa.bflag= 0; // no need of weights
	if (msa.bflag == 0) // no sequence weights
		for (i= 1; i < msa.K; ++i)
			for (j= i + 1; j <= msa.K; ++j) // for all sequence combinations
				msa.scale[i][j]= 1; // init scale factor with one (no scaling)

	if ((msa.bflag != 0) || (eflag != 0)) { // if weight factors or epsilons needed
		if (msa.oflag != 0)
			System.err.println("Calculating pairwise alignments.");
		msa.primer(); // init pairwise alignments
	}

	if (msa.bflag != 0) { // weight factors needed
		if (msa.oflag != 0)
			System.err.println("Calculating weights.");
		msa.bias(); // ~beeinflussen
	}

	if (eflag != 0) { // epsilons needed
		if (msa.oflag != 0)
			System.err.println("Calculating epsilons.");
		try {
			msa.ecalc(2 * len - 2); // calculate epsilons
		} catch (Exception e) {
			System.err.println(e);
		}
		if (msa.oflag != 0) {
			System.err.println("----Estimated epsilons----");
			for (i= 1; i < msa.K; ++i)
				for (j= i + 1; j <= msa.K; ++j) // all combinations of sequences
					System.err.println("I = " + i + "  J = " + j + "  epsilon = " + msa.epsi[i][j]);
		}
	}

	if (mflag != 0) { // optimal multiple alignment required
		if (msa.oflag != 0)
			System.err.println("Calculating pairwise projection costs.");
		msa.faces(); // projections of k-space = faces
		// FREE vv, hh, dd = Garbage.collector?
		if (msa.delta < 0)
			for (msa.delta= 0, i= 1; i < msa.K; ++i)
				for (j= i + 1; j <= msa.K; ++j) // all combinations of sequences
					msa.delta += (msa.scale[i][j] * msa.epsi[i][j]);
		msa.Upper= msa.Lower + msa.delta; // determine upper bound of costs
		for (i= 0; i < msa.K; ++i)
			for (j= i + 1; j <= msa.K; ++j) {
				msa.proj[i][j]= 0; // init half of the proj[][] with 0
				msa.scale[j][i]= msa.scale[i][j]; // make symmetric for msa()
			}

		if (msa.oflag != 0)
			System.err.println("Calculating multiple alignment.");
		msa.display(msa.msa()); // display result by edge computed by msa()
	}

		// running time
	endtime= System.currentTimeMillis();
	deltatime= endtime - starttime;
	if (msa.oflag!= 0)
		System.out.println("Elapsed time = " + deltatime);

}
/**
 * creates nodes for the tree.
 * <bias.c>
 * Erstellungsdatum: (22.01.2002 15:42:05)
 * @return align.Node
 * @param sqn int sequence number or -1 for internal nodes
 * @param ltop float distance to parent node
 * @param lt align.Node pointers to child nodes
 * @param rt align.Node
 */
protected Node makenode(int sqn, float ltop, Node lt, Node rt) {

	Node ptr;

	ptr= new Node();
	Node.vlist[Node.vcount++] = ptr;
	ptr.ltop = ltop;
	ptr.lt = lt;
	ptr.rt = rt;
	ptr.sqn = sqn;
	if (sqn<0) {
		lt.bro = rt;
		rt.bro = lt;
		lt.par = rt.par = ptr;
	}
	return(ptr);

}
/**
 * three argument minimum.
 * (from msa.c)
 * Creation date: (19.01.2002 19:40:46)
 * @return int
 * @param x int
 * @param y int
 * @param z int
 */
protected int min3(int x, int y, int z) {
	if (x<y)
		if (x<z)
			return x;
		else
			return z;
	else
		if (y<z)
			return y;
		else
			return z;
}
/**
 * selects the optimum pair of nodes from amongst the nodes
 * currently pointed to by node_index to make neighbors,
 * joins them through an internal node, and repeats the
 * process until only two nodes remain.
 * <bias.c>
 * Creation date: (19.01.2002 22:36:25)
 */
protected void minimize_Sij() {

	int i, j, min_i=0, min_j=0;
	float tmp, min = BIG +1;

	for (i=0; i<indexlen; i++)
		for (j=i+1; j<indexlen; j++) {	// for all combinations of Nodes
			tmp = compute_S(i,j);		// compute distance
 	 		if (tmp < min) {			// determine Minimum
	     		min_i= i;
				min_j= j;
				min= tmp;
			}
        }

	coalesce(min_i,min_j);				// create new internal node joining i and j
//	System.out.println("coalescing:"+min_i+","+min_j);
	if (indexlen > 2)
		minimize_Sij();	// repeat recursively until only 2 nodes are left
}
/**
 * ************************ Multiple Sequence Alignment ************************** <br>
 * <br>
 * Version 2.0     November 22, 1994<br>
 * Program (version 1.0, June 1989)<br>
 * by  JD Kececioglu, SF Altschul, DJ Lipman & R Miner<br>
 * <br>
 * Improvements (from 1.0 to 2.0) by SK Gupta, AA Schaffer,
 * with some guidance from JD Kececioglu<br>
 * <br>
 * Please cite papers 6 and 7 below if you wish to cite the software
 * package MSA.<br>
 * <br>
 * See:<br>
 *    1. Carrillo & Lipman, "The Multiple Sequence Alignment Problem in Biology",
 * 		SIAM J. Appl. Math. 48 (1988) 1073-1082;<br>
 *    2. Altschul & Lipman, "Trees, Stars, and Multiple Biological Sequence
 *		Alignment", SIAM J. Appl. Math. 49 (1989) 197-209;<br>
 *    3. Altschul, "Gap Costs for Multiple Sequence Alignment",
 *		J. Theor. Biol. 138 (1989) 297-309;<br>
 *    4. Altschul, Carroll & Lipman, "Weights for Data Related by a Tree",
 *		J. Molec. Biol. 207 (1989) 647-653;<br>
 *    5. Altschul, "Leaf Pairs and Tree Dissections",<br>
 *		SIAM J. Discrete Math. 2 (1989) 293-299;<br>
 *    6. Lipman, Altschul & Kececioglu, "A Tool for Multiple Sequence Alignment",
 *		Proc. Natl. Acad. Sci. USA 86 (1989) 4412-4415.<br>
 *    7. Gupta, Kececioglu & Schaffer, "Improving the Practical Time<br>
 *                and Space Efficiency of the Shortest-Paths Approach to<br>
 *                Sum-of-Pairs Multiple Sequence Alignment", J. Computational<br>
 *                Biology 2(1995) 459-472.<br>
 * <br>
 * Computes an optimal multiple alignment within a lattice defined by
 * the join of vertices from two-dimensional path graphs.  There
 * is one such path graph for each pair of sequences, and all vertices
 * contained in paths whose cost is within an epsilon of the minimal
 * are included.  Features include:<br>
 * <br>
 * The program computes a minimal SP alignment (Sum of Pairwise costs).
 * Each pairwise cost may be given a different weight.  These weights
 * may be calculated from an evolutionary tree using either of two
 * rationales (Altschul et al., "Weights for Data Related by a Tree",
 * J. Molec. Biol. 208) or equal weights may be specified.  The
 * evolutionary tree is estimated from pairwise distances using the
 * Neighbor Joining method (Saitou & Nei, Mol. Biol. Evol. 4:406-425).<br>
 * <br>
 * Epsilons for each sequence pair may be input by the user or
 * estimated.  An heuristic multiple alignment is computed and
 * the epsilons are taken to be the difference between the
 * imposed and minimal pairwise costs.<br>
 * <br>
 * Gap costs are of the affine type generalized for multiple sequence
 * alignments (Altschul, "Gap Costs for Multiple Sequence Alignment",
 * J. Theor. Biol. 138:297-309).  The user may specify whether terminal
 * gaps are counted.<br>
 * <br>
 * The user may specify residues in any of the sequences to be
 * forced into alignment.<br>
 * <br>
 * ******************************************************************************* <br>
 * compute multiple sequence alignment.
 * <main>
 * Erstellungsdatum: (18.01.2002 17:35:45)
 * @return align.Edge
 */
protected Edge msa() throws CancelException {

	int I, J;							// seq# loop counter
	int inc;							// measure for time progress '*' proceed
	int ccost= 0;						// also for time measuring (when is a '*' displayed)
	Vertex s, t;						// container for source and sink vertex
	Edge e, f;							// adjacent edges u(p) -e-> v(q) -f-> w(r)
	Vertex v, w;						// f-terminating edges; 	w unused !?!!
//  int[NUMBER+1] p, q, r are 'super-local' (static!) ??!!
	int difference;						// f.D - e.D
	char[] C= new char[NUMBER + 1];		// a column of chars (spanned over an edge)
	int[] delta0= new int[NUMBER + 1];	// offset[all Dimensions] for edge e
	int[] delta1= new int[NUMBER + 1];	// offset[] for edge f
	int[] ends= new int[NUMBER];		// seq# which are @ start/end position
	int endcount, endindex;				// max Pointer, loop Iterator for ends[]
	int d;								// (successively) sum up cost for each dim
	boolean mbreak= true; 							// dummy for break multiple loops ('goto')
	inc= 0;												// UNUSED if (oflag!=0), for compiler necessary ***
	endcount= 0;										// UNUSED if (oflag!=0), for compiler necessary ***
	int[] p= new int[NUMBER+1];					// 'super-local' from msa()
	int[] q= new int[NUMBER+1];					// 'super-local' from msa()
	int[] r= new int[NUMBER+1];					// 'super-local' from msa()



				// compute shortest paths to vertices in intersected region of lattice ??!
	s= source();						// create the source
	t= sink();							// create the sink (lastMan)
	h= heap(Upper);						// create the heap

	presource= create_vertex(null);		// dummy vertex before first vertex
	e= create_edge(presource, s);		// dummy edge leading to 1st vertex
	e.dist= 0;
	e.refer++; 								// make sure edge does not get freed
	e.backtrack= null;

	insert(e, h);						// put this 1st edge into the heap

	if (oflag != 0) {
		System.out.println("....1....2....3....4....5....6....7....8....9....0");
		inc= 1 + Upper / 50;			// measure for graphical time process inc
	}
	if (proxy!= null)
		proxy.setMaximum(49);


	while ((e= extract())!= null && (v= e.head) != t) {	// while(!heap.isEmpty && !lastManStanding)

		if ((proxy!= null)&& (e.dist > ccost))
			proxy.increase();
		if (oflag != 0 && e.dist > ccost) {
			System.out.print("*");		// fill up the 50 places of the scala
			ccost+= inc;				// inc by 1.1/50 of Upper
		}

		if (e.dist <= Upper) { 		// Homing costs heuristic criterium

			coord(e.tail, p); 			// put coordiantes of tail into p array,..
			safe_coord(e.head, q); 		// .. and coordinates of head of edge into q

			for (I= 1; I <= K; I++) 	// calculate offset[all Dimensions] of e (s. Fig.3, p 13)
				delta0[I]= q[I] - p[I];		// by substracting tail-point coordinates from head-point coordinates

			if (gflag != 0) {			// if end gaps are NOT penalized
				endcount= 0;
				for (I= 1; I <= K; I++)					// look in all dimensions
					if (q[I] == 0 || q[I] == N[I]) {	// if we are @ the start or the end of the sequence
						delta0[I]= 2;						// correct offset to 'endPosition' for T (cost-field)
						ends[endcount++]= I;				// save seqs which are @ start/end, increment counter
					}
			}

			for (I= 2; I <= K; I++)		// for (all combinations)
				for (J= 1; J < I; J++) {
					Tpointer[I][J][0]= T[delta0[I]][delta0[J]][0];	// converts sub-arrays of T with already
					Tpointer[I][J][1]= T[delta0[I]][delta0[J]][1];	// determined precursor positions of seq I, J
					Tpointer[I][J][2]= T[delta0[I]][delta0[J]][2];	// to seq-orientated TPointer[I][J][to][do]
				}

			if (v.out == null) 			// if first time visiting v (= e.head)
				adjacent(e, q);				// generate outgoing edges list

			for (f= v.out; f!= null; f= f.next) {	// iterate outgoing list
				mbreak= true;										// dummy for breaking loop
				difference= f.dist - e.dist;		// difference between adjacent edges
				if (mbreak&& difference > 0) {			// must be > 0 ???!
					safe_coord(f.head, r); 			// get coordinates of next vertex into r

					for (I= 1; I <= K; I++) {		// spanning over edge f:
						C[I]= (r[I] > q[I] ? S[I][r[I]] : getDASH()); 	// get Column of chars
						delta1[I]= r[I] - q[I];						// get offset vector
					}

					if (gflag != 0)		// if end gaps are NOT penalized (treated in another way)
						for (endindex= 0; endindex < endcount; endindex++)
							delta1[ends[endindex]]= 2;	// set all end indices of previous edge e here also
/*						for (I=1;I<=K;I++)
							if (q[I]==0||q[I]==N[I])				// seems to be more ineffective ??!
								delta1[I] = 2;						// (though more enlightening)
*/

					d= 0;						// init sum
					d += scale[1][2] * 			// sum for the first two seqs
						(sub(C[1],C[2]) + Tpointer[2][1][delta1[2]][delta1[1]]);	// D[C[1]][C[2]]
					for (I= K; I>=3; I--) {		// for rest of sequences (wo #1 and #2)
						if (d >= difference) {		// e.D + COST(f,e) < f.D; s. Fig.4
							mbreak= false; 								// goto nextedge; --> breaker mbreak (try with goto <label> ?!)
							break;										// prevents I-- in circular for loop
						}
						for (J= 1; J < I; J++)		// successively increase the sum (prob. time efficiency !?!)
							d += scale[I][J] *
								(sub(C[I],C[J]) + Tpointer[I][J][delta1[I]][delta1[J]]);	// D[C[I]][C[J]]

					}

					if (mbreak && d < difference) {	// check 'e.D + COST(f,e) < f.D' for last iteration (sum up) OR pre-matural break of preceeding for-loop
						delete(f, h);					// if succeeds, (temporarily) remove edge f from its old position in heap
						e.refer++;						// e has a new backtracking edge in some optimal path pointing to it
/*						if (f.backtrack != null)						// Garbage collector ??!
							if (--f.backtrack.refer== 0)
								free_edge(f.backtrack);
*/
						f.dist= d + e.dist;				// f.D= f.E + COST(f,e); s. Fig. 4
						f.backtrack= e;					// set the new backtracking edge to e
						insert(f, h);					// re-insert edge f in the heap

					} // if (d< difference)

				} // if (difference > 0)

				continue;												// <nextedge>; try with goto <label> ??!
			} // for(iterate adjacent edges list)
		} // if(e.dist<= Upper); HOMING COSTS

/*		if (e.refer == 0)												// Garbage collector ??!
			free_edge(e);
*/
	} // while (edges in the heap && !lastManStanding)

	if (oflag!= 0)
		System.err.println();

	return e;						// return null || lastManStanding

}
/**
 * Subroutine to calculate the optimal alignment of a new segment and
 * a multiple alignment<br>
 * <ecalc>
 * FIXED: A3= SS[]+ low inited originaly in align(), now ported here
 * Creation date: (20.01.2002 00:53:07)
 * @globalParam	 OAL[][] read already determined aligning pos relations
 * @globalReturn NAL[][] write into new aligning pos relations
 * @return int new length of partial layout now comprising seq I; either enlarged or same as before
 * @param I int	seq# of segment to be aligned now
 * @param low int start pos of alignment in segment I
 * @param n int length of alignment in segment I
 * @param m int length of (partial) layout of already aligned segments
 * @throws <code>Exception</code> if 'BIG not big enough'.
 * Also <code>ArrayIndexOutOfBoundsException</code> are upcasted if reducing BIG to 99,999.
 * @see #BIG
 */
protected int newal(int I, int low, int n, int m) throws Exception {

	int 	i,j,
			k,x,q,qq,
			dg,vg,hg,
			sum,	// sum of scale factors of all seqs Ord[x] < seq Ord[I] (x to I scale)
			sum2;	// sum only if start pos is no gap ?!
	int[]	sc= new int[NUMBER+1];	// scale factors of seq x to seq I (or vice versa)
	int[]	d= new int[2*LENGTH];
	int[] 	v= new int[2*LENGTH];
	int[]	h= new int[2*LENGTH];
	int[]	dn= new int[2*LENGTH];
	int[]	vn= new int[2*LENGTH];
	int[] 	hn= new int[2*LENGTH];

	for (sum=sum2=0,k=1;k<I;++k) {
		sum+= sc[k] = Ord[k]<Ord[I] ?	// init sc and update sum
			scale[Ord[k]][Ord[I]] : 		// scale[][] not mirrored on diagonal
			scale[Ord[I]][Ord[k]];
		if (OAL[k][1]!=0) 				// Con+1->NAL->OAL - gap (-1) is now (0)
			sum2+= sc[k];
//		System.out.println("sc"+sc[k]);
	}
//	System.out.println("sum= "+sum+", sum2= "+sum2);

	for (j=0;j<=m;j++) 					// init layout in size of already obtained layout length m
		hn[j]=dn[j]=vn[j]=BIG;

	for (i=0;i<=n;i++) {				// for (all positions in to aligning segment)
		for (j=0;j<=m;j++) {
			h[j]=hn[j]; d[j]=dn[j]; v[j]=vn[j];
		}
		if (i==0) {
			dn[0]=0;
			vn[0]=sum* (low!=0?G:GG);
			hn[0]=sum2*(low!=0?G:GG);
		}
		else {
			vn[0]=v[0]+sum*sub(A3[i+low],getDASH()); // D[A3[i+low]][DASH]
			hn[0]=dn[0]=BIG;
		}
		vv[i][0]=dd[i][0]=hh[i][0]=1;

			// Calculate optimal alignment in the case that terminal gap costs
			// are different than internal gap costs

		if (gflag!=0)
			for (j=1;j<=m;j++) {
				dg=d[j-1]; hg=h[j-1];
				for (x=0,k=1;k<I;++k) {
					q=OAL[k][j-1];
					qq=OAL[k][j];
					if (q<NN[k]) {
						dg += sc[k] * T[1][q>0?1:0][1][qq>0?1:0];
						if (q!=0 || i+low>1) hg +=
							sc[k]*T[0][q>0?1:0][1][qq>0?1:0];
					}
					x+= sc[k] * sub(SS[k][qq],A3[i+low]); // D [SS[k][qq]] [A3[i+low]]
				}
//				System.out.println("1.dn="+dg+","+hg+","+v[j-1]+"+"+x);
				dn[j]	= x + (k= min3(dg,hg,v[j-1]));
				dd[i][j]= (k==dg) ? DIAG : (k==hg ? HORZ : VERT);

				dg=d[j]; hg=h[j];
				for (k=1;k<I;++k)
					if ((qq=OAL[k][j])<NN[k]) {
						dg += sc[k] * T[1][qq>0?1:0][1][0];
						if (qq!=0 || i+low>1)
							hg+=sc[k] * T[0][qq>0?1:0][1][0];
					}
//				System.out.println("1.vn="+dg+","+hg+","+v[j]+"+"+(sum*sub(DASH,A3[i+low])));
				vn[j]=sum*sub(getDASH(),A3[i+low])+(k=min3(dg,hg,v[j])); // D[DASH][A3[i+low]]
                   // A. A. Schaffer fixed bug reported by A. Ropelewski
				if (k > BIG) {
				  throw new Exception("BIG is not big enough");
//				  System.err.print("\n BIG in defs.h is not big enough");
//				  return 0;		// System.exit(1);
				}

				vv[i][j]= ((k==dg) ? DIAG : ((k==hg) ? HORZ : VERT));

				dg=dn[j-1]; vg=vn[j-1]; hg=hn[j-1];
				for (x=0,k=1;k<I;++k) {
					q=OAL[k][j-1];
					qq=OAL[k][j];
					if (qq>1)
						hg+=sc[k] * T [0] [q>0?1:0] [0] [1];
					if (low+i<NN[I]) {
						dg += sc[k] * T [1] [q>0?1:0] [0] [qq>0?1:0];
						vg += sc[k] * T [1] [0]   [0] [qq>0?1:0];
					}
					x+= sc[k] * sub(SS[k][qq],getDASH());	// D[SS[k][qq]] [DASH]
				}
//				System.out.println("1.hn="+dg+","+hg+","+vg);
				hn[j]= x+(k=min3(dg,hg,vg));
				// A. A. Schaffer fixed bug reported by A. Ropelewski
				if (k > BIG){
					throw new Exception("BIG is not big enough");
//			 		System.err.print("\n BIG in defs.h is not big enough");
//			  		return 0;		// System.exit(1);
				}
				hh[i][j]= k==dg ? DIAG : (k==hg ? HORZ : VERT);
			}

			// Calculate optimal alignment in the case that terminal gap costs
			// are the same as internal gap costs

		else
			for (j=1;j<=m;j++) {	// for (all positions in layout)

									// determine DIAGONAL COSTS
				dg=d[j-1]; hg=h[j-1];
				for (x=0,k=1;k<I;++k) {
					q=OAL[k][j-1];
					qq=OAL[k][j];
					dg+= sc[k] * T [1] [q>0?1:0] [1] [qq>0?1:0];
					hg+= sc[k] * T [0] [q>0?1:0] [1] [qq>0?1:0];
					x+= sc[k] * sub(SS[k][qq],A3[i+low]);			// sum of mis-/matches; D [SS[k][qq]] [A3[i+low]]
				}
//				System.out.println("dn="+dg+","+hg+","+v[j-1]+"+"+x);
				dn[j]	= x+ (k= min3(dg,hg,v[j-1]));
				dd[i][j]= (k==dg) ? DIAG : (k==hg ? HORZ : VERT);	// preference: DIAG before HORZ before VERT (moving in I is important!)

									// determine VERTICAL COSTS
				dg=d[j]; hg=h[j];
				for (k=1;k<I;++k) {
					qq=OAL[k][j];
					dg+= sc[k] * T [1] [qq>0?1:0] [1] [0];	// moving in seq I, but not in seq k
					hg+= sc[k] * T [0] [qq>0?1:0] [1] [0];	// T[][][1][0]
				}
//				System.out.println("vn="+dg+","+hg+","+v[j]+"+"+(sum*sub(DASH,A3[i+low])));
				vn[j]	= sum*sub(getDASH(),A3[i+low]) + (k= min3(dg,hg,v[j]));	// D[DASH][A3[i+low]]
				vv[i][j]= (k==dg) ? DIAG : (k==hg ? HORZ : VERT);

									// determine HORIZONTAL COSTS
				dg=dn[j-1]; vg=vn[j-1]; hg=hn[j-1];
				for (x=0,k=1;k<I;++k) {
					q=OAL[k][j-1];
					qq=OAL[k][j];
					dg+= sc[k] * T [1] [q>0?1:0] [0] [qq>0?1:0];	// not moving in seq I T[][][0][]
					vg+= sc[k] * T [1] [0]   	 [0] [qq>0?1:0];
					hg+= sc[k] * T [0] [q>0?1:0] [0] [qq>0?1:0];
					x+= sc[k] * sub(SS[k][qq],getDASH());				// sum up gap costs for gap in seq I, why not sum2* D[..][-] ?!; D [SS[k][qq]] [DASH]
				}
//				System.out.println("hn="+dg+","+hg+","+vg);
				hn[j]	= x + (k=min3(dg,hg,vg));
				hh[i][j]= (k==dg) ? DIAG : (k==hg ? HORZ : VERT);
			} // end for (all pos in current layout of segment 1..m)
		} // end for (all pos of segment in seq I 1..n)


			// Traceback to reconstruct optimal alignment

		j=0;					// init counter for (evtl. extended) length of segment layout
		k= (dn[m]<=vn[m]) ? 	// init last step's direction
			(dn[m]<=hn[m]?DIAG:HORZ) :
			(vn[m]<=hn[m]?VERT:HORZ);
//		System.out.println("m"+m+"n"+n);
		while ( n!=0 || m!=0)	// while there are characters (of new segment I aligning) OR positions (of already aligned segments) left
			if (k==DIAG) {				// optimum was DIAGONAL
				for (i=1;i<I;++i) 			// for (all other already aligned seqs)
					NAL[i][j]= OAL[i][m];	// update newly aligned pos in other seq-arrays
				NAL[I][j++]=low+n;			// update newly aligned pos in new aligned seq (I)
				k=dd[n--][m--];			// traceback diagonal
			}
			else
				if (k==VERT) {			// optimum was VERTICAL
					for (i=1;i<I;++i)
						NAL[i][j]=0;
					NAL[I][j++]= low+ n;
					k=vv[n--][m];
				} else {				// optimum was HORIZONTAL
					for (i=1;i<I;++i)
						NAL[i][j]= OAL[i][m];
					NAL[I][j++]= 0;
					k=hh[n][m--];
				}

		return(j);	// return new length of (partial) segment layout
}
/**
 * Subroutine to choose an order in which to construct a progressive
 * multiple alignment.<br>
 * <ecalc>
 * Creation date: (20.01.2002 00:34:48)
 * @param oldp1 int
 * @param p1 int
 * @globalParam 	S[],N[],Con[][][]
 * @globalReturn	SS[]	re-ordered seqs
 * @globalReturn	NN[]	re-ordered |seqs|
 * @globalReturn	Ord[]	order of seqs (Ord[1]= seq x; S[x]= SS[1])
 */
protected void order(int oldp1, int p1) throws CancelException {

	int I,J,
		i,		// pointer in seq I, ranging from PI to end;;new average lowest scoring seq
		j,		// difference between aligned pos in J of last and actual pos in I
		M,		// loop counter (adapted j)
		i1,j1,	// actual lowest penalty pair of seqs I,J
		PI,PJ,	// position pointer (upcounting from start) for seq I,J
		end,	// stopping position in seq I
		t,		// last consistency value @ (i-1) from I with J;; dissimilarty of segments (sum of penaties)
		mind,	// minimum of dissimilarity between a seq-pair I,J;; min of 'average distance' to the most similar 2 seqs
		p;		// position pointer in aligned segments (of seq I,J)

	int[] a1= new int[LENGTH+1];				// aligned segments of seqs I,J
	int[] a2= new int[LENGTH+1];
	int[] test= new int[NUMBER+1];				// flag array showing already ordered (=0) and not yet ordered (=1) seqs
	int[][] dis= new int[NUMBER+1][NUMBER+1];	// dissimilarity scores segments I/J (sum of penalties over segments);
	int[] sum= new int[NUMBER+1];				// sum of dissimilarity scores of each seq to the ones already ordered

	i1= j1= i= 0;				// pre-init, not necessary
	a1[0]=a2[0]=getDASH();			// e-pos of aligned segments= '-'
	mind=BIG;					// init with max for finding minimum


				//	Calculate all pairwise costs for the segments in question

	for (I=1;I<K;++I)
		for (J=I+1;J<=K;++J) {	// for (all combinations of seqs)
			PJ=Con[1][J][oldp1];	// reference pos in seq J (aligned with seq1 @ pos oldp1)
			PI=Con[1][I][oldp1];	// starting @ pos in seq I which was aligned with seq[1][oldpos]
			end=Con[1][I][p1];		// stopping @ pos in seq I which was aligned with seq[1][p1]
			for (p=0,i=PI+1;i<=end;++i) {	// for all pos's in seq I (corresponding to [oldpos..p1] in seq1)
				j=(t=Con[I][J][i-1])>=0 ? 		// if (t=)last aligned pos between seqs I,J was no gap
					Con[I][J][i]-t:				// 		(j=) difference between pos's in seq J following up pos's in I were aligned
					(Con[I][J][i]>0?1:0);		// else: if (actual pos!= gap) j= 1, else (2 following gaps) j= 0
				for (M=j>0?j:1;
							// BUG FIX by micha:
							// for some reason M= j= Con[][][] is bigger than sequence length
							// if setting LENGTH to a bigger value
							// patched by (p< (a1.length- 1))
							(p< (a1.length- 1))&& (M!=0);--M) { 	// for all segment positions

					a1[++p]= M==1 ? S[I][++PI] : getDASH();	// write chars in segment alignment arrays
					a2[p]  = M<=j ? S[J][++PJ] : getDASH();	// (M<=j) only while loop proceeding (ie j=2..)
				}
			}


			for (t=0,i=1;i<p;++i) 			// for (all pos's in aligned segment)
				t+=	sub((char)a1[i],(char)a2[i])+	// add up match/mismatch/gap penalty; D[a1[i]][a2[i]]
					T[a1[i-1]!=getDASH()?1:0]		// + opening/continuing gap costs in both seqs
					 [a2[i-1]!=getDASH()?1:0] 		// (bool -> int)
					 [a1[i]!=getDASH()?1:0]
					 [a2[i]!=getDASH()?1:0];
			if (t<mind) { 					// check whether I,J is the lowest penalty score
				mind=t; 						// update new minimum penalty
				i1=I; 							// save I, J
				j1=J;
			}
			dis[I][J]=dis[J][I]=t;			// save dissimilarity (=segment penalty sum) between I and J
		}
	// end of double loop: for (all combinations of sequences do segment alignment, determine dis-score, find minimum penalty pair)


				//	Make lowest cost pair the first two segments in the order

	for (I=1;I<=K;++I) 	// init flag-array (1= not yet re-ordered)
		test[I]=1;

	Ord[1]=i1; SS[1]=S[i1]; NN[1]=N[i1]; test[i1]=0; // lowest penalty pair are
	Ord[2]=j1; SS[2]=S[j1]; NN[2]=N[j1]; test[j1]=0; // first two seqs in order

	for (I=1;I<=K;++I) 	// init sum[]
		sum[I]=dis[I][i1]+dis[I][j1];	// sum of dissimilarity of each seq to the lowest scoring ones

				//	Fill out the order using average distances

	for (j=3;j<=K;++j) {	// for all seqs but the two lowest scoring ones
		mind=BIG;
		for (I=1;I<=K;++I)	// find new min average penalty
			if (test[I]!=0 && sum[I]<mind) 	// exclude the already ordered seqs; find new min
				mind=sum[i=I];				// save new min in mind, seq which produced it in i
		Ord[j]=i; SS[j]=S[i]; NN[j]=N[i]; test[i]=0;	// re-order next lowest penalty seq
		for (I=1;I<=K;++I) 	// correct average scores to newly ordered seq
			sum[I]+=dis[I][i];	// add dissimilarity to just ordered seq
	}

// 	test out for order
/*	System.out.println("order");
	for (j=1; j<=K; j++)
		System.out.print(Ord[j]+", ");
	System.out.println();
*/

}
/**
 * output multiple	sequence alignment rows.
 * <main>
 * Erstellungsdatum: (21.01.2002 10:55:17)
 */
protected void output() {

	int	k, c;

	if (C==0)
		  return;
	for (k=1;k<=K;k++) {
		for (c=0;c<C;c++)
			System.out.print(M[k][c]);
		System.out.println();
	}
	System.out.println();
	C=0;
}


/**
 * Subroutine to calculate imposed cost for any pair of sequences.<br>
 * <ecalc>
 * Creation date: (20.01.2002 01:19:25)
 * @return int
 * @param I,J int sequences to calculate the optimal aligning costs from
 * @param b int ending position of layout
 * @globalParam AL[][]= complete alignment Layout
 */
protected int pcost(int I, int J, int b) throws CancelException {

	int i;
	int s=0;

	for (i=1;i<=b;++i)
		s+=sub((char)AL[I][i],(char)AL[J][i])+T[AL[I][i-1]!=getDASH()?1:0]	// D[AL[I][i]][AL[J][i]]
		[AL[J][i-1]!=getDASH()?1:0] [AL[I][i]!=getDASH()?1:0] [AL[J][i]!=getDASH()?1:0];

	i=1;
	while (AL[I][i]==getDASH() && AL[J][i]==getDASH())
		++i;
	if (AL[I][i]==getDASH() || AL[J][i]==getDASH())
		s+=(GG-G);

	i=b;
	while (AL[I][i]==getDASH() && AL[J][i]==getDASH())
		--i;
	if (AL[I][i]==getDASH() || AL[J][i]==getDASH())
		s+=(GG-G);

	return(s);
}
/**
 * function of primer.c
 * calculates pairwise costs (using current costfile)
 * and pairwise distances based on number of identities.
 * costs are used to estimate epsilons (ecalc()),
 * distances are used to calculate evolutionary tree (bias()).
 * Creation date: (19.01.2002 17:20:58)
 */
protected void primer() throws CancelException {

	char[] A, B;
	// act. 'super-local' for common use with faces()
	int I,J,i,j,Gi,Gj,n,m,q;
	int[] C= new int[LENGTH+1];		// forced alignment positions

	if (oflag!= 0) {
		for (i=(K*K-K)/2;i!= 0;--i)	// k*(k-1)/2 combinations
			System.err.print(".");
		System.err.println();
	}
	if (proxy!= null) {
		proxy.setMaximum((K*K-K)/2- 1);
	}


	if (fflag== 0)					// no user forced positions
		for (i=1;i<=LENGTH;++i)
			C[i]= 0;				// init forced positions with 0

	for (I=1;I<=K;I++)				// init Consistency for each seq with itself
		for (i=N[I];i>=0;i--)
			Con[I][I][i]= i;		// always aligned in same position (w itself)


	for (I=1;I<K;I++)				// all combinations of sequences
		for (n=N[I],A=S[I],J=I+1;J<=K;J++) {	// load comparing seq in A, length in n
			if (fflag!=0)
				fix(I,J,C,n);		// fix if user specified
			m=N[J];					// length of comparED seq
			B=S[J];					// comparED seq

							// compute distance form <0,0> to <i,j>

			dd[0][0]= 0;			// init first diagonal edge with 0
			hh[0][0]= vv[0][0]= GG;	// first vert & horiz edge penalized (G) or not (0)
			Con[I][J][0]= Con[J][I][0]= 0;	// pos 0 in each seq is e (empty string)
			for (j=1;j<=m;j++) {			// init upper end edges
				vv[0][j]= dd[0][j]= BIG;	// coming from nirvana
				hh[0][j]= hh[0][j-1]+ sub(getDASH(),B[j]);	// gap opening costs (from e) + multiples of gap penalizes
				Con[J][I][j]= -1;			// init comparing/compared row over complete length with (-1)= not aligned position in consistency
			}
			for (i=1;i<=n;i++) {			// init left end edges
				hh[i][0]= dd[i][0]= BIG;	// coming from nirvana
				vv[i][0]= vv[i-1][0]+ ((C[i]!=0) ? BIG : sub(A[i],getDASH()));	// inf. for fixed position OR gap opening costs (from e) + multiples of gap penalizes (since last fixed pos)
				Con[I][J][i]= -1;			// init comparing/compared row over complete length with (-1)= not aligned position in consistency
			}

			for (i=1;i<=n;i++)				// main double loop over 2D forward matrix (a-Matrix of both compared sequences)
				for (q=C[i], Gi= (i==n)?GG:G, j=1;j<=m;j++) {	// q marks fixed pos, Gi knows wheather this is end pos in seq 1 (I,i,n)
					Gj=(j==m)?GG:G;								// Gj determines terminal/normal gap penalize in 2nd seq (J,j,m)
					dd[i][j]= min3(dd[i-1][j-1],hh[i-1][j-1],vv[i-1][j-1])	// take min result leading to tail of diagonal edge
								+ (((q!=0)&&(q!=j))?BIG:sub(A[i],B[j]));		// add mis-/match cost or init with BIG if fixed pos (to support other incoming edge)
					hh[i][j]= min3(dd[i][j-1]+Gi,hh[i][j-1],vv[i][j-1]+Gi)	// take min result leading to tail of diagonal edge, add gap opening cost (terminal or normal) if previous is vert or diagonal (~ new opening gap in one this seq)
								+ sub(getDASH(),B[j]);							// add add gap penalize (for moving horizontal)
					vv[i][j]= min3(dd[i-1][j]+Gj,hh[i-1][j]+Gj,vv[i-1][j])	// take min result leading to tail of diagonal edge, add gap opening cost (term/normal,other seq pos determines) if path direction is changing now
								+ ((q!=0)?BIG:sub(A[i],getDASH()));				// add gap penalize OR init with BIG if fixed pos (to support other incoming edge)
				}

			costs[I][J]= min3(dd[n][m], hh[n][m], vv[n][m]);
			scale[J][I]= convert(I,J,n,m);
//			System.out.println("scale["+J+"]["+I+"]= "+scale[J][I]);

			if (oflag!= 0)
				System.err.print("*");		// reports 'combination computed'
			if (proxy!= null)
				proxy.increase();
/*	ifdef BIAS2
			if (scale[J][I]<=0) scale[J][I]=1;	// BIAS2 commented out
*/
		} // end of inner for loop (of the double loop iterating combinations)

	if (oflag!= 0)
		System.err.println();
}
/**
 * compute projected cost of edge <q,r> preceded by <p,q>.
 * <main>
 * Creation date: (20.01.2002 17:38:21)
 * @param p int[]
 * @param q int[]
 * @param r int[]
 */
protected void project(int[] p, int[] q, int[] r) throws CancelException {

	char[] C= new char[NUMBER+1];
	int I, J;
	int[][] t= new int[2][NUMBER+1];

	for (I=1;I<=K;I++) {
		C[I] = (r[I]>q[I] ? S[I][r[I]] : getDASH());
		t[0][I] = q[I] - p[I];
		t[1][I] = r[I] - q[I];
	}
	if(gflag!= 0)
          for (I=1;I<=K;I++)
            if (q[I]==0||q[I]==N[I])
              t[0][I]=t[1][I]=2;
	for (I=1;I<K;I++)
          for (J=I+1;J<=K;J++)
	    proj[I][J]+=sub(C[I],C[J])+T[t[0][I]][t[0][J]][t[1][I]][t[1][J]];	// D[C[I]][C[J]]

}
/**
 * recursively traverse the tree to compute the sum of path
 * lengths between A and B.
 * <bias.c>
 * Creation date: (19.01.2002 21:47:05)
 * @return int
 * @param A align.Node pointers to nodes
 * @param B align.Node
 */
protected int rdist(Node A, Node B) {

	if (A.sqn< 0) {			// if A is an internal node
		++count;
        return(rdist(A.lt, B) + rdist(A.rt, B));	// recurse on sons
	}
	else if (B.sqn < 0) {	// if B is an internal node
		++count;
		return(rdist(A, B.lt) + rdist(A, B.rt));	// recurse on sons
	}
							// if neither is an internal node
	return(Dij[A.sqn][B.sqn]);						// return distance between seqs

}
/**
 * traverses the tree in preorder form, writes node out to
 * a file and optionally the screen
 * Creation date: (19.01.2002 22:54:57)
 * @param A align.Node pointer to the root node
 */
protected void rpt(Node A) {

	if (oflag== 0)
		return;

	if (A.sqn>= 0) 			// leaf (with seq#)
		System.out.println("Leaf #"+(1+A.sqn)+"        Distance to parent = "+A.ltop);
	else {
		if (A.sqn== ROOT) 	// root (recognized by seq#= -9)
			System.out.println("---------------- Tree given from ancestor ----------------");
		else 				// internal node (seq#= -1)
			System.out.println("Internal Node  Distance to parent = "+A.ltop);
							// root and internal nodes have children
		System.out.print("On the left:   ");
		rpt(A.lt);			// iterate recursively
		System.out.print("On the right:  ");
		rpt(A.rt);			// for both sons
	}
}
/**
 * compute lattice coordinate of vertex.
 * <main>
 * Creation date: (20.01.2002 17:31:22)
 * @param v align.Vertex
 * @param p int[]
 */
protected void safe_coord(Vertex v, int[] p) {

	int	i;
	CoordinateValues a;		// *a (!)
	Coordinate b;			// *b (!)

	for (i=K, a=v.prev_coord_val; i>=1; i--, a=b.prev_coord_val) {
		b = a.curr_coord;
		p[i] = a.value;		// KICKED: p[i] = a - b.coord_vals + b.lo;
	}

}
/**
 * create sink vertex of lattice.
 * Creation date: (20.01.2002 16:55:43)
 * @return align.Vertex
 */
protected Vertex sink() {
	int	i;
	int[] p= new int[NUMBER+1];
	int[] index= new int[LENGTH+1];
	Coordinate a;
	CoordinateValues f;

	for (i=1;i<=K;i++)		// last point in universe
          p[i] = N[i];
          			// A2[1] still inited from source(); container ??!
	for (i=2,a=A2[1];i<=K;i++) {					// A= char[], A2= Coordinate[]
		f = a.coord_vals[p[i-1] - a.lo];			// array-pointer
		a = f.next_coord = create_coordinate(index,intersect(p,i,index),f);
	}
	f = a.coord_vals[p[K] - a.lo];					// array-pointer

				// (VERTEX *) (f.next_coord = (COORDINATE *)create_vertex(f))
	return (Vertex) (f.next_coord= (Coordinate) create_vertex(f));

}
/**
 * create source vertex of lattice.
 * Creation date: (20.01.2002 13:46:34)
 * @return align.Vertex
 */
protected Vertex source() {

	int[] p= new int[NUMBER+1];
	int i;
	int[] index= new int[LENGTH+1];
	Coordinate a;

	for (i=1;i<=K;i++)			// create Point [0,0,..,0]
          p[i] = 0;
	for (i=N[1];i>=0;i--)  		// index pos. over length of seq 1
	  index[i] = i;

	for (a=A2[1]=create_coordinate(index,N[1]+1,null), 		// A2[1] mis-used as container ??!
		i=2; i<=K; i++)			// for all other dimensions 2..K
		a=a.coord_vals[0].next_coord
		    = create_coordinate(index,intersect(p,i,index),a.coord_vals[0]);

    A2[1].refer++;							// Make sure coordinate does not get freed
    											// for re-use in sink() ??!

    return (Vertex) (a.coord_vals[0].next_coord= (Coordinate) create_vertex(a.coord_vals[0]));	// KICKED: a.coord_vals.next_coord= (COORDINATE *)create_vertex(a.coord_vals)

}
/**
 * associative array (with indices 'A' through 'Z')
 * Erstellungsdatum: (18.01.2002 15:10:02)
 * @return short
 * @param a char
 * @param b char
 */
protected short sub(char a, char b) throws CancelException {

//	return D[((byte) a)-'-'][((byte) b)-'-'];
	if (cancel)
		throw new CancelException();
	return D[(byte) a][(byte) b];
}
/**
 * Die Beschreibung der Methode hier eingeben.
 * Erstellungsdatum: (18.01.2002 15:10:02)
 * @return short
 * @param a char
 * @param b char
 */
protected void sub(char a, char b, short newValue) {

/*
	byte a_= (byte) (((byte) a)- (byte) '-');
	byte b_= (byte) (((byte) b)- (byte) '-');
*/
	D[a][b]= newValue;
	D[b][a]= newValue;
}
/**
 * recursively traverses the tree to compute the sum of path
 * lengths from a given node out to the leaves
 * <bias.c>
 * Creation date: (19.01.2002 21:54:06)
 * @return float
 * @param A align.Node a pointer to a node
 * @param total float cumulative path length
 */
protected float subdist(Node A, float total) {

	if (A.sqn>= 0) 			// if A isn't an internal Node
		return(total + A.ltop);			// sum up distance and return it

	count2++;				// else A IS internal Node
	return(subdist(A.lt,A.ltop+ total)	// send summed up distance to
		+ subdist(A.rt,A.ltop+ total)); // both sons recursively

}
/**
 * trace is a recursive function that traverses tree to find all leaf pairs
 * whose first member is a given leaf  (Rationale-1)
 * <bias.c>
 * Creation date: (19.01.2002 23:43:25)
 * @param prod float
 * @param no align.Node parent of a Node
 * @param sis align.Node sister (brother) of a Node
 */
protected void trace(float prod, Node no, Node sis, int i) {

	if (no.sqn> INOD) 			// if parent is leaf (corresponds to real seq#)
		B2[pN[i].sqn][no.sqn]= prod;	// pointer passed from whts()
	else
		if (sis==null) {		// if no brothers
			trace(prod/2, no.rt, null, i);
			trace(prod/2, no.lt, null, i);
		} else
			if (no.sqn!= ROOT) {	// if root
				trace(prod/2, sis, null, i);
				trace(prod/2, no.par, no.bro, i);
			} else 				// internal node with brothers, no root
				trace(prod, sis, null, i);

/* ifdef BIAS2		// (Rationale-2)

	if (no->sqn > INOD) B[(*pN)->sqn][no->sqn] = sum*prod;
	else if (sis==NULL) {
		trace(prod * no->lt->W, sum + no->rt->ltop, no->rt, NULL);
		trace(prod * no->rt->W, sum + no->lt->ltop, no->lt, NULL);
	}
	else {
		trace(prod * no->V, sum + sis->ltop, sis, NULL);
		if (no->sqn != ROOT)
			trace(prod * sis->W, sum + no->ltop, no->par, no->bro);
	}
*/
}
/**
 * create a vectro of integers.
 * <main>
 * Erstellungsdatum: (21.01.2002 19:01:37)
 * @return int[] of new length with values (as far/short as) provided
 * @param a int[] input array providing values for new array
 * @param n int new length of array
 */
protected int[] vector(int[] a, int n) {

	int[] v;
	int i;

	v= new int[n];
	for (i= 0; i< n; i++)
		v[i]= a[i];

	return v;
}
/**
 * ************************************************************************* <br>
 * Program written by SF Altschul<br>
 * Version 1.0   June 27, 1989<br>
 * <br>
 * Program to calculate pair weights given an evolutionary tree.<br>
 * <br>
 * See:	Altschul, "Leaf Pairs and Tree Dissections",<br>
 *		SIAM J. Discrete Math. 2 (1989).<br>
 *	Altschul, Carrol & Lipman, "Weights for Data Related by a Tree",<br>
 * 		J. Molec. Biol. 208 (1989). <br>
 * <br>
 * ************************************************************************* <br>
 * <bias.c>
 * Creation date: (19.01.2002 23:08:11)
 */
protected void whts() {

	int i,j;
	Node no;	// container for iterating node array
	float sm;	// smallest

				// Rationale-1 weights

	for (pN=Node.vlist,i=0; (no= pN[i]).sqn> INOD; ++i) // pointer translation: (pN=vlist; (no= pN).sqn> INOD; ++pN), vcount inited 0-based !
		trace((float)1.0, no.par, no.bro, i);			// iterate trace() until first internal node is found, pointer passing to subroutine

	for (sm=BIG,j=1;j<K;++j)
		for (i=0;i<j;++i) 			// all reverse combinations
			if (B2[i][j]<sm) 		// find smallest (value given by trace())
				sm=B2[i][j];
	for (i=0;i<K-1;++i)
		for (j=i+1;j<K;++j) { 		// all (reverse) combinations
			scale[i+1][j+1]=(int) (B2[i][j]/sm+0.5);	// middle values to sm, is (int) cast correct ??!
//			System.out.println("2scale["+(i+1)+"]["+(j+1)+"]= "+scale[i+1][j+1]);
		}

/*
#ifdef BIAS2			// Rationale-2 weights

						// Calculate the weights of all trees hanging from all internal nodes
	for (pN=vlist; (no= *pN)->sqn > INOD; ++pN) {
		no->w = 1.0;
		no->W = no->ltop;
	}
	for (; (no= *pN)->sqn > ROOT; ++pN) {
		no->w = no->lt->w * no->rt->W + no->lt->W * no->rt->w;
		no->W = no->ltop  * no->w     + no->lt->W * no->rt->W;
	}
	no->V = 1;
	no->v = 0;
	do {
		no= *(--pN);
		no->v = no->par->v * no->bro->W + no->par->V * no->bro->w;
		no->V = no->ltop   * no->v      + no->par->V * no->bro->W;
	}
	while (pN != vlist);

						// Calculate weights for leaf pairs using precomputed subtree weights
	for (; (no= *pN)->sqn >INOD; ++pN)
		trace(1.0,no->ltop,no->par,no->bro);

						// Scale pair weights so that smallest is about 8
	sm=1.0E+30;
	for (j=1;j<K;++j) for (i=0;i<j;++i) if (B[i][j]<sm) sm=B[i][j];
	sm /= 7.9;
	for (i=0;i<K-1;++i) for (j=i+1;j<K;++j) scale[i+1][j+1]=B[i][j]/sm+0.5;

#else					// Rationale-1 weights
*/
}
}
