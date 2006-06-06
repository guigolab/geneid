package qalign.algo.msa;

/**
 * <bias.c>
 * nodes for the tree.
 * (alias tn)
 * Creation date: (19.01.2002 20:50:34)
 * @author:
 */
public class Node {
	public int sqn;		// sequence number or -1 for internal nodes
	public float ltop;		// distance to parent node
	public float w;
	public float W;
	public float v;
	public float V;
	public Node lt,rt;		// pointers to child nodes
	public Node par;		// pointer to parent node
	public Node bro;		// pointer to brother node

	public static int vcount;							// <bias> inited in MSA.init_data()
	public static Node[] vlist= new Node[2*MSA.NUMBER];// <bias> exported from msa-globals

/**
 * Die Beschreibung der Methode hier eingeben.
 * Erstellungsdatum: (22.01.2002 15:44:12)
 */
public Node() {}
/**
 * creates nodes for the tree.
 * Creation date: (19.01.2002 21:28:22)
 * @param sqn  - sequence number or -1 for internal nodes
 * @param dtop - distance to parent node
 * @param lt,rt - pointers to child nodes
 */
public Node(int sqn, float ltop, Node lt, Node rt) {

    vlist[vcount++]= this;		// save in list

	this.ltop = ltop;				// init attributes
	this.lt = lt;
	this.rt = rt;
	this.sqn = sqn;

	if (sqn<0) {				// if internal node
		this.lt.bro = rt;				// chain children
		this.rt.bro = lt;				// as brothers
		this.lt.par= rt.par= this;		// and let them know their father
	}

}
}
