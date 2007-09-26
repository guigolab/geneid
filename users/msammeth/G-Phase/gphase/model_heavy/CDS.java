/*
 * Created on Mar 3, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model_heavy;

import java.io.Serializable;

/**
 * 
 * 
 * @author micha
 * @deprecated not yet used
 */
public class CDS implements Serializable {

	Exon exon= null;	// array?
	Transcript[] transcripts= null;
	
	int start= -1;
	int stop= -1;
	int frame= -1;
}
