/*
 * Created on Mar 8, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import gphase.Constants;
import gphase.db.EnsemblDBAdaptor;

/**
 * 
 * 
 * @author msammeth
 */
public class GraphHandler {

	public static final String getGraphAbsPath() {
		return Constants.HOME_DIR+ File.separator+ GRAPH_SUBDIR+ File.separator+ GRAPH_FNAME; 
	}

	public static final String getGraphAbsPath(Species spec) {
		return getGraphAbsPath(new Species[] {spec});
	}
	
	public static final String getGraphAbsPath(Species[] spec) {
		
		String result= "";
		
		for (int i = 0; i < spec.length; i++) 
			result+= spec[i].getCommonName()+"_";
		
		return Constants.HOME_DIR+ File.separator+ GRAPH_SUBDIR+ File.separator+ result+ GRAPH_FNAME; 
	}

	public static Graph readIn(String fName) {
	
		Graph g= null;
		try {
			ObjectInputStream oi= new ObjectInputStream(new FileInputStream(fName));
			g= (Graph) oi.readObject();
			oi.close();
			System.out.println("Graph read "+g.countGenesTranscriptsExons()+"---");
		} catch (FileNotFoundException e) {
			System.err.println("Graph does not exist!");
			//System.err.println(e);
		} catch (Exception e) {
			System.err.println("Error reading graph: "+ e.getMessage());
			e.printStackTrace();
		}		
		
		return g;
	}

	public static void writeOut(Graph g) {
		writeOut(g, getGraphAbsPath(g.getSpecies()));
	}
	
	public static void writeOut(Graph g, String fName) {
		
			// write graph
		try {
			ObjectOutputStream oo= new ObjectOutputStream(new FileOutputStream(fName)); 
			oo.writeObject(g);
			oo.flush();
			oo.close();
		} catch (Exception e) {
			System.err.println("Error writing graph: "+ e.getMessage());
			e.printStackTrace();
		}		
	}

	public static final String GRAPH_ENCODE_FNAME = "graph_enc.oos";
	public static final String GRAPH_FNAME = "graph.oos";
	public static final String GRAPH_SUBDIR = "graph";

}
