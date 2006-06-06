/*
 * Created on Mar 13, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.io;

import gphase.db.EnsemblDBAdaptor;
import gphase.gui.Paparazzi;
import gphase.model.AbstractRegion;
import gphase.model.DefaultRegion;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.util.Vector;

import com.p6spy.engine.logging.appender.StdoutLogger;

/**
 * 
 * 
 * @author msammeth
 */
public class Tyler {

	final static String TYLER_SUBDIR= "tyler";
	final static String TYLER_OUTPUT_SUBDIR= "pics";
	final static String TYLER_DEFAULT_INPUT_FNAME= "u12intronList.txt";
	String fName= null;
	String[] geneNames= null;
	DefaultRegion[] u12Introns= null;
	
	public static void main(String[] args) {
		
		EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
		Graph g= adaptor.getGraphAllGenes(new Species(EnsemblDBAdaptor.SPECIES_SMALL[0]));
		
		Tyler myTyler= new Tyler(TYLER_SUBDIR+ File.separator+ TYLER_DEFAULT_INPUT_FNAME);
		String[] names= myTyler.getGeneNames();
		Gene[] genes= new Gene[names.length];
		for (int i = 0; i < names.length; i++) {	// convert names to genes (lookup)
			genes[i]= g.getGene(names[i]);
			if (genes[i]== null)
				System.err.println("Gene not found "+ names[i]);
		}
		Paparazzi.burst(genes, myTyler.getU12Introns(), new File(TYLER_SUBDIR+File.separator+TYLER_OUTPUT_SUBDIR));
	}	
	
	public Tyler(String newFName) {
		this.fName= newFName;
		readIn();
	}
	
	boolean readIn() {

		Vector geneNamesV= new Vector();
		Vector iStartV= new Vector();
		Vector iEndV= new Vector();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fName));
			while (buffy.ready()) {
				String s= buffy.readLine().trim();
				StringTokenizer toki= new StringTokenizer(s, "\t");
				toki.nextElement();	// internal number
				geneNamesV.add(toki.nextElement());
				toki.nextElement();	// chromosome name
				iStartV.add(toki.nextElement());
				iEndV.add(toki.nextElement());
			}						
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		
		geneNames= new String[geneNamesV.size()];
		u12Introns= new DefaultRegion[geneNamesV.size()];
		for (int i = 0; i < geneNames.length; i++) {
			geneNames[i]= (String) geneNamesV.elementAt(i);
			u12Introns[i]= new DefaultRegion(
					Integer.parseInt((String) iStartV.elementAt(i)), 
					Integer.parseInt((String) iEndV.elementAt(i))
			);
		}
		
		return true;
	}
	
	private void kissGina() {
		System.out.println("mua!");
	}
	public String[] getGeneNames() {
		return geneNames;
	}
	public DefaultRegion[] getU12Introns() {
		return u12Introns;
	}
}
