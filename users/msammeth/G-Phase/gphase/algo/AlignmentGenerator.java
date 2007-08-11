/*
 * Created on Apr 8, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.algo;

import java.awt.Point;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.ensembl.driver.compara.HomologyAdaptor;

import gphase.Constants;

import gphase.Toolbox;
import gphase.db.EnsemblDBAdaptor;
import gphase.ext.ClustalWrapper;
import gphase.ext.DialignWrapper;
import gphase.ext.MLaganWrapper;
import gphase.ext.MuscleWrapper;
import gphase.ext.TCoffeeWrapper;
import gphase.io.ALNWrapper;
import gphase.io.GAFWrapper;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.GeneHomology;
import gphase.model.Graph;
import gphase.model.Species;
import qalign.OSChecker;
import qalign.algo.CostTable;
import qalign.tools.FASTAWrapper;
import qalign.tools.MSFWrapper;

/**
 * 
 * 
 * @author micha
 */
public class AlignmentGenerator {

	public static final String ALIGNMENT_SUBDIR= "ali"; 
	Graph graph;	
	
	static String getDate() {
		return "["+ new Date(System.currentTimeMillis()).toString()+ "]";
	}
	
	/**
	 * 
	 */
	public AlignmentGenerator(Graph newGraph) {
		graph= newGraph;
	}

	public static void main(String[] args) {
		
			// load graph 
		EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
		Graph g= adaptor.getGraphAllHomologs();
//		Logger.println("Graph filtered "+g.countGenesTranscriptsExons()+"---", false);
		
			// generate alignments
		AlignmentGenerator ag= new AlignmentGenerator(g);
		ag.alignHomologs(g);
		// ag.analyzeAlignments();
		// ag.writeOutTimes();
	}
	
	public void alignHomologs(Graph gr) {
		
		Gene[] baseGenes= gr.getSpecies()[0].getGenes();
		Vector blacklist= new Vector();
		int counter= 0;
		for (int j = 0; j < baseGenes.length; j++) {		// along genes

			int x;
			for (x = 0; x < blacklist.size(); x++)			// check whether already iterated 
				if (blacklist.elementAt(x)== baseGenes[x]) {
					blacklist.remove(x);	// can happen only once
					break;
				}
			if (x< blacklist.size())
				continue;

			GeneSet set= new GeneSet(gr);	// add
			set.add(baseGenes[j]);
			//align(set.getSet());
			++counter;
			
			Gene[] adds= set.getGeneParalogs(baseGenes[j]);	// update blacklist
			for (int i = 0; i < adds.length; i++) 
				blacklist.add(adds[i]);
		}
		System.out.println(counter+ " sets to align.");
	}
	
	public void writeOutTimes() {
		
		PrintStream oldOut= System.out;
		try {
			 System.setOut(new PrintStream(new FileOutputStream(Constants.DATA_DIR+ File.separator+ "times.tst")));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String[] cList= new File(Constants.DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ File.separator+ "clustal").list();
		Arrays.sort(cList);
		String[] dList= new File(Constants.DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ File.separator+ "dialign").list();
		Arrays.sort(dList);
		String[] mList= new File(Constants.DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ File.separator+ "mlagan").list();
		Arrays.sort(mList);

		
		for (int i = 0; i < cList.length; i++) {
			if (Arrays.binarySearch(dList, cList[i])< 0)
				continue;
			if (Arrays.binarySearch(mList, cList[i])< 0)
				continue;

			System.out.print(cList[i]+ " ");
			MSFWrapper msf= new MSFWrapper(
					Constants.DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ 
					File.separator+ "clustal"+ File.separator+ cList[i]
			);
			try {msf.read();} 
			catch (Exception e) {
				String messy= e.getMessage();
				if (messy.indexOf("min")< 0) {
					System.out.println();
					continue;
				}
				System.out.print(messy.substring(0, messy.indexOf("min"))+ " ");
//				Pattern patty= Pattern.compile("(\\d+)\\D*");
//				Matcher matty= patty.matcher(e.getMessage());
//				System.out.print(matty.group(1)+ " ");
				System.out.flush();
			}

			msf= new MSFWrapper(
					Constants.DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ 
					File.separator+ "dialign"+ File.separator+ cList[i]
			);
			try {msf.read();} 
			catch (Exception e) {
				String messy= e.getMessage();
				if (messy.indexOf("min")< 0) {
					System.out.println();
					continue;
				}
				System.out.print(messy.substring(0, messy.indexOf("min"))+ " ");
				System.out.flush();
			}
			
			msf= new MSFWrapper(
					Constants.DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ 
					File.separator+ "mlagan"+ File.separator+ cList[i]
			);
			try {msf.read();} 
			catch (Exception e) {
				String messy= e.getMessage();
				if (messy.indexOf("min")< 0) {
					System.out.println();
					continue;
				}
				System.out.print(messy.substring(0, messy.indexOf("min"))+ " ");
				System.out.flush();
			}

			System.out.println();
		}
		System.out.flush();
		System.out.close();			
	}

	public void analyzeAlignments() {
			
			PrintStream oldOut= System.out;
			String[] baseList= new File(Constants.DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR).list();
			for (int i = 0; i < baseList.length; i++) {
				if (baseList[i].contains("tba")
					|| baseList[i].contains("clustal")
					|| baseList[i].contains("dialign"))
					continue;	// skip
				try {
					 System.setOut(new PrintStream(new FileOutputStream(Constants.DATA_DIR+ File.separator+ "benchOut."+ baseList[i])));
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
				String[] list= new File(Constants.DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ 
						File.separator+ baseList[i]).list();
				if (list== null)
					continue;		// sort for reproducible recover after crash
				Arrays.sort(list);
				for (int j = 0; j < list.length; j++) {
					MSFWrapper msf= new MSFWrapper(
							Constants.DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ 
							File.separator+ baseList[i]+ File.separator+ list[j]
					);
					try {msf.read();} 
					catch (Exception e) {
//						e.printStackTrace();
//						System.err.println(
//								msf.getAbsFileName().substring(msf.getAbsFileName().lastIndexOf(File.separator))+ "\n"
//								+ e.getMessage());					
//						System.err.flush();
						if (msf.getSeqNames()== null|| msf.getSequences()== null)
							continue;
					}
					
	//				System.out.println(
	//						"["+ new Date(System.currentTimeMillis())+ "]\t"+
	//						msf.getAbsFileName().substring(msf.getAbsFileName().lastIndexOf(File.separator)));
					analyzeAli(msf.getSeqNames(), msf.getSequences());
	
				}
				System.out.flush();
				System.out.close();			
			}
		}
	
	void analyzeAli(String[] names, String[] layout) {
		
		Gene[] genes= new Gene[names.length];
		for (int i = 0; i < genes.length; i++) 
			genes[i]= graph.getGene(names[i]);
		
			// count tuples
		for (int p = 0; p < genes.length; p++) {
			for (int q = (p+1); q < genes.length; q++) {
				 
				// TODO filter for overlapping exons
				int i= genes[p].getStart();
				int j= genes[q].getStart();
				int k= 0; 
				
				int xalign= 0;
				int ralign= 0;
				while (k< layout[p].length()) {
					
					if (Constants.IS_GAP(layout[p].charAt(k))) {
						if (!Constants.IS_GAP(layout[q].charAt(k))) 
							j++;						
					} else {
						if (!Constants.IS_GAP(layout[q].charAt(k))) {	// 2positions
							Exon e1= genes[p].getExon(i);
							Exon e2= genes[q].getExon(j);
							if (e1== null^ e2== null)
								xalign++;
							else if (e1!= null&& e2!= null) {
								if ((e1.getHomolog(e2.getGene())== e2)||
										(e1.getHomolog(e2.getGene())== null&&
										e2.getHomolog(e1.getGene())== null))
									ralign++;	// allow free alignment if no ubrh was found
							}
							j++;
						}
						i++;
					}
					k++;
				}

				float SI= ((float) ralign)/ ((float) (ralign+ xalign));				
				System.out.print(SI+" ");
//				System.out.println(genes[p]+ " x "+ genes[q]+": "+SI+"("+ralign+", "+xalign+")");
			}
		}
		System.out.println();
		
	}

	
	void generateAlignments() {

		Species[] sp= graph.getSpecies();
		Vector sVec= new Vector();
		for (int i = 0; i < sp.length; i++) 
			sVec.add(sp[i]);
			
		Vector s;
		Vector t;
		for (int i = 0; i < sVec.size(); i++) {
			t= new Vector();
			s= ((Vector) sVec.clone());
			t.add(s.remove(i));
			generateAlignments(s,1,2,3,t);
		}
	}
	
	void generateAlignments(Vector sVec, int depth, int start, int end, Vector t) {

			// recursion
		Vector s;
		for (int i = 0; i < sVec.size(); i++) {
			s= ((Vector) sVec.clone());
			t.add(s.remove(i));
			generateAlignments(s,++depth,start,end,t);
		}
		
			// alignment
		if (depth>= start&& depth<= end) 
			align(t);
	}

	void align(Vector speciesList) {
		
			// find min
		int min= -1;
		int minVal= -1;
		for (int i = 0; i < speciesList.size(); i++)  
			if (((Species) speciesList.elementAt(i)).getGeneNb()> minVal) {
				min= i;
				minVal= ((Species) speciesList.elementAt(i)).getGeneNb();
			}
		
			// iterate genes and intersect
		Iterator it= ((Species) speciesList.elementAt(min)).getGeneIterator();
		while (it.hasNext()) {
			Gene g= (Gene) it.next();
			Vector rest= (Vector) speciesList.clone();
			rest.remove(min);
			Species[] rr= new Species[rest.size()];
			for (int i = 0; i < rr.length; i++) 
				rr[i]= (Species) rest.elementAt(i);

			if (g.hasHomologyWith(rr)) {
				Species[] rrr= new Species[rr.length+ 1];
				for (int i = 0; i < rr.length; i++) 
					rrr[i]= rr[i];
				rrr[rrr.length- 1]= ((Species) speciesList.elementAt(min));
//				align(g, rrr);
			}
		}
	}
	

	
	public String[] alignTCoffee(String[] names, String[] concatSeqs) {
		
			// align sequences
		String fName= writeOutTemp(names, concatSeqs);
		TCoffeeWrapper wrap= new TCoffeeWrapper(fName);
		wrap.setOutputMSF(true);
		wrap.runTCoffee(fName, "");
		
			// get result
		MSFWrapper msfReader= new MSFWrapper(fName.substring(0, fName.lastIndexOf('.'))+ ".msf_aln");
		try {
			msfReader.read();		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		String[] layout= msfReader.getSequences();
		String[] mapNames= msfReader.getSeqNames();
		String[] mapSeqs= new String[layout.length];
		for (int i= 0; i < mapSeqs.length; i++) {
			mapSeqs[i]= layout[i];
		}
		for (int i= 0; i < mapNames.length; i++) {
			if (!mapNames[i].equalsIgnoreCase(names[i]))
				for (int j= 0; j < names.length; j++) {
					if (mapNames[i].equalsIgnoreCase(names[j]))
						mapSeqs[i]= layout[j];
				}
		}
		
		return mapSeqs;
	}

	public static AlignmentWrapper alignClustal(String fName) {
	
		ClustalWrapper clustal= new ClustalWrapper(fName);
		clustal.setCostType(CostTable.DNARNA);
		clustal.startGrid();
		
		return clustal;
	}
	
	public void align(Gene[][] set) {

		String[] names= new String[set.length];
		String[] seqs= new String[set.length];
		
		align(names, seqs, set, 0);	// start recursion
	}
	
	public void align(String[] names, String[] seqs, Gene[][] set, int depth) {
		
		if (depth> set.length) {	// recursion abort -> align
			align(names, seqs);
			return;
		}
		
		for (int i = 0; i < set[depth].length; i++) {		// all possibilities for this depth
			names[depth]= set[depth][i].getStableID();
			seqs[depth]= graph.readSequence(set[depth][i]);
			align(names, seqs, set, depth+1);				// recursion
		}
		
	}
	
	public void align(String[] names, String[] seqs) {

		String fName= writeOutTemp(names, seqs);
		String bName= "";
		for (int i = 0; i < names.length; i++) 
			bName+= names[i]+ "_";
		bName= bName.substring(0, bName.length()- 1);
		System.out.println(bName);
		String tName;
		
			// mlagan
		alignMLagan(names, seqs, bName);
//			// clustal
//		tName= Constants.DATA_DIR+ File.separator+
//			ALIGNMENT_SUBDIR+ File.separator+
//			"clustal"+ File.separator+
//			bName+ ".msf";
//		if (!new File(tName).exists()) {
//			try{new File(tName).createNewFile();}	// block
//			catch (IOException e) {;}
//			AlignmentWrapper wrapper= alignClustal(fName);
//			System.out.println(getDate()+" mlagan started: "+ fName);
//		} else 
//			System.out.println(getDate()+" mlagan skipped: "+ fName);
//		
//		
//			// dialign
//		tName= Constants.DATA_DIR+ File.separator+
//		ALIGNMENT_SUBDIR+ File.separator+
//		"dialign"+ File.separator+
//		bName+ ".msf";
//		if (!new File(tName).exists()) {
//			AlignmentWrapper wrapper= alignDialign(fName);
//			System.out.println(getDate()+" dialign started: "+ fName);
//		} else
//			System.out.println(getDate()+" dialign skipped: "+ fName);
//	
		
			// muscle
//		tName= Constants.DATA_DIR+ File.separator+
//		ALIGNMENT_SUBDIR+ File.separator+
//		"muscle"+ File.separator+
//		bName+ ".msf";
//		if (!new File(tName).exists()) {
//			long t0= System.currentTimeMillis();
//			AlignmentWrapper wrapper= alignMuscle(fName);
//			long t1= System.currentTimeMillis();
//			MSFWrapper msf= new MSFWrapper(tName);
//			msf.setDescription(((t1-t0)/ 60000)+ " min");
//			msf.setSeqNames(names);
//			msf.setSequences(wrapper.getLayout());
//			try {msf.write();} catch(Exception e){ 
//				try {new FileWriter(tName).write(((t1-t0)/ 60000)+ " min"); }	// fail
//				catch (IOException ex) {;} 
//			}
//			Logger.println("\tmuscle aligned", false);
//		} else
//			Logger.println("\tfile exists, skipped muscle alignment.", false);
	
			// tcoffee
//		tName= Constants.DATA_DIR+ File.separator+
//		ALIGNMENT_SUBDIR+ File.separator+
//		"tcoffee"+ File.separator+
//		bName+ ".msf";
//		if (!new File(tName).exists()) {
//			long t0= System.currentTimeMillis();
//			AlignmentWrapper wrapper= alignTCoffee(fName);
//			long t1= System.currentTimeMillis();
//			MSFWrapper msf= new MSFWrapper(tName);
//			msf.setDescription(((t1-t0)/ 60000)+ " min");
//			msf.setSeqNames(names);
//			msf.setSequences(wrapper.getLayout());
//			try {msf.write();} catch(Exception e){ 
//				try {new FileWriter(tName).write(((t1-t0)/ 60000)+ " min"); }	// fail
//				catch (IOException ex) {;} 
//			}
//			Logger.println("\ttcoffee aligned", false);
//		} else
//			Logger.println("\tfile exists, skipped tcoffee alignment.", false);
		
		new File(fName).delete();
	}
	
	public static AlignmentWrapper alignTCoffee(String fName) {
		
			TCoffeeWrapper wrapper= new TCoffeeWrapper(fName);
			wrapper.start();
			
			return wrapper; //wrapper.getScore();
	}
	 
	public static AlignmentWrapper alignDialign(String fName) {
		
			DialignWrapper wrapper= new DialignWrapper();
			wrapper.runDialign(fName);

			return wrapper;
	}	
	
	public static AlignmentWrapper alignMuscle(String fName) {
		
			MuscleWrapper wrapper= new MuscleWrapper(fName, fName+ ".mus");
			wrapper.execute();

			return wrapper;
	}	
	
	public void alignMLagan(String[] names, String[] seqs, String bName) {
		
		String tName= Constants.DATA_DIR+ File.separator+
		ALIGNMENT_SUBDIR+ File.separator+
		"mlagan"+ File.separator+
		bName+ ".msf";
		String tName2= Constants.DATA_DIR+ File.separator+
		ALIGNMENT_SUBDIR+ File.separator+
		"mlagan"+ File.separator+
		bName+ ".mfas";
		if (new File(tName).exists()|| new File(tName2).exists()) 
			System.out.println(getDate()+" mlagan skipped: "+ bName);
		else {
//			try{new File(tName).createNewFile();}	// block
//			catch (IOException e) {;}

			String[] nnames= new String[names.length];		// add trivial names
			for (int i = 0; i < names.length; i++) 
				nnames[i]= names[i]+ "| "+ Species.getCommonNameForPrefix(Gene.getSpeciesPfx(names[i]));
			String[] sfNames= new String[names.length];
			for (int i = 0; i < sfNames.length; i++)		// write out single files 
				sfNames[i]= writeOutTemp(new String[]{names[i]}, new String[]{seqs[i]});
			
			final String MLAGAN_SH= "mlagan.sh"; 
			String sh= "export LAGAN_DIR="+ Constants.HOME_DIR+ File.separator+ Constants.MLAGAN_SUBDIR+ "_$SGE_ARCH"+ 
				Constants.MLAGAN_SUBSUBDIR+ "\n";			
			sh+= "time "+Constants.HOME_DIR+ File.separator+ Constants.MLAGAN_SUBDIR+ "_$SGE_ARCH"+ 
				File.separator+ Constants.MLAGAN_EXEC+ " "; 
			for (int i = 0; i < sfNames.length; i++) 
				sh+= sfNames[i]+" ";
			sh+= "-tree \""+Constants.getTree(3)+"\""+ 
				" -out "+Constants.getMlaganOutDir()+ File.separator+ bName+ ".mfas 2>&1\n";	// only can mfasta
			sh+= "echo $HOSTNAME >&2\n"+
				"mv -vi "+Constants.getMlaganOutDir()+ File.separator+ MLAGAN_SH+ ".o$JOB_ID "
					+ Constants.getMlaganOutDir()+ File.separator+ bName+ ".log\n";
			sh+= "mv -vi "+Constants.getMlaganOutDir()+ File.separator+ MLAGAN_SH+ ".e$JOB_ID "
					+ Constants.getMlaganOutDir()+ File.separator+ bName+ ".inf";
			
			try {
				BufferedWriter buffy= new BufferedWriter(new FileWriter(
						Constants.getScratchDir()+ File.separator+ MLAGAN_SH));
				buffy.write(sh);
				buffy.flush();
				buffy.close();
			} catch (IOException e) {
				System.err.println("Problems writing shell script "+ MLAGAN_SH);
			}

			try {
				Runtime.getRuntime().exec(Constants.ALIGNMENT_SUBDIR+ File.separator+ "qsub_mfas.sh");
			} catch (IOException e) {
				System.err.println("Problems running shell script "+ MLAGAN_SH);
			}

			System.out.println(getDate()+" mlagan started: "+ bName);
		}
	}	

	public static String writeOutTemp(String[] names, String[] seqs) {
		
			// init file
		File tFile= null;
		try {
			tFile= File.createTempFile("tmp", ".tfas", new File(Constants.getScratchDir()));
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
		String tName= tFile.getAbsolutePath();
		tFile.deleteOnExit();

			// write
		FASTAWrapper fasta= new FASTAWrapper(tName);
		try {	
			fasta.writeFASTA(names, seqs);
		} catch (Exception e) {
			return null;
		}
		
		return tName;
	}

	void writeGraph() {
		
			// write graph
		try {
			ObjectOutputStream oo= new ObjectOutputStream(new FileOutputStream(
					Constants.DATA_DIR+ File.separator+ GRAPH_FNAME)); 
			oo.writeObject(graph);
			oo.flush();
			oo.close();
			System.out.println("["+ new Date(System.currentTimeMillis())+ "]: Graph wrote ---");
		} catch (Exception e) {
			System.err.println("Error writing graph: "+ e.getMessage());
			e.printStackTrace();
		}		
	}

	static final String GRAPH_FNAME = "graph.oos";
}
