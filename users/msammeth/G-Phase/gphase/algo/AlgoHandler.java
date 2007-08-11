package gphase.algo;

import gphase.Constants;
import gphase.io.gtf.GTFChrReader;
import gphase.model.Gene;
import gphase.model.Species;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.HashMap;

import sun.security.action.GetLongAction;

public class AlgoHandler {

	public static String getGencodeAnnotation() {
		File annotDir= new File("annotation");
		String[] dir= annotDir.list();
		for (int i = 0; i < dir.length; i++) {
			if (dir[i].toUpperCase().contains("GENCODE"))
				return new File(annotDir+ File.separator+ dir[i]).getAbsolutePath();
		}
		return null;
	}
	public static Method[] getMethods(Class target, String[] mNames) {
		Class[] sig= new Class[] {Gene.class, PrintStream.class, HashMap.class};
		Method[] m= new Method[mNames.length];
		for (int i = 0; i < mNames.length; i++) {
			try {
				m[i]= target.getMethod(mNames[i], sig);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return m;
	}
	
	public static void _00_mainLoopChromosomes(String fName, Method[] m) {
		try {			
			GTFChrReader reader= new GTFChrReader(fName);
			reader.setChromosomeWise(true);
			reader.setReadGene(true);
			reader.read();
			Gene[] genes= reader.getGenes();
			HashMap localHash= new HashMap();
			while (genes!= null) {
				for (int i = 0; i < genes.length; i++) 
					for (int j = 0; j < m.length; j++) 
						m[j].invoke(null, new Object[] {genes[i], null, localHash});
				reader.read();
				genes= reader.getGenes();
			}			
			
//				output
			for (int i = 0; i < m.length; i++) {
				String baseName= fName+"_"+m[i].getName();
				m[i].invoke(null, new Object[]{null, baseName, localHash});
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * target method must comply with signature (Graph, Printstream)
	 * @param target
	 */
	public static void _00_mainLoop(File[] files, Method[] target) {
		
		String[] fList= new File("annotation"+ File.separator+ "ensembl").list();
		
		for (int i = 0; i < files.length; i++) {
			
			System.out.println(files[i]);
			try {
				GTFChrReader wrapper= new GTFChrReader(files[i].getAbsolutePath());
				if ((!files[i].getName().contains("_norman_"))&& (!wrapper.isApplicable()))
					wrapper.reformatFile();
				Gene[] ge= wrapper.read();
				HashMap localHash= new HashMap();
				while (ge!= null) {
					for (int j = 0; j < ge.length; j++) {
						for (int k = 0; k < target.length; k++) {
							try {
								target[k].invoke(null, new Object[] {ge[j], null, localHash});
							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}
					ge= wrapper.read();
				} 
				for (int k = 0; k < target.length; k++) {
					try {
						PrintStream p= new PrintStream(new FileOutputStream(files[i].getAbsolutePath()+"."+target[k].getName(), true));
						target[k].invoke(null, new Object[] {null, p, localHash});
						p.flush();
						p.close();
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
					
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	/**
				 * target method must comply with signature (Graph, Printstream)
				 * @param target
				 */
				public static void _00_mainLoopHuman(Method[] target) {
					
					File[] files= new File[] {
							new File(getGencodeAnnotation()),
							new File(Constants.getLatestUCSCAnnotation("human", "RefSeq", null)),
							new File(Constants.getLatestUCSCAnnotation("human", "Ensembl", null))
					};
					_00_mainLoop(files, target);
				}

	/**
	 * target method must comply with signature (Graph, Printstream)
	 * @param target
	 */
	public static void _00_mainLoopEnsemblChrom(String[] speNames, String startSpe, Method[] target) {
		
		String[] fList= new File("annotation"+ File.separator+ "ensembl").list();
		boolean skip= false;
		if (startSpe!= null)
			skip= true;
		
		for (int i = 0; i < speNames.length; i++) {
			if (skip&& !speNames[i].equals(startSpe))
				continue;
			skip= false;
			
			System.out.println(Species.SP_NAMES_METAZOA[i]);
			try {
				String iname= null;
				for (int x = 0; x < fList.length; x++) {
					if (fList[x].startsWith(Species.SP_NAMES_METAZOA[i])&&
							fList[x].endsWith(".gff")&&
							fList[x].contains("_norman_")) {
						iname= "annotation"+ File.separator+ "ensembl"+ File.separator+ fList[x];
						break;
					}
				}
				for (int x = 0; iname== null&& x < fList.length; x++) {
					if (fList[x].startsWith(Species.SP_NAMES_METAZOA[i])&&
							fList[x].endsWith(".gff")) {
						iname= "annotation"+ File.separator+ "ensembl"+ File.separator+ fList[x];
						break;
					}
				}
				if (iname== null) {
					System.err.println("Annotation for "+speNames[i]+" not found");
					continue;
				}
					
				
				GTFChrReader wrapper= new GTFChrReader(iname);
				if ((!iname.contains("_norman_"))&& (!wrapper.isApplicable()))
					wrapper.reformatFile();
				Gene[] ge= wrapper.read();
				HashMap localHash= new HashMap();
				while (ge!= null) {
					for (int j = 0; j < ge.length; j++) {
						for (int k = 0; k < target.length; k++) {
							try {
								target[k].invoke(null, new Object[] {ge[j], null, localHash});
							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}
					ge= wrapper.read();
				} 
				for (int k = 0; k < target.length; k++) {
					try {
						PrintStream p= new PrintStream(new FileOutputStream(iname+"."+target[k].getName(), true));
						target[k].invoke(null, new Object[] {null, p, localHash});
						p.flush();
						p.close();
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
					
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
