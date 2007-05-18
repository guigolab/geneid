package gphase;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Vector;

import gphase.algo.ASAnalyzer;
import gphase.algo.Analyzer;
import gphase.db.MapTable;
import gphase.io.TabDelimitedFormatWrapper;
import gphase.io.gtf.EncodeWrapper;
import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.ProgressiveIOWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.tools.Arrays;

public class Elza {
	public static void _01_allEventOut() {
		String comSpeName = "human";
		String[] specAnnos = Species.SPECIFIC_ANNOTATIONS[Species
				.getSpeciesNumber(comSpeName)];
		for (int i = specAnnos.length - 1; i >= 0; --i) {
			System.out.println(specAnnos[i]);
			String fName = Constants.getLatestUCSCAnnotation(comSpeName,
					specAnnos[i], null);
			EncodeWrapper gtf = new EncodeWrapper(fName);
			Graph g = gtf.getGraph(false);
			ASVariation[][] asVars = g
					.getASVariations(ASMultiVariation.FILTER_NONE);
			Arrays.sort2DFieldRev(asVars);

			try {
				PrintStream p = new PrintStream(new File("elza"
						+ File.separator + comSpeName + Constants.separator
						+ specAnnos[i] + Constants.separator + "allEvents"));
				ASAnalyzer.outputVariationCoordinates(asVars, p);
				p.flush();
				p.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

	}

	public static void main(String[] args) {
		// _01_allEventOut();
		// _01_allEventOut_new();
		 //_01_allEventOut_asd();
		// MapTable.addPDBInfo();
		_02_landscapePDB();
	}

	public static void _01_allEventOut_asd() {
		String comSpeName = "human";
		String[] specAnnos = Species.SPECIFIC_ANNOTATIONS[Species
				.getSpeciesNumber(comSpeName)];

		GTFChrReader gtf = new GTFChrReader(ASAnalyzer.INPUT_ASD);
		// gtf.reformatFile();
		gtf.transcript_id_tag = "Parent";
		Vector v = Analyzer.getASVariationCoordinates(gtf,
				ASMultiVariation.FILTER_NONE, -1);
		try {
			PrintStream p = new PrintStream(new File("elza" + File.separator
					+ comSpeName + Constants.separator + gtf.getFileName()
					+ Constants.separator + "allEvents"));
			for (int j = 0; j < v.size(); j++)
				p.println(v.elementAt(j));
			p.flush();
			p.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static void _01_allEventOut_new() {
		String comSpeName = "human";
		String[] specAnnos = Species.SPECIFIC_ANNOTATIONS[Species
				.getSpeciesNumber(comSpeName)];

		for (int i = specAnnos.length - 1; i >= 0; --i) {
			System.out.println(specAnnos[i]);
			String fName = Constants.getLatestUCSCAnnotation(comSpeName,
					specAnnos[i], null);

			GTFChrReader gtf = new GTFChrReader(fName);
			Vector v = Analyzer.getASVariationCoordinates(gtf,
					ASMultiVariation.FILTER_NONE, -1);
			try {
				PrintStream p = new PrintStream(new File("elza"
						+ File.separator + comSpeName + Constants.separator
						+ specAnnos[i] + Constants.separator + "allEvents"));
				for (int j = 0; j < v.size(); j++)
					p.println(v.elementAt(j));
				p.flush();
				p.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

	}

	public static void _02_landscapePDB() {

		String commonSpeName = "human";
		HashMap pdbMap= new HashMap();
		
		// read in additional maptable
		try {
			TabDelimitedFormatWrapper reader = new TabDelimitedFormatWrapper(
					new File("maptables" + File.separator + "elza.maptable")
							.getAbsolutePath());
			reader.read();
			String[] ids = null, pdb = null;
			pdb = reader.getColumn(0);
			ids = reader.getColumn(1);
			for (int i = 1; i < ids.length; i++) {	// skip first line
				ids[i] = ids[i].substring(ids[i].indexOf("via") + 4, ids[i]
						.length());
				pdbMap.put(pdb[i], ids[i]);
			}
			System.out.println("read " + ids.length + " additional ids.");
		} catch (Exception e) {
			e.printStackTrace();
		}

			// readin IDs
		String[] ids= null;
		try {
			TabDelimitedFormatWrapper reader = new TabDelimitedFormatWrapper(
					new File("elza" + File.separator + "pdbs_nrdb2_all_nr")
							.getAbsolutePath());
			reader.read();
			ids = reader.getColumn(0);
			Vector v= new Vector();
			int cntMap= 0;
			for (int i = 0; i < ids.length; i++) {
				String pdb= ids[i].toUpperCase();
				Object o= pdbMap.get(pdb);
				if (o== null)
					v.add(pdb);
				else {
					v.add(o);
					++cntMap;
				}
			}
			System.out.println("query "+ids.length+" ("+cntMap+" mappable)");
			ids= (String[]) Arrays.toField(v);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
			// get genes
		String[] notFoundIDs = ids;
		Vector v= new Vector();
		String[] specAnnos = Species.SPECIFIC_ANNOTATIONS[Species
				.getSpeciesNumber(commonSpeName)];
		String subdir= "elza";
		for (int k = 0; specAnnos != null && notFoundIDs != null
				&& notFoundIDs.length > 0 && k < specAnnos.length; k++) {

			String[] specIDs = MapTable.getPrimaryGeneID("human", specAnnos[k],
					notFoundIDs);
			notFoundIDs = MapTable.getNotFoundIDs();

			String fName = Constants.getLatestUCSCAnnotation(commonSpeName,
					specAnnos[k], null);
			EncodeWrapper gtf = new EncodeWrapper(fName);
			try {
				// PrintStream err= System.err;
				// System.setErr(new PrintStream("delme.tmp"));
				gtf.read(specIDs);
				// new File("delme.tmp").delete();
				// System.setErr(err);
			} catch (Exception e) {
				e.printStackTrace();
			}
			GTFObject[] obj = gtf.getGtfObj();
			v = (Vector) gphase.tools.Arrays.addAll(v, obj);
			int x = 0;
			if (specIDs != null)
				x = specIDs.length;
			int y = 0;
			if (obj != null)
				y = obj.length;
			System.out.println("Retrieved " + x + " new gene ids, " + y
					+ " new gtf objects. Total " + v.size());
		}
		if (notFoundIDs != null && notFoundIDs.length > 0) {
			System.out.print("Not found: ");
			for (int k = 0; k < notFoundIDs.length; k++)
				System.out.print(notFoundIDs[k] + " ");
			System.out.println();
		}

		// get reference gene set
		GTFObject[] objs = (GTFObject[]) gphase.tools.Arrays.toField(v);
		v= null;
		System.gc();
		System.out.println("analyzing " + objs.length + " gtf objects..");
		Gene[] ge = ProgressiveIOWrapper.assemble(objs);
		objs= null;
		System.gc();
		HashMap chrHash = new HashMap();
		for (int k = 0; k < ge.length; k++) {
			Vector geV = (Vector) chrHash.get(ge[k].getChromosome());
			if (geV == null) {
				geV = new Vector();
				chrHash.put(ge[k].getChromosome(), geV);
			}
			geV.add(ge[k]);
		}

		// overlap with reference annotations
		for (int k = 0; specAnnos != null && k < specAnnos.length; k++) {
			ProgressiveIOWrapper gtf = new ProgressiveIOWrapper(Constants
					.getLatestUCSCAnnotation(commonSpeName, specAnnos[k]));
			try {
				gtf.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			Species sp = gtf.assemble();
			ge = sp.getGenes();
			System.out.println("Analyzing reference annotation " + specAnnos[k]
					+ ", " + ge.length + " genes.");

			File dir = new File(subdir + File.separator + commonSpeName + "_"
					+ specAnnos[k] + "_landscape"+File.separator+".");
			if (dir.exists())
				dir.delete();
			dir.mkdir();

			System.out.print("retrieving variations..");
			System.out.flush();
			ASVariation[][] vars = sp
					.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			gphase.tools.Arrays.sort2DFieldRev(vars);
			Structurator.writeHTML(vars, dir.getAbsolutePath());
			System.out.println("done");

			Vector geneV = new Vector(ge.length / 2);
			for (int m = 0; m < ge.length; m++) {
				Vector geV = (Vector) chrHash.get(ge[m].getChromosome());
				if (geV == null) {
					sp.remove(ge[m], true);
					continue;
				}
				int n = 0;
				for (n = 0; n < geV.size(); n++)
					if (ge[m].overlaps((Gene) geV.elementAt(n)))
						break;
				if (n == geV.size())
					sp.remove(ge[m], true);
			}
			ge = sp.getGenes();
			System.out.println("Found " + ge.length + " genes.");

			dir = new File(subdir + File.separator + "allPDB_"
					+ specAnnos[k] + "_landscape"+File.separator+".");
			if (dir.exists())
				dir.delete();
			dir.mkdir();

			vars = sp.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			gphase.tools.Arrays.sort2DFieldRev(vars);
			Structurator.writeHTML(vars, dir.getAbsolutePath());
		}

	}
}
