package gphase;

import gphase.algo.ASAnalyzer;
import gphase.io.gtf.EncodeWrapper;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.GTFWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.AbstractRegion;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.EncodeRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.model.Translation;
import gphase.tools.ENCODE;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JFrame;

import org.freehep.graphicsio.emf.SetStretchBltMode;

import qalign2.algo.sequence.DistanceMatrixModel;
import qalign2.algo.tree.NJTreeModel;
import qalign2.algo.tree.NJWrapper;
import qalign2.gui.powerPane.treeView.TreePanel;
import qalign2.tree.RootedTree;

public class SpliceSiteConservationComparator {

	final static String ALI_DIR = "encode/msa/MFA/";

	String encRegion;

	String[] names;

	int[] pos;
	
	int[][][] ssNr= new int[2][][]; 	// AD, AA, CD, CA
	int[][][][][] ssColNr= new int[2][][][][]; // UTR/CDS - AD, AA, CD, CA, All - genes - SSs/1 - sites/column 

	int delta;

	static void output5UTRspliced(boolean alternative) {
		Graph g = getGraph();
		g.filterNonCodingTranscripts();

		Gene[] ge = g.getGenes();
		Vector w = new Vector(); // constitutive
		Vector u = new Vector(); // alternative
		for (int i = 0; i < ge.length; i++) {
			ASVariation[] as = ge[i]
					.getASVariations(ASMultiVariation.FILTER_NONE);
			Vector v = new Vector();
			for (int k = 0; as != null && k < as.length; k++) { // collect AS
																// transcripts
																// in UTR
				if (as[k].isTouching5UTR()) {
					int m;
					for (m = 0; m < v.size(); m++)
						if (v.elementAt(m) == as[k].getTranscript1())
							break;
					if (m == v.size())
						v.add(as[k].getTranscript1());
					for (m = 0; m < v.size(); m++)
						if (v.elementAt(m) == as[k].getTranscript2())
							break;
					if (m == v.size())
						v.add(as[k].getTranscript2());
				}

			}
			Transcript[] tc = ge[i].getTranscripts(); // get constitutive UTR
														// transcripts
			for (int j = 0; j < tc.length; j++) {
				int k;
				for (k = 0; k < v.size(); k++)
					if (v.elementAt(k) == tc[j])
						break;
				if (k == v.size())
					w.add(tc[j]);
				else
					u.add(tc[j]);
			}
		}

		try {
			String id;
			Vector y;
			if (alternative) {
				id = "alt";
				y = u;
			} else {
				id = "const";
				y = w;
			}
			BufferedWriter buffy = new BufferedWriter(new FileWriter("5UTR_"
					+ id + "_spliced"));
			for (int i = 0; i < y.size(); i++) {
				boolean skip = false;
				Transcript t = (Transcript) y.elementAt(i);
				DirectedRegion[] reg = t.get5UTRRegion(false);
				if (reg == null)
					continue;
				String region = "";
				int min = 0, max = 0;
				for (int j = 0; j < reg.length; j++) {
					if (j == 0)
						min = Math.abs(reg[j].get5PrimeEdge());
					if (j == reg.length - 1)
						max = Math.abs(reg[j].get3PrimeEdge());
					DirectedRegion regSav = (DirectedRegion) reg[j].clone();
					reg[j] = ENCODE.convertToEncodeCoord(reg[j]);
					if (reg[j] == null) {
						System.err
								.println("outside encode, skipping transcript "
										+ t.getTranscriptID());
						skip = true;
						break;
					}
					String s = getSubstring(reg[j].getID(), "human", Math
							.abs(reg[j].getStart()),
							Math.abs(reg[j].getEnd()) + 1); // end excl
					if (!reg[j].isForward())
						s = gphase.tools.Arrays.reverseComplement(s);
					region += s;
				}
				if (skip)
					continue;
				if (min > max) {
					int h = min;
					min = max;
					max = h;
				}
				String ss = ">" + t.getChromosome() + "_" + min + "_" + max
						+ "_" + t.getStrand() + " " + t.getTranscriptID()
						+ "\n" + region + "\n";
				buffy.write(ss);
			}
			buffy.flush();
			buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	// TODO: check redundancy!
	static void atgPhase(PrintStream pr) {
		Graph g = getGraph();
		g.filterNonCodingTranscripts();

		Gene[] ge = g.getGenes();
		Vector u = new Vector();
		Vector v = new Vector();
		int c0 = 0, c1 = 0, c2 = 0, c3 = 0;
		for (int i = 0; i < ge.length; i++) {

			// get AS events in 5UTR
			ASVariation[] as = ge[i]
					.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			for (int k = 0; as != null && k < as.length; k++) { // collect AS
																// transcripts
																// in UTR
				if (as[k].isCompletelyIn5UTR()) {
					u.add(as[k]);
				}
			}

			// get const 5UTR
			// Transcript[] tr= ge[i].getTranscripts();
			// Vector tmp= new Vector();
			// for (int j = 0; j < tr.length; j++) {
			// DirectedRegion utr= tr[j].get5UTRRegion(true)[0];
			// SpliceSite[] ss= ge[i].getSpliceSites();
			// boolean not= false;
			// for (int k = 0; k < ss.length; k++)
			// if (utr.contains(ss[k].getPos()))
			// if (!ss[k].isConstitutive()) {
			// not= true;
			// break;
			// }
			// if(!not)
			// tmp.add(tr);
			// }
			// if (tmp.size()> 0) { // find longest utr
			// int max= 0;
			// int maxNr= -1;
			// for (int j = 0; j < tmp.size(); j++)
			// if (((Transcript) tmp.elementAt(j)).getLength()> max) {
			// max= ((Transcript) tmp.elementAt(j)).getLength();
			// maxNr= j;
			// }
			// v.add(tmp.elementAt(maxNr));
			// }
		}

		int offsetDown = 10;
		int offsetUp = 10;
		Vector chkTrans = new Vector();
		HashMap redHash = new HashMap();
		for (int i = 0; i < u.size(); i++) {
			ASVariation as = (ASVariation) u.elementAt(i);
			DirectedRegion[] varRegs = as.getVariableRegion();
			for (int j = 0; j < varRegs.length; j++) {
				DirectedRegion regSav = (DirectedRegion) varRegs[j].clone();
				varRegs[j] = ENCODE.convertToEncodeCoord(varRegs[j]);
				String s = getSubstring(varRegs[j].getID(), "human", Math
						.abs(varRegs[j].getStart()), Math.abs(varRegs[j]
						.getEnd()));
				// if (!as.getGene().isForward())
				// s= Arrays.reverseComplement(s);
				s = s.toUpperCase();

				int index = -1;
				int ind2 = 0;
				String q1 = "ATG";
				String[] q2 = new String[] { "TGA", "TAG", "TAA" };
				if (!as.getGene().isForward()) {
					q1 = "CAT";
					// q2= new String[] {"TCA", "CTA", "TTA"};
				}

				int c0Sav = c0;
				while ((ind2 = s.indexOf(q1, (index + 1))) > index) {

					// retrieve genomic position (directed)
					int genPos = ind2 + Math.abs(regSav.getStart());
					if (!as.getGene().isForward())
						genPos += 2; // not the C of cat..
					if (!as.getGene().isForward())
						genPos = -genPos;

					// check redundancy
					String id = as.getTranscript1().getTranscriptID()
							.substring(
									0,
									as.getTranscript1().getTranscriptID()
											.lastIndexOf("-"))
							+ "-" + genPos;
					Object o = redHash.get(id);
					if ((redHash.get(id) != null)
							&& ((Integer) redHash.get(id)).intValue() == Math
									.abs(genPos)) {
						index = ind2;
						continue;
					} else
						redHash.put(id, new Integer(Math.abs(genPos)));

					// check phase
					Transcript tr = null;
					if (as.getTranscript1().isExonic(genPos)) {
						if (as.getTranscript2().isExonic(genPos))
							System.err
									.println("assertion for as region failed");
						else
							tr = as.getTranscript1();
					} else {
						if (as.getTranscript2().isExonic(genPos))
							tr = as.getTranscript2();
						else
							System.err
									.println("assertion for as region failed, null");
					}
					int dist = tr.getDistFromATG(genPos);
					int mod = dist % 3;

					// cross-check for real atg
					Transcript[] trr = tr.getGene().getTranscripts();
					boolean skip = false;
					for (int k = 0; k < trr.length; k++)
						if (trr[k].getTranslations()[0].get5PrimeEdge() == genPos) {
							skip = true;
							break;
						}
					if (skip) {
						index = ind2;
						continue; // filter no real new atg (eg
									// RP11-126K1.4-003)
					}

					// look for stops
					String utr5 = tr.getSequence(Transcript.REGION_5UTR);
					int ph0 = 0;
					int phX = 0;
					int exOffset = tr.getExonicPosition(regSav.get5PrimeEdge());
					if (tr.isForward())
						exOffset += ind2;
					else
						exOffset += s.length() - ind2 - 2; // 2 for start of
															// ATG
					for (int k = 0; k < q2.length; k++) {
						int idx1 = exOffset, idx2 = 0;
						while ((idx2 = utr5.indexOf(q2[k], (idx1 + 1))) > idx1) {
							int xx = idx2 - exOffset;
							if ((xx) % 3 == 0)
								ph0++;
							else
								phX++;
							idx1 = idx2;
						}
					}

					// write
					pr.println(regSav.getChromosome() + " "
							+ as.getGene().getStrand() + " "
							+ tr.getTranscriptID() + " " + Math.abs(genPos)
							+ " " + mod + " " + ph0 + " " + phX + " "
							+ utr5.substring(exOffset));
					c0++;
					if (mod == 0) {
						c1++;
						if (ph0 == 0)
							c2++;
					} else if (ph0 == 1)
						c3++;

					int gPos = Math.abs(varRegs[j].get5PrimeEdge());
					if (tr.isForward())
						gPos += ind2;
					else
						gPos += s.length() - ind2 - 2; // 2 for start of ATG
					if (mod == 0 && ph0 == 0)
						System.out.println(regSav.getChromosome() + " "
								+ as.getGene().getStrand() + " "
								+ tr.getTranscriptID() + " " + Math.abs(genPos)
								+ " " + varRegs[j].getID() + " " + gPos + " "
								+ mod + " " + ph0 + " " + phX + " "
								+ utr5.substring(exOffset));

					// update checked transcripts
					int k;
					for (k = 0; k < chkTrans.size(); k++)
						if (chkTrans.elementAt(k) == tr)
							break;
					if (k == chkTrans.size())
						chkTrans.add(tr);

					index = ind2;
				}
			}
		}

		// get alt trancripts wo atg
		Vector altTransWOalt = new Vector();
		for (int i = 0; i < u.size(); i++) {
			ASVariation as = (ASVariation) u.elementAt(i);

			Transcript tr = as.getTranscript1();
			int j;
			for (j = 0; j < chkTrans.size(); j++)
				if (tr == chkTrans.elementAt(j))
					break;
			if (j == chkTrans.size())
				altTransWOalt.add(tr);

			tr = as.getTranscript2();
			for (j = 0; j < chkTrans.size(); j++)
				if (tr == chkTrans.elementAt(j))
					break;
			if (j == chkTrans.size())
				altTransWOalt.add(tr);
		}

		// count const transcripts wo atg
		for (int i = 0; i < v.size(); i++) {

		}

		System.out.println(chkTrans.size() + " w 1+ ATG / "
				+ altTransWOalt.size() + " wo ATG in AS transcripts");
		System.out.println("ATG total: " + c0 + ", inPh: " + c1
				+ ", inPh+noStop: " + c2 + ", outPh+stop:" + c3);
	}

	public static void proteomeVSutrDiversity() {

		Graph g = getGraph();
		g.filterNonCodingTranscripts();
		g.filterSingleTranscriptGenes();

		Gene[] ge = g.getGenes();
		float[] cdsVar = new float[ge.length];
		float[] utrVar = new float[ge.length];
		float cdstt = 0f, utrtt = 0f;
		int ctr = 0;
		for (int i = 0; i < ge.length; i++) {
			ASVariation[] vars = ge[i]
					.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			if (vars == null)
				continue;
			DirectedRegion totcds = ge[i].getRealCDS();
			DirectedRegion totutr5 = ge[i].getReal5UTR();
			int cdsDif = 0, utr5Dif = 0;
			for (int j = 0; j < vars.length; j++) {
				if (totcds.contains(vars[j].getRegion()))
					cdsDif += vars[j].getLengthDiff(true);
				else if (totutr5.contains(vars[j].getRegion()))
					utr5Dif += vars[j].getLengthDiff(true);
			}
			cdsVar[i] = (float) cdsDif / totcds.getLength();
			utrVar[i] = (float) utr5Dif / totutr5.getLength();
			int a = utr5Dif;
			int b = totutr5.getLength();
			cdstt += cdsVar[i];
			if (!Float.isNaN(utrVar[i]))
				utrtt += utrVar[i];
			System.out.println(cdsVar[i] + ",\t" + utrVar[i] + "\t\t" + cdstt
					+ "," + utrtt);
			++ctr;
		}
		System.out.println("==> " + cdstt + "," + utrtt + "\t\t"
				+ (cdstt / ctr) + "," + (utrtt / ctr));

	}

	public static void proteomeVSutrDiversity2() {

		Graph g = getGraph();
		g.filterNonCodingTranscripts();
		g.filterSingleTranscriptGenes();

		Gene[] ge = g.getGenes();
		float[] cdsVar = new float[ge.length];
		float[] utrVar = new float[ge.length];
		float cdstt = 0f, utrtt = 0f;
		int ctr = 0;
		Comparator compi = new ASVariation.StructureComparator();
		System.out.println("reg\tcommon\tdiff");
		for (int i = 0; i < ge.length; i++) {
			ASVariation[] vars = ge[i]
					.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			if (vars == null)
				continue;
			HashMap usedVars = new HashMap();
			for (int j = 0; j < vars.length; j++) {
				if (vars[j].isPartiallyCoding())
					continue;

				int diff = vars[j].getDiffLength();
				// int common= vars[j].getCommonLength();

				String id = "";
				int base = 0;
				if (vars[j].isProteinCoding()) {
					id = "prot";
					// int b1= vars[j].getTranscript1().getCDSLength(true);
					// int b2= vars[j].getTranscript2().getCDSLength(true);
					// base= Math.min(b1,b2);
					DirectedRegion[] dir = DirectedRegion.unite(new DirectedRegion[][] {vars[j]
							.getTranscript1().getCDSRegions(), vars[j]
							.getTranscript2().getCDSRegions()});
					for (int k = 0; k < dir.length; k++)
						base += dir[k].getLength();
				} else if (vars[j].isNotAtAllCoding()) {
					if (!vars[j].isCompletelyIn5UTR())
						continue;
					id = "utr";
					// int b1= vars[j].getTranscript1().get5UTRLength(true);
					// int b2= vars[j].getTranscript2().get5UTRLength(true);
					// base= Math.min(b1,b2);
					DirectedRegion[] dir = DirectedRegion.unite(new DirectedRegion[][] {vars[j]
							.getTranscript1().get5UTRRegion(false), vars[j]
							.getTranscript2().get5UTRRegion(false)});
					for (int k = 0; k < dir.length; k++)
						base += dir[k].getLength();
				}
				System.out.println(id + "\t" + base + "\t" + diff);
			}
		}
	}

	static private SpliceSite[] extractSS(ASVariation[] vars) {
		Comparator compi= new SpliceSite.PositionComparator();
		SpliceSite[] ss= new SpliceSite[0];
		for (int j = 0; j < vars.length; j++) {
			SpliceSite[] su= vars[j].getSpliceUniverse();
			for (int k = 0; k < su.length; k++) {	// retrieve SSs
				int p= Arrays.binarySearch(ss, su[k], compi);
				if (p< 0)
					ss= (SpliceSite[]) gphase.tools.Arrays.insert(ss, su[k], p);	// every SS only once
			}
		}
		return ss;
	}

	public static void fullATGAnalysis() {

		Graph g = getGraph();
		g.filterNonCodingTranscripts();

		Vector v = new Vector();
		v = getATG_UTR_alt(g, v);
		// v= getATG_UTR(g, v);
		// v= getATG_Start(g, v);

		GTFObject[] gtfObj = (GTFObject[]) gphase.tools.Arrays.toField(v);
		GTFWrapper gtfWrap = new GTFWrapper(new File("ATG_analysis.gtf")
				.getAbsolutePath());
		gtfWrap.setGtfObj(gtfObj);
		try {
			gtfWrap.write();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	static Vector getATG_Start(Graph g, Vector gtfV) {
		Gene[] ge = g.getGenes();

		Vector u = new Vector();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] tr = ge[i].getTranscripts();
			for (int j = 0; j < tr.length; j++) {

				// get translation start
				Translation tl = tr[j].getTranslations()[0];
				if (tl == null)
					continue;
				int genPos = tl.get5PrimeEdge();

				// get sequence
				String exSeq = tr[j]
						.getSequence(Transcript.REGION_COMPLETE_TRANSCRIPT);
				if (exSeq == null) // outside encode
					continue;
				// determine boundaries, sequence
				int offsetDown = 30;
				int offsetUp = 30;

				int offset = tr[j].getExonicPosition(genPos);
				String atg = exSeq.substring(offset, offset + 3);
				int deltaUp = 0, deltaDown = 0;
				if (offset - offsetDown < 0) {
					deltaDown = offsetDown - offset;
					offsetDown = offset;
				}
				if (offset + 3 + offsetUp > exSeq.length()) {
					deltaUp = offset + 3 + offsetUp - exSeq.length();
					offsetUp = exSeq.length() - (offset + 3);
				}
				String seq = exSeq.substring(offset - offsetDown, offset + 3
						+ offsetUp); // 3 for ATG
				for (int k = 0; k < deltaDown; k++)
					seq = "-" + seq;
				for (int k = 0; k < deltaUp; k++)
					seq += "-";

				// build gtf obj
				int min = genPos;
				int max = Math.abs(genPos + 2);
				min = Math.abs(min); // not before ;)
				if (min > max) {
					int h = min;
					min = max;
					max = h;
				}
				GTFObject gtf = new GTFObject();
				gtf.setSeqname("chr" + ge[i].getChromosome());
				gtf.setSource("encodehg17");
				try {
					gtf.setFeature("ATG");
				} catch (Exception e) {
					e.printStackTrace();
				}
				gtf.setStart(min);
				gtf.setEnd(max);
				try {
					gtf.setStrand((genPos > 0) ? "+" : "-");
				} catch (Exception e) {
					e.printStackTrace();
				}
				gtf.addAttribute("type", "ATG_start");
				gtf.addAttribute("seq", seq);
				gtf.addAttribute("transcript", tr[j].getTranscriptID());

				// check redundancy
				int k;
				for (k = 0; k < u.size(); k++)
					if (u.elementAt(k).equals(gtf))
						break;
				if (k == u.size())
					u.add(gtf);
			}
		}

		for (int i = 0; i < u.size(); i++)
			gtfV.add(u.elementAt(i));
		return gtfV;
	}

	static Vector getATG_UTR(Graph g, Vector gtfV) {
		Gene[] ge = g.getGenes();

		Vector u = new Vector();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] tr = ge[i].getTranscripts();
			for (int j = 0; j < tr.length; j++) {

				// get sequence
				String utr5 = tr[j].getSequence(Transcript.REGION_5UTR);
				if (utr5 == null) // outside encode
					continue;
				int index = -1;
				int ind2 = 0;
				String q1 = "ATG";
				while ((ind2 = utr5.indexOf(q1, (index + 1))) > index) {

					// determine boundaries, sequence
					int offsetDown = 30;
					int offsetUp = 30;
					int offset = ind2;
					String atg = utr5.substring(offset, offset + 3);
					int deltaUp = 0, deltaDown = 0;
					if (offset - offsetDown < 0) {
						deltaDown = offsetDown - offset;
						offsetDown = offset;
					}
					if (offset + 3 + offsetUp > utr5.length()) {
						deltaUp = offset + 3 + offsetUp - utr5.length();
						offsetUp = utr5.length() - (offset + 3);
					}
					String seq = utr5.substring(offset - offsetDown, offset + 3
							+ offsetUp); // 3 for ATG
					for (int k = 0; k < deltaDown; k++)
						seq = "-" + seq;
					for (int k = 0; k < deltaUp; k++)
						seq += "-";

					// cross-check for real atg
					int min = tr[j].getGenomicPosition(ind2);
					Transcript[] trr = tr[j].getGene().getTranscripts();
					boolean skip = false;
					for (int k = 0; k < trr.length; k++)
						if (trr[k].getTranslations()[0].get5PrimeEdge() == min) {
							skip = true;
							break;
						}
					if (skip) {
						index = ind2;
						continue; // filter no real new atg (eg
									// RP11-126K1.4-003)
					}

					// build gtf obj
					if (min == tr[j].getTranslations()[0].get5PrimeEdge())
						continue; // filter no real new atg
					int max = Math.abs(min + 2);
					min = Math.abs(min); // not before max ;)
					if (min > max) {
						int h = min;
						min = max;
						max = h;
					}
					GTFObject gtf = new GTFObject();
					gtf.setSeqname("chr" + ge[i].getChromosome());
					gtf.setSource("encodehg17");
					try {
						gtf.setFeature("ATG");
					} catch (Exception e) {
						e.printStackTrace();
					}
					gtf.setStart(min);
					gtf.setEnd(max);
					try {
						gtf.setStrand((tr[j].isForward()) ? "+" : "-");
					} catch (Exception e) {
						e.printStackTrace();
					}
					gtf.addAttribute("type", "ATG_utr");
					gtf.addAttribute("seq", seq);
					gtf.addAttribute("transcript", tr[j].getTranscriptID());

					// check redundancy
					int k;
					for (k = 0; k < gtfV.size(); k++) {
						GTFObject o = (GTFObject) gtfV.elementAt(k);
						if (o.getStart() == gtf.getStart()
								&& o.getEnd() == gtf.getEnd())
							break;
					}
					if (k == gtfV.size())
						u.add(gtf);

					index = ind2;
				}

			}
		}

		for (int i = 0; i < u.size(); i++)
			gtfV.add(u.elementAt(i));
		return gtfV;
	}

	// chr1:148,111,718-148,111,797
	static Vector getATG_UTR_alt(Graph g, Vector gtfV) {

		Gene[] ge = g.getGenes();
		Vector u = new Vector();

		// get AS events in 5UTR
		for (int i = 0; i < ge.length; i++) {
			ASVariation[] as = ge[i]
					.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			for (int k = 0; as != null && k < as.length; k++) { // collect AS
																// transcripts
																// in UTR
				if (as[k].isCompletelyIn5UTR()) {
					u.add(as[k]);
				}
			}
		}

		for (int i = 0; i < u.size(); i++) {
			ASVariation as = (ASVariation) u.elementAt(i);
			DirectedRegion[] varRegs = as.getVariableRegion();
			for (int j = 0; j < varRegs.length; j++) {
				DirectedRegion regSav = (DirectedRegion) varRegs[j].clone();
				varRegs[j] = ENCODE.convertToEncodeCoord(varRegs[j]);
				String s = getSubstring(varRegs[j].getID(), "human", Math
						.abs(varRegs[j].getStart()), Math.abs(varRegs[j]
						.getEnd()));
				s = s.toUpperCase();

				int index = -1;
				int ind2 = 0;
				String q1 = "ATG";
				if (!as.getGene().isForward()) {
					q1 = "CAT";
				}
				while ((ind2 = s.indexOf(q1, (index + 1))) > index) {

					// retrieve genomic position (directed)
					int genPos = ind2 + Math.abs(regSav.getStart());
					if (!as.getGene().isForward())
						genPos += 2; // not the C of cat..
					if (!as.getGene().isForward())
						genPos = -genPos;

					// get corresponding UTR
					Transcript tr = null;
					if (as.getTranscript1().isExonic(genPos)) {
						if (as.getTranscript2().isExonic(genPos))
							System.err
									.println("assertion for as region failed");
						else
							tr = as.getTranscript1();
					} else {
						if (as.getTranscript2().isExonic(genPos))
							tr = as.getTranscript2();
						else
							System.err
									.println("assertion for as region failed, null");
					}

					String utr5 = tr.getSequence(Transcript.REGION_5UTR); // already
																			// reversed
																			// for
																			// neg
																			// strand

					// determine boundaries, sequence
					int offsetDown = 30;
					int offsetUp = 30;

					int offset = tr.getExonicPosition(regSav.get5PrimeEdge());
					if (tr.isForward())
						offset += ind2;
					else
						offset += s.length() - ind2 - 2; // 2 for start of
															// ATG
					String atg = utr5.substring(offset, offset + 3);
					int deltaUp = 0, deltaDown = 0;
					if (offset - offsetDown < 0) {
						deltaDown = offsetDown - offset;
						offsetDown = offset;
					}
					if (offset + 3 + offsetUp > utr5.length()) {
						deltaUp = offset + 3 + offsetUp - utr5.length();
						offsetUp = utr5.length() - (offset + 3);
					}
					String seq = utr5.substring(offset - offsetDown, offset + 3
							+ offsetUp); // 3 for ATG
					for (int k = 0; k < deltaDown; k++)
						seq = "-" + seq;
					for (int k = 0; k < deltaUp; k++)
						seq += "-";

					// cross-check for real atg
					Transcript[] trr = tr.getGene().getTranscripts();
					boolean skip = false;
					for (int k = 0; k < trr.length; k++)
						if (trr[k].getTranslations()[0].get5PrimeEdge() == genPos) {
							skip = true;
							break;
						}
					if (skip) {
						index = ind2;
						continue; // filter no real new atg (eg
									// RP11-126K1.4-003)
					}

					// build gtf obj
					int min = Math.abs(genPos);
					int max = Math.abs(genPos + 2);
					if (Math.abs(min) > Math.abs(max)) {
						int h = min;
						min = max;
						max = h;
					}
					GTFObject gtf = new GTFObject();
					gtf.setSeqname("chr" + regSav.getChromosome());
					gtf.setSource("encodehg17");
					try {
						gtf.setFeature("ATG");
					} catch (Exception e) {
						e.printStackTrace();
					}
					gtf.setStart(min);
					gtf.setEnd(max);
					try {
						gtf.setStrand((genPos > 0) ? "+" : "-");
					} catch (Exception e) {
						e.printStackTrace();
					}
					gtf.addAttribute("type", "ATG_UTR_alt");
					gtf.addAttribute("seq", seq);
					gtf.addAttribute("transcript", tr.getTranscriptID());

					// check redundancy
					int k;
					for (k = 0; k < gtfV.size(); k++)
						if (gtfV.elementAt(k).equals(gtf))
							break;
					if (k == gtfV.size())
						gtfV.add(gtf);

					index = ind2;
				}
			}
		}

		return gtfV;
	}

	static void ATGphase(boolean alternative) {
		Graph g = getGraph();
		g.filterNonCodingTranscripts();

		Gene[] ge = g.getGenes();
		Vector w = new Vector(); // constitutive
		Vector u = new Vector(); // alternative
		for (int i = 0; i < ge.length; i++) {
			ASVariation[] as = ge[i]
					.getASVariations(ASMultiVariation.FILTER_NONE);
			Vector v = new Vector();
			for (int k = 0; as != null && k < as.length; k++) { // collect AS
																// transcripts
																// in UTR
				if (as[k].isTouching5UTR()) {
					int m;
					for (m = 0; m < v.size(); m++)
						if (v.elementAt(m) == as[k].getTranscript1())
							break;
					if (m == v.size())
						v.add(as[k].getTranscript1());
					for (m = 0; m < v.size(); m++)
						if (v.elementAt(m) == as[k].getTranscript2())
							break;
					if (m == v.size())
						v.add(as[k].getTranscript2());
				}

			}
			Transcript[] tc = ge[i].getTranscripts(); // get constitutive UTR
														// transcripts
			for (int j = 0; j < tc.length; j++) {
				int k;
				for (k = 0; k < v.size(); k++)
					if (v.elementAt(k) == tc[j])
						break;
				if (k == v.size())
					w.add(tc[j]);
				else
					u.add(tc[j]);
			}
		}

		try {
			String id;
			Vector y;
			if (alternative) {
				id = "alt";
				y = u;
			} else {
				id = "const";
				y = w;
			}
			int ctrNull = 0, ctrInTot = 0, ctrOutTot = 0, ctrTot = 0, ctr1Tot = 0, ctr2Tot = 0;
			for (int i = 0; i < y.size(); i++) {
				boolean skip = false;
				Transcript t = (Transcript) y.elementAt(i);
				DirectedRegion[] reg = t.get5UTRRegion(false);
				String region = "";
				int min = 0, max = 0;
				for (int j = 0; j < reg.length; j++) {
					if (j == 0)
						min = Math.abs(reg[j].get5PrimeEdge());
					if (j == reg.length - 1)
						max = Math.abs(reg[j].get3PrimeEdge());
					DirectedRegion regSav = (DirectedRegion) reg[j].clone();
					reg[j] = ENCODE.convertToEncodeCoord(reg[j]);
					if (reg[j] == null) {
						// System.err.println("outside encode, skipping
						// transcript "+t.getTranscriptID());
						skip = true;
						break;
					}
					String s = getSubstring(reg[j].getID(), "human", Math
							.abs(reg[j].getStart()),
							Math.abs(reg[j].getEnd()) + 1); // end excl
					if (!reg[j].isForward())
						s = gphase.tools.Arrays.reverseComplement(s);
					region += s;
				}
				if (skip)
					continue;
				int index = -1;
				int ind2 = 0;
				int ctrIn = 0, ctrOut = 0, ctr1 = 0, ctr2 = 0, ctr = 0;
				while ((ind2 = region.indexOf("ATG", (index + 1))) > index) {
					int dist = region.length() - ind2;
					if (dist % 3 == 0)
						ctrIn++;
					else {
						ctrOut++;
						if (dist % 3 == 1)
							ctr1++;
						else
							ctr2++;
					}
					ctr++;
					index = ind2;
				}
				if (ctr > 0)
					System.out.println("(" + ctr + "):\t" + ctrIn + " in\t"
							+ ctrOut + " out (" + ctr1 + "," + ctr2 + ")");
				else
					ctrNull++;
				ctrTot += ctr;
				ctrInTot += ctrIn;
				ctrOutTot += ctrOut;
				ctr1Tot += ctr1;
				ctr2Tot += ctr2;
			}
			System.out.println("tot " + ctrTot + " in " + y.size() + " utrs");
			System.out.println("(" + ctrNull + ") without ATG: "
					+ ((float) ctrNull / (ctrNull + y.size())));
			System.out.println("(" + ctrInTot + ") in phase: "
					+ ((float) ctrInTot / ctrTot));
			System.out.println("(" + ctrOutTot + ") out phase: "
					+ ((float) ctrOutTot / ctrTot) + "("
					+ +((float) ctr1Tot / ctrTot) + ","
					+ ((float) ctr2Tot / ctrTot) + ")");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	static Graph getGraph() {
		// "encode/44regions_genes_CHR_coord.gtf"
		// "encode/gencode_races.gtf"
		// "encode/EnsemblGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_fromUCSC_inENCODEonly.gtf"
		// "encode/RefSeqGenes_fromUCSC.inENCODE.gtf"
		// "encode/RefSeqGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_all_fromENSEMBL.gtf"
		String fName = "encode/44regions_genes_CHR_coord.gtf";
		EncodeWrapper myWrapper = new EncodeWrapper(new File(fName)
				.getAbsolutePath()); // testGTF.gtf
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		boolean encode = false;
		if (fName.startsWith("encode/44regions_genes_CHR_coord"))
			encode = true;

		Graph g = myWrapper.getGraph(encode); // <===== check ENCODE here !!!
		return g;
	}

	static void checkATG() {
		// "encode/44regions_genes_CHR_coord.gtf"
		// "encode/gencode_races.gtf"
		// "encode/EnsemblGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_fromUCSC_inENCODEonly.gtf"
		// "encode/RefSeqGenes_fromUCSC.inENCODE.gtf"
		// "encode/RefSeqGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_all_fromENSEMBL.gtf"
		String fName = "encode/44regions_genes_CHR_coord.gtf";
		EncodeWrapper myWrapper = new EncodeWrapper(new File(fName)
				.getAbsolutePath()); // testGTF.gtf
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		boolean encode = false;
		if (fName.startsWith("encode/44regions_genes_CHR_coord"))
			encode = true;

		Graph g = myWrapper.getGraph(encode); // <===== check ENCODE here !!!
		g.filterNonCodingTranscripts();
		// g.filterCodingTranscripts();
		ASVariation[][] as = g
				.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
		int total = 0, pos = 0, totLen = 0;
		int[] distr = new int[30];
		for (int i = 0; i < distr.length; i++)
			distr[i] = 0;
		for (int i = 0; i < as.length; i++) {
			// if (!as[i][0].toString().equals("(1^ // 2^)"))
			if (!as[i][0].toString().equals("(1= // 2=)"))
				// if (!as[i][0].toBitString().contains("DBA")&&
				// !as[i][0].toBitString().contains("DCA")) // AD
				// if (!as[i][0].toBitString().contains("ABD")&&
				// !as[i][0].toBitString().contains("ACD")) // AA
				// if (!as[i][0].toBitString().contains("DBD")&&
				// !as[i][0].toBitString().contains("DCD")) // IR
				continue;
			for (int j = 0; j < as[i].length; j++) {
				// if (!as[i][j].isCompletelyIn5UTR())
				if (!as[i][j].isCompletelyInCDS())
					continue;
				DirectedRegion[] varRegs = as[i][j].getVariableRegion();
				// DirectedRegion[] varRegs= as[i][j].getExonicOnlyRegion();
				if (varRegs == null)
					continue;
				DirectedRegion[] savVarRegs = new DirectedRegion[varRegs.length];
				for (int k = 0; k < savVarRegs.length; k++)
					savVarRegs[k] = (DirectedRegion) varRegs[k].clone();
				for (int k = 0; k < varRegs.length; k++) {
					varRegs[k] = ENCODE.convertToEncodeCoord(varRegs[k]);
					if (varRegs[k] == null)// outside encode
						continue;
					String varStr = getSubstring(varRegs[k].getID(), "human",
							Math.abs(varRegs[k].getStart()), Math
									.abs(varRegs[k].getEnd()) + 1); // end excl
					if (!as[i][j].getGene().isForward())
						varStr = gphase.tools.Arrays.reverseComplement(varStr);
					int ctr = 0, index = -1;
					int ind2 = 0;
					while ((ind2 = varStr.indexOf("ATG", (index + 1))) > index) {
						++ctr;
						index = ind2;
					}
					// System.out.println(ctr);
					try {
						++distr[ctr];
					} catch (ArrayIndexOutOfBoundsException e) {
						int diff = ctr - distr.length + 1;
						int[] distrNew = new int[distr.length + diff];
						for (int m = 0; m < distr.length; m++)
							distrNew[i] = distr[i];
						++distrNew[ctr];
						distr = distrNew;

					}
					if (ctr > 0)
						++pos;
					++total;
					totLen += varStr.length();
				}
			}
		}
		for (int i = 0; i < distr.length; i++)
			if (distr[i] != 0)
				System.out.println(i + "\t" + distr[i]);

		System.out.println(pos + "/" + total + "=" + ((float) pos / total));
		System.out.println(pos + "/" + totLen + "=" + ((float) pos / totLen));
	}

	static void checkATG2() {
		// "encode/44regions_genes_CHR_coord.gtf"
		// "encode/gencode_races.gtf"
		// "encode/EnsemblGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_fromUCSC_inENCODEonly.gtf"
		// "encode/RefSeqGenes_fromUCSC.inENCODE.gtf"
		// "encode/RefSeqGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_all_fromENSEMBL.gtf"
		String fName = "encode/44regions_genes_CHR_coord.gtf";
		EncodeWrapper myWrapper = new EncodeWrapper(new File(fName)
				.getAbsolutePath()); // testGTF.gtf
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		boolean encode = false;
		if (fName.startsWith("encode/44regions_genes_CHR_coord"))
			encode = true;

		Graph g = myWrapper.getGraph(encode); // <===== check ENCODE here !!!
		g.filterNonCodingTranscripts();
		// g.filterCodingTranscripts();
		ASVariation[][] as = g
				.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
		int asym = 0, tot = 0;
		Vector vAlt = new Vector();
		Vector vCon = new Vector();
		for (int i = 0; i < as.length; i++) {
			// if (!as[i][0].toString().equals("(1^ // 2^)"))
			// if (!as[i][0].toString().equals("(1= // 2=)"))
			// if (!as[i][0].toString().equals("( // 1^2=)"))
			// if (!as[i][0].toBitString().startsWith("DBA")&&
			// !as[i][0].toBitString().startsWith("DCA")) // AD
			// if (!as[i][0].toBitString().startsWith("DABD")&&
			// !as[i][0].toBitString().startsWith("DACD")) // AA
			// if (!as[i][0].toBitString().contains("DBD")&&
			// !as[i][0].toBitString().contains("DCD")) // IR
			// continue;
			for (int j = 0; j < as[i].length; j++) {
				if (!as[i][j].isCompletelyIn5UTR())
					// if (!as[i][j].isCompletelyInCDS())
					continue;
				DirectedRegion[] varRegs = as[i][j].getVariableRegion();
				DirectedRegion[] varRegs2 = as[i][j].getExonicOnlyRegion();
				if (varRegs == null)
					continue;
				int ctr1 = 0;
				if (varRegs != null && varRegs[0] != null) {
					varRegs[0] = ENCODE.convertToEncodeCoord(varRegs[0]);
					if (varRegs[0] == null)
						continue;
					String varStr1 = getSubstring(varRegs[0].getID(), "human",
							Math.abs(varRegs[0].getStart()), Math
									.abs(varRegs[0].getEnd()) + 1); // end excl
					if (!as[i][j].getGene().isForward())
						varStr1 = gphase.tools.Arrays.reverseComplement(varStr1);
					vAlt.add(varStr1);
					int index = -1;
					int ind2 = 0;
					while ((ind2 = varStr1.indexOf("ATG", (index + 1))) > index) {
						++ctr1;
						index = ind2;
					}
				}
				int ctr2 = 0;
				if (varRegs2 != null && varRegs2[0] != null) {
					varRegs2[0] = ENCODE.convertToEncodeCoord(varRegs2[0]);
					if (varRegs2[0] == null)
						continue;
					String varStr2 = getSubstring(varRegs[0].getID(), "human",
							Math.abs(varRegs2[0].getStart()), Math
									.abs(varRegs2[0].getEnd()) + 1); // end
																		// excl
					if (!as[i][j].getGene().isForward())
						varStr2 = gphase.tools.Arrays.reverseComplement(varStr2);
					int index = -1;
					int ind2 = 0;
					while ((ind2 = varStr2.indexOf("ATG", (index + 1))) > index) {
						++ctr2;
						index = ind2;
					}
				}
				System.out.println(ctr1 + "\t" + ctr2);
				// if ((ctr1== 0&& ctr2!= 0)|| (ctr1!= 0&& ctr2== 0))
				// if (ctr1!= 0&& ctr2== 0)
				if (ctr1 > 0)
					// if (ctr1> 0|| ctr2> 0)
					asym++;
				tot++;
			}
		}
		System.out.println(asym + "/" + tot);
		// for (int i = 0; i < vAlt.size(); i++) {
		// System.out.println(vAlt.elementAt(i));
		// }
	}

	static void checkATG1() {
		// "encode/44regions_genes_CHR_coord.gtf"
		// "encode/gencode_races.gtf"
		// "encode/EnsemblGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_fromUCSC_inENCODEonly.gtf"
		// "encode/RefSeqGenes_fromUCSC.inENCODE.gtf"
		// "encode/RefSeqGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_all_fromENSEMBL.gtf"
		String fName = "encode/44regions_genes_CHR_coord.gtf";
		EncodeWrapper myWrapper = new EncodeWrapper(new File(fName)
				.getAbsolutePath()); // testGTF.gtf
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		boolean encode = false;
		if (fName.startsWith("encode/44regions_genes_CHR_coord"))
			encode = true;

		Graph g = myWrapper.getGraph(encode); // <===== check ENCODE here !!!
		g.filterNonCodingTranscripts();
		// g.filterCodingTranscripts();
		ASVariation[][] as = g
				.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
		int total = 0, pos = 0, totLen = 0;
		int[] distr = new int[20];
		for (int i = 0; i < distr.length; i++)
			distr[i] = 0;
		for (int i = 0; i < as.length; i++) {
			for (int j = 0; j < as[i].length; j++) {
				if (!as[i][j].is5UTR())
					continue;
				AbstractRegion[][] varRegs = as[i][j].getVariableRegions();
				AbstractRegion[][] savVarRegs = new AbstractRegion[2][];
				for (int k = 0; k < savVarRegs.length; k++)
					// debug
					for (int m = 0; m < savVarRegs.length; m++)
						savVarRegs[k][m] = (AbstractRegion) varRegs[k].clone();
				for (int k = 0; k < varRegs[0].length; k++) {
					if (k == 0) {
						if (as[i][j].getSpliceChain1()[0].isDonor())
							varRegs[0][0].setStart(as[i][j].getTranscript1()
									.getExonicPosition(
											varRegs[0][0].getStart(), -1));
						else
							varRegs[0][0].setStart(as[i][j].getTranscript1()
									.getExonicPosition(
											varRegs[0][0].getStart(), -2));
					}
					varRegs[k] = ENCODE.convertToEncodeCoord(varRegs[0][k]);
					String varStr = getSubstring(varRegs[k].getID(), "human",
							Math.abs(varRegs[k].getStart()), Math
									.abs(varRegs[k].getEnd()) + 1); // end excl
					if (!as[i][j].getGene().isForward())
						varStr = Arrays.reverseComplement(varStr);
					int ctr = 0, index = -1;
					int ind2 = 0;
					while ((ind2 = varStr.indexOf("ATG", (index + 1))) > index) {
						++ctr;
						index = ind2;
					}
					// System.out.println(ctr);
					++distr[ctr];
					if (ctr > 0)
						++pos;
					++total;
					totLen += varStr.length();
				}
			}
		}
		for (int i = 0; i < distr.length; i++)
			if (distr[i] != 0)
				System.out.println(i + "\t" + distr[i]);

		System.out.println(pos + "/" + total + "=" + ((float) pos / total));
		System.out.println(pos + "/" + totLen + "=" + ((float) pos / totLen));
	}

	static void treeConstruct() {
		// "encode/44regions_genes_CHR_coord.gtf"
		// "encode/gencode_races.gtf"
		// "encode/EnsemblGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_fromUCSC_inENCODEonly.gtf"
		// "encode/RefSeqGenes_fromUCSC.inENCODE.gtf"
		// "encode/RefSeqGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_all_fromENSEMBL.gtf"
		// "encode/Sequences_mapped_HAVANA_136.gtf"
		String fName = "encode/44regions_genes_CHR_coord.gtf";
		EncodeWrapper myWrapper = new EncodeWrapper(new File(fName)
				.getAbsolutePath()); // testGTF.gtf
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		boolean encode = false;
		if (fName.startsWith("encode/44regions_genes_CHR_coord"))
			encode = true;

		Graph g = myWrapper.getGraph(encode); // <===== check ENCODE here !!!
		g.filterNonCodingTranscripts();

		// splice sites
		SpliceSite[][] ss2 = g.getSpliceSites(Gene.REGION_REAL_5UTR);
		String[] species = { "human", "chimp", "macaque", "baboon", "galago",
				"marmoset" };
		DefaultRegion[] eRegs = ENCODE.getEncodeRegions();
		SpliceSiteConservationComparator conserv = new SpliceSiteConservationComparator(
				null, species);
		double[][] pwDist = conserv.getPWDistances(ss2[0]);
//		for (int i = 0; i < pwDist.length; i++) {
//			for (int j = 0; j < pwDist.length; j++) {
//				System.out.print(pwDist[i][j]+" ");
//			}
//			System.out.println();
//		}
		DistanceMatrixModel dist = new DistanceMatrixModel(pwDist, species);
		NJWrapper nj = null;
		try {
			nj= new NJWrapper(dist, species);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		RootedTree r= nj.getNjTree();	
		
		TreePanel pan= new TreePanel(r, TreePanel.DENDO_TREE);
//		JFrame frame= new JFrame();
//		frame.setSize(300,300);
//		frame.getContentPane().add(pan);
//		frame.pack();
//		frame.setVisible(true);
	}

	public static void analyzeStops(Graph g) {
		
			// "analyze_stop_constSS_maxCDS_statistics.txt"
			// "analyze_stop_constSS_maxCDS_donors.txt"
			// "analyze_stop_constSS_maxCDS_acceptors.txt"
		PrintStream p= null;
		try {
			p= new PrintStream("analyze_stop_constSS_maxCDS_donors.txt");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		String[] stops= new String[] {"TAA", "TGA", "TAG"};
		g.filterNonCodingTranscripts();
		g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
		SpliceSite[] ss= g.getSpliceSites(SpliceSite.CONSTITUTIVE_SS, Gene.REGION_MAX_CDS);
		DefaultRegion[] reg= ENCODE.getEncodeRegions();
		int ctrStopDon= 0, ctrStopAcc= 0, ctrNoStopDon= 0, ctrNoStopAcc= 0;
		Vector seqV= new Vector(ss.length);
		for (int i = 0; i < reg.length; i++) {
			Vector v= new Vector();
			for (int j = 0; j < ss.length; j++) 
				if (reg[i].contains(ss[j]))
					v.add(ss[j]);
			SpliceSite[] regSS= (SpliceSite[]) gphase.tools.Arrays.toField(v);
				
			SpliceSiteConservationComparator comp= 
				new SpliceSiteConservationComparator(reg[i].getID(), new String[] {"human"});
			
			for (int j = 0; regSS!= null&& j < regSS.length; j++) {
				int pos1, pos2, delta= 10;
				if ((regSS[j].getGene().getStrand()> 0&& regSS[j].isDonor()) ||
						regSS[j].getGene().getStrand()< 0&& regSS[j].isAcceptor()) {
					pos1= Math.abs(regSS[j].getPos())- reg[i].getStart()+ 1- 2;
					pos2= pos1+ delta;
				} else {
					pos2= Math.abs(regSS[j].getPos())- reg[i].getStart()+ 2;	// d kill +1
					pos1= pos2- delta;
				}
				String s= comp.getSubstring("human", pos1, pos2);
				if (regSS[j].getGene().getStrand()< 0)
					s= gphase.tools.Arrays.reverseComplement(s);
				System.currentTimeMillis();
				
				Exon[] e= regSS[j].getExons();
				int k;
				for (k = 0; k < e.length; k++) {
						// 	skip nc SSs
					if (regSS[j].isDonor()&& !e[k].isCoding3Prime())
						continue;
					if (regSS[j].isAcceptor()&& !e[k].isCoding5Prime())
						continue;
					
					int fr= e[k].getFrame();	// for acceptors, this directly corresponds to frame of pos#2
					if (regSS[j].isDonor()) {
						fr+= (3- (Math.abs(e[k].get5PrimeEdge()- e[k].get5PrimeCDS())+ 1- 1)% 3);	// end frame, position#2 in string s
						fr= --fr< 0?3+ fr:fr;	// pos#2 -> pos#3
						fr= fr> 2?fr- 3:fr;
					} 
					
					int backNt= (3-fr== 3)?0:3-fr;	// frame -> nt before					
					int pos= 2- backNt;	// 3rd pos in string - Ns read before
					while (pos+ 2< s.length()) {
						String cod= s.substring(pos, pos+ 3);
						int m;
						for (m = 0; m < stops.length; m++) 
							if (stops[m].equalsIgnoreCase(cod))
								break;
						if (m< stops.length)
							break;
						pos+= 3;
					}
					if (pos+ 2< s.length())
						break;
				}
				if (k< e.length) {
					if (regSS[j].isDonor()) {
						++ctrStopDon; 
						seqV.add(s);
					} else {
						++ctrStopAcc;
					}
				} else {
					if (regSS[j].isDonor()) {
						++ctrNoStopDon; 
						seqV.add(s);
					} else {
						++ctrNoStopAcc;
					}
				}
			}
		}
		
//		p.println("Donor "+ctrStopDon+ " stop, "+ ctrNoStopDon+ " nostop");
//		p.println("Actor "+ctrStopAcc+ " stop, "+ ctrNoStopAcc+ " nostop");
		for (int i = 0; i < seqV.size(); i++) {
			p.println(">const_don_realCDS");
			p.println(seqV.elementAt(i));
		}
		
	}
	
	
	public static void outputRegions(Graph g, String fName) {
		Gene[] ge= g.getGenes();
		Vector v= new Vector(ge.length* 3);
		GTFObject o;
		DirectedRegion r;
		for (int i = 0; i < ge.length; i++) {
			r= ge[i].getReal5UTR();
			if (r!= null&& Math.abs(r.getStart())<= Math.abs(r.getEnd())) {
				o= new GTFObject();
				o.setSeqname(r.getChromosome());
				try {
					o.setFeature("5UTR");
					o.setStrand(ge[i].getStrand());
				} catch (Exception e) {e.printStackTrace(); }
				o.setSource(ge[i].getGeneID());
				o.setStart(Math.abs(r.getStart()));
				o.setEnd(Math.abs(r.getEnd()));
				o.addAttribute("5UTR", "real");
				v.add(o);
			}
			
			r= ge[i].getRealCDS();
			if (r!= null&& Math.abs(r.getStart())<= Math.abs(r.getEnd())) {
				o= new GTFObject();
				o.setSeqname(r.getChromosome());
				try {
					o.setFeature("CDS");
					o.setStrand(ge[i].getStrand());
				} catch (Exception e) {e.printStackTrace(); }
				o.setSource(ge[i].getGeneID());
				o.setStart(Math.abs(r.getStart()));
				o.setEnd(Math.abs(r.getEnd()));
				o.addAttribute("CDS", "real");
				v.add(o);
			}
		}
		
		fName= Toolbox.checkFileExists(fName);
		GTFWrapper wrapper= new GTFWrapper(fName);
		wrapper.setGtfObj((GTFObject[]) gphase.tools.Arrays.toField(v));
		try {wrapper.write();
		} catch (Exception e) {e.printStackTrace();}
	}
	
	public static void main(String[] args) {
		
		Graph g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);

		analyzeStops(g);
		
		// "conserv_analysis_full.txt"
//		String outName= "regions.gtf";
//		outputRegions(g, outName);
		
//		PrintStream p= null;
//		try {
//			p = new PrintStream(outName);
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
		
//		double[][] cons= calcConservation(encode, myWrapper.getGtfObj(), p);
//		gphase.tools.Arrays.output(cons, System.out);
		
		// calcConservationChk(encode, myWrapper.getGtfObj(), p);
		
			// calcConservation 27.10.06
		//calcConservation(encode, myWrapper.getGtfObj(), p);
		
		//checkBackground();
		
		// checkATG2();
		//checkSSConservation();
		// output5UTRspliced(false);
		// ATGphase(false);

			// for evolution
		//treeConstruct();

		// try {
		// PrintStream pr= new PrintStream("ATG_POSITIONS");
		// atgPhase(pr);
		// pr.flush(); pr.close();
		// } catch (Exception e) {
		// e.printStackTrace();
		// }
		//		
		// fullATGAnalysis();
		// proteomeVSutrDiversity2();
	}

	static void checkSSConservation() {
		// "outputSS_encode_coding.encode_coord.gtf"
		// "test.gtf"
		long t0 = System.currentTimeMillis();
		// "outputSS_encode_coding_NEW.encode_coord.gtf"
		String inputFile = "SSout_new.encode_coord.gtf";
		String[] species = { "human", "chimp", "macaque", "baboon", "galago",
				"marmoset" };
		int delta = 7;
		GTFWrapper gtf = new GTFWrapper(new File(inputFile).getAbsolutePath());
	
		// collect enc regions
		Vector v = new Vector();
		try {
			gtf.read();
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		GTFObject[] gtfObj = gtf.getGtfObj();
		for (int i = 0; i < gtfObj.length; i++) {
			int j;
			for (j = 0; j < v.size(); j++)
				if (gtfObj[i].getSeqname().equals(v.elementAt(j)))
					break;
			if (j == v.size())
				v.add(gtfObj[i].getSeqname());
		}
		String[] encReg = (String[]) gphase.tools.Arrays.toField(v);
	
		// iterate regions
		Vector resultV = new Vector();
		for (int i = 0; i < encReg.length; i++) {
			v = new Vector();
			for (int j = 0; j < gtfObj.length; j++)
				if (gtfObj[j].getSeqname().equalsIgnoreCase(encReg[i]))
					v.add(gtfObj[j]);
			GTFObject[] gtfEnc = (GTFObject[]) gphase.tools.Arrays.toField(v);
	
			SpliceSiteConservationComparator sp = new SpliceSiteConservationComparator(
					encReg[i], species);
			GTFObject[] tmpObj = sp.calcConservation(gtfEnc, delta);
			for (int j = 0; j < tmpObj.length; j++)
				resultV.add(tmpObj[j]);
		}
	
		gtfObj = (GTFObject[]) gphase.tools.Arrays.toField(resultV);
		gtf.setGtfObj(gtfObj);
		gtf.setSortAttributes(new String[] { "type", "region", "splicing",
				"chimp", "macaque", "ss_length" });
		gtf.addFileSuffix("_conserv_profile");
		try {
			gtf.write();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println((System.currentTimeMillis() - t0) / 1000 + " sec");
	
	}

	static double[][] calcConservation(boolean encode, GTFObject[] gtfObj, PrintStream p) {
		
		Vector v= new Vector(gtfObj.length);
		for (int i = 0; i < gtfObj.length; i++) {
			Object o= ENCODE.convertToEncodeCoord(gtfObj[i]);
			if (o!= null)
				v.add(o);
		}
		//System.out.println(gtfObj.length+" -> "+v.size());
		gtfObj= (GTFObject[]) gphase.tools.Arrays.toField(v);
		
		// "outputSS_encode_coding.encode_coord.gtf"
		// "test.gtf"
		long t0 = System.currentTimeMillis();
		// "outputSS_encode_coding_NEW.encode_coord.gtf"
		String[] species = { "human", "chimp", "macaque", "baboon", "galago",
				"marmoset" };
		v = new Vector();
		for (int i = 0; i < gtfObj.length; i++) {
			int j;
			for (j = 0; j < v.size(); j++)
				if (gtfObj[i].getSeqname().equals(v.elementAt(j)))
					break;
			if (j == v.size())
				v.add(gtfObj[i].getSeqname());
		}
		String[] encReg = (String[]) gphase.tools.Arrays.toField(v);
	
		// iterate regions
		Vector[][] result= new Vector[2][];	// 5UTR, CDS
		for (int i = 0; i < result.length; i++) {
			result[i]= new Vector[5];	// AD, AA, CD, CA
			for (int j = 0; j < result[i].length; j++) 
				result[i][j]= new Vector();
		}
		
		int[][] ssNr= new int[2][];		// count nrSSs
		for (int i = 0; i < ssNr.length; i++) {
			ssNr[i]= new int[4];
			for (int j = 0; j < ssNr.length; j++) {
				ssNr[i][j]= 0;
			}
		}
		Vector[][] ssColProfV= new Vector[2][];
		for (int i = 0; i < ssColProfV.length; i++) {
			ssColProfV[i]= new Vector[5];
			for (int j = 0; j < ssColProfV[i].length; j++) {
				ssColProfV[i][j]= new Vector();
			}
		}
		for (int i = 0; i < encReg.length; i++) {
			v = new Vector();
			for (int j = 0; j < gtfObj.length; j++)
				if (gtfObj[j].getSeqname().equalsIgnoreCase(encReg[i]))
					v.add(gtfObj[j]);
			GTFObject[] gtfEnc = (GTFObject[]) gphase.tools.Arrays.toField(v);
	
			SpliceSiteConservationComparator sp = new SpliceSiteConservationComparator(
					encReg[i], species);
			Graph g= EncodeWrapper.assemble(encode, gtfEnc);
			Vector[][] r= sp.calcConservation(g);
			for (int j = 0; j < r.length; j++) 
				for (int k = 0; k < r[j].length; k++) 
					result[j][k].addAll(r[j][k]);
			int[][][] ss= sp.getSsNr();	// join ss-Nr
			for (int j = 0; j < ss.length; j++) {	// utr/ cds
				for (int k = 0; k < ss[j].length; k++) {	// genes
					for (int m = 0; m < ss[j][k].length; m++) {	// ad, aa, cd, ca
						ssNr[j][m]+= ss[j][k][m];
					}
				}
			}
			int[][][][][] colProf= sp.getSsColNr();
			for (int j = 0; j < colProf.length; j++) {	// utr/ cds
				for (int k = 0; k < colProf[j].length; k++) {
					for (int m = 0; m < colProf[j][k].length; m++) {
						for (int n = 0; colProf[j][k][m]!= null&& n < colProf[j][k][m].length; n++) {
							ssColProfV[j][m].add(colProf[j][k][m][n]);
						}
					}
				}
			}
		}
	
		// calc avg
		double[] resUTR5= new double[result[0].length];
		for (int i = 0; i < result[0].length; i++) {
			double sum= 0d;
			for (int j = 0; j < result[0][i].size(); j++) 
				sum+= ((Double) result[0][i].elementAt(j)).doubleValue();
			resUTR5[i]= sum/ result[0][i].size(); 
		}
		double[] resCDS= new double[result[1].length];
		for (int i = 0; i < resCDS.length; i++) {
			double sum= 0d;
			for (int j = 0; j < result[1][i].size(); j++) 
				sum+= ((Double) result[1][i].elementAt(j)).doubleValue();
			resCDS[i]= sum/ result[1][i].size(); 
		}
		
		//System.out.println((System.currentTimeMillis() - t0) / 1000 + " sec");
		double[][] res= new double[][] {resUTR5, resCDS};
		gphase.tools.Arrays.output(res, p);
		p.println("N:");
		for (int i = 0; i < ssNr.length; i++) {
			for (int j = 0; j < ssNr[i].length; j++) 
				p.print(ssNr[i][j]+"\t");
			p.println();
		}
		String[] reg= new String[] {"5UTR", "3UTR"};
		String[] str= new String[] {"AD", "AA", "CD", "CA", "All"};
		for (int i = 0; i < ssColProfV.length; i++) {
			p.println(" == "+reg[i]+" == ");
			for (int j = 0; j < ssColProfV[i].length; j++) {
				p.println(str[j]);
				for (int k = 0; k < ssColProfV[i][j].size(); k++) {
					int[] prof= (int[]) ssColProfV[i][j].elementAt(k);
					if (prof!= null) {
						for (int m = 0; m < prof.length; m++) {
							p.print(prof[m]+" ");
						}
						p.println();
					}
				}
				p.println();
			}
			p.println("----------------------------------------------");
		}
		
		return res;
	}

	static double checkBackground() {
		
		String[] names = { "human", "chimp", "macaque", "baboon", "galago",
				"marmoset" };
		
		String[] encReg= ENCODE.getEncodeRegionNames();
		double avg= 0d;
		int cols= 0;	// nb of columns
		for (int i = 0; i < encReg.length; i++) {
			// extract sequences (always convert to forward for SS visual control)
			String[] seqBlock = new String[names.length];
			for (int j = 0; j < seqBlock.length; j++) 
				seqBlock[j] = getSequence(encReg[i], names[j]);
			
				// calculate conservation
			int x = 0;
			for (x = 0; x < names.length; x++)
				if (names[x].equalsIgnoreCase("human"))
					break;
			for (int k = 0; k < seqBlock[x].length(); k++) {
				if (isQuestionable(seqBlock[x].charAt(k)))
					continue;
				char c= Character.toUpperCase(seqBlock[x].charAt(k));
				int base= 0;
				int csvd= 0;
				for (int j = 0; j < seqBlock.length; j++) {
					if (j== x)
						continue;
					if (!isQuestionable(seqBlock[j].charAt(k))) {
						++base;
						if (Character.toUpperCase(seqBlock[j].charAt(k))== c)
							++csvd;
					}
				}
				if (base> 0) { 
					avg+= csvd/ ((double) base);
					++cols;
				}
			}
		}
		avg/= cols;
		System.out.println(avg);
		return avg;
	}

	static void calcConservationChk(boolean encode, GTFObject[] gtfObj, PrintStream p) {
		
		Vector v= new Vector(gtfObj.length);
		for (int i = 0; i < gtfObj.length; i++) {
			Object o= ENCODE.convertToEncodeCoord(gtfObj[i]);
			if (o!= null)
				v.add(o);
		}
		//System.out.println(gtfObj.length+" -> "+v.size());
		gtfObj= (GTFObject[]) gphase.tools.Arrays.toField(v);
		
		// "outputSS_encode_coding.encode_coord.gtf"
		// "test.gtf"
		long t0 = System.currentTimeMillis();
		// "outputSS_encode_coding_NEW.encode_coord.gtf"
		String[] species = { "human", "chimp", "macaque", "baboon", "galago",
				"marmoset" };
		v = new Vector();
		for (int i = 0; i < gtfObj.length; i++) {
			int j;
			for (j = 0; j < v.size(); j++)
				if (gtfObj[i].getSeqname().equals(v.elementAt(j)))
					break;
			if (j == v.size())
				v.add(gtfObj[i].getSeqname());
		}
		String[] encReg = (String[]) gphase.tools.Arrays.toField(v);
	
		int[][] result= new int[2][];	// UTR, CDS
		for (int j = 0; j < result.length; j++) {
			result[j]= new int[2];		// AD, AA
			for (int i = 0; i < result[j].length; i++) 
				result[j][i]= 0;
		}
		for (int i = 0; i < encReg.length; i++) {
			v = new Vector();
			for (int j = 0; j < gtfObj.length; j++)
				if (gtfObj[j].getSeqname().equalsIgnoreCase(encReg[i]))
					v.add(gtfObj[j]);
			GTFObject[] gtfEnc = (GTFObject[]) gphase.tools.Arrays.toField(v);
	
			SpliceSiteConservationComparator sp = new SpliceSiteConservationComparator(
					encReg[i], species);
			Graph g= EncodeWrapper.assemble(encode, gtfEnc);
			int[][] r= sp.countEvents(g);
			for (int j = 0; j < r.length; j++) {
				for (int k = 0; k < r[j].length; k++) 
					result[j][k]+= r[j][k];
			}
		}
	
		System.out.println(result[0][0]+"\t"+result[0][1]);
		System.out.println(result[1][0]+"\t"+result[1][1]);
	}

	static void conservationAnalysis() {
		// "outputSS_encode_coding.encode_coord.gtf"
		// "test.gtf"
		long t0 = System.currentTimeMillis();
		// "outputSS_encode_coding_NEW.encode_coord.gtf"
		String inputFile = "SSout_new.encode_coord.gtf";
		String[] species = { "human", "chimp", "macaque", "baboon", "galago",
				"marmoset" };
		int delta = 7;
		GTFWrapper gtf = new GTFWrapper(new File(inputFile).getAbsolutePath());
	
		// collect enc regions
		Vector v = new Vector();
		try {
			gtf.read();
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		GTFObject[] gtfObj = gtf.getGtfObj();
		for (int i = 0; i < gtfObj.length; i++) {
			int j;
			for (j = 0; j < v.size(); j++)
				if (gtfObj[i].getSeqname().equals(v.elementAt(j)))
					break;
			if (j == v.size())
				v.add(gtfObj[i].getSeqname());
		}
		String[] encReg = (String[]) gphase.tools.Arrays.toField(v);
	
		// iterate regions
		Vector resultV = new Vector();
		for (int i = 0; i < encReg.length; i++) {
			v = new Vector();
			for (int j = 0; j < gtfObj.length; j++)
				if (gtfObj[j].getSeqname().equalsIgnoreCase(encReg[i]))
					v.add(gtfObj[j]);
			GTFObject[] gtfEnc = (GTFObject[]) gphase.tools.Arrays.toField(v);
	
			SpliceSiteConservationComparator sp = new SpliceSiteConservationComparator(
					encReg[i], species);
			GTFObject[] tmpObj = sp.calcConservation(gtfEnc, delta);
			for (int j = 0; j < tmpObj.length; j++)
				resultV.add(tmpObj[j]);
		}
	
		gtfObj = (GTFObject[]) gphase.tools.Arrays.toField(resultV);
		gtf.setGtfObj(gtfObj);
		gtf.setSortAttributes(new String[] { "type", "region", "splicing",
				"chimp", "macaque", "ss_length" });
		gtf.addFileSuffix("_conserv_profile");
		try {
			gtf.write();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println((System.currentTimeMillis() - t0) / 1000 + " sec");
	
	}

	/**
	 * 
	 * @param relPos
	 *            0-based
	 * @param seq
	 * @return
	 */
	public int getAbsPos(int relPos, String seq) {
		int ctr = 0;
		int i;
		for (i = 0; i < seq.length(); i++) {
			if (ctr == relPos)
				break;
			if (seq.charAt(i) == '-')
				continue;
			ctr++;
		}

		if (i == seq.length())
			return -1;
		return i;
	}

	public SpliceSiteConservationComparator(String newEncRegion,
			String[] species) {
		this.encRegion = newEncRegion;
		this.names = species;
		// readIn(newEncRegion, species);
	}

	void readIn(String encRegion, String[] species) {

		// ali= new String[species.length];
		// names= new String[species.length];
		// for (int i = 0; i < species.length; i++) {
		// String absPath= new
		// File(ALI_DIR+encRegion).getAbsolutePath()+".mfa"+"_"+species[i];
		// FASTAWrapper fasta= new FASTAWrapper(absPath);
		// try {
		// fasta.read();
		// } catch (Exception e) {
		// e.printStackTrace();
		// }
		// ali[i]= fasta.getSequences()[0];
		// names[i]= fasta.getSeqNames()[0];
		// }
	}

	String getSubstring(String spec, int pos1, int pos2) {
		try {
			RandomAccessFile raf = new RandomAccessFile(new File(ALI_DIR
					+ encRegion + ".mfa.proj_" + spec), "r");
			String bline = "";
			while (!bline.startsWith(">"))
				bline = raf.readLine().trim();
	
			int lineLength = 60; // dirty, but only 5min time !
			int offset = bline.length() + 1 // 1st line
					+ pos1 + (pos1 / lineLength);
			raf.seek(offset);
	
			String line = raf.readLine();
			int diff = pos2 - pos1;
			if (line.length() > diff)
				return line.substring(0, diff);
			while (line.length()< diff)
				line += raf.readLine();
			return line.substring(0, diff);
		} catch (Exception e) {
			e.printStackTrace();
		}
	
		return null;
	}

	static String getSequence(String encRegion, String spec) {
		try {
			File f= new File(ALI_DIR
					+ encRegion + ".mfa.proj_" + spec);
			BufferedReader raf = new BufferedReader(new FileReader(f));
			String bline = "";
			while (!bline.startsWith(">"))
				bline = raf.readLine().trim();

			StringBuffer result= new StringBuffer((int) f.length());	// hoping its not bigger than int
			String line = raf.readLine();
			while (line!= null&& line.length()> 0) {
				result.append(line);
				line= raf.readLine();
			}
			return result.toString();
		} catch (Exception e) {
			e.printStackTrace();
		}

		return null;
	}

	// pos1 incl, pos2 excl
	public static String getSubstring(String encRegName, String specName,
			int pos1, int pos2) {
		try {
			RandomAccessFile raf = new RandomAccessFile(new File(ALI_DIR
					+ encRegName + ".mfa.proj_" + specName), "r");
			String bline = "";
			while (!bline.startsWith(">"))
				bline = raf.readLine().trim();

			int lineLength = 60; // dirty, but only 5min time !
			int offset = bline.length() + 1 // 1st line
					+ pos1 + (pos1 / lineLength);
			raf.seek(offset);

			String line = raf.readLine().trim();
			int diff = pos2 - pos1;
			while (line.length() <= diff)
				line += raf.readLine();
			return line.substring(0, diff);
		} catch (Exception e) {
			e.printStackTrace();
		}

		return null;
	}

	public GTFObject[] calcConservation(GTFObject[] gtfObj, int newDelta) {
	
		int x = 0;
		for (x = 0; x < names.length; x++)
			if (names[x].equalsIgnoreCase("human"))
				break;
		for (int i = 0; i < gtfObj.length; i++) {
			PrintStream pr = null;
			try {
				pr = new PrintStream(new BufferedOutputStream(
						new FileOutputStream(gtfObj[i].getSeqname(), true)));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
	
			// don: 3 + 2 + 4 = 9
			// acc: 6 + 2 + 1 = 9
			// don: 12 + 2 + 13 = 27
			// acc: 15 + 2 + 10 = 27
			int relPos = -1;
			int delta = 37;
			int offset = 0;
			int donShift = -12;
			int accShift = 10;
			int symShift = 10;
			if (gtfObj[i].isStrand()) {
				if (gtfObj[i].getAttribute("type").startsWith("don")) {
					relPos = gtfObj[i].getStart() + offset + donShift;
				} else if (gtfObj[i].getAttribute("type").startsWith("acc")) {
					relPos = gtfObj[i].getEnd() + 1 - offset + accShift;
					delta = -delta;
				} else {
					relPos = gtfObj[i].getStart() - symShift;
					delta = 2 * symShift
							+ (gtfObj[i].getEnd() - gtfObj[i].getStart());
				}
			} else {
				if (gtfObj[i].getAttribute("type").startsWith("don")) {
					relPos = gtfObj[i].getEnd() + 1 - offset - donShift;
					delta = -delta;
				} else if (gtfObj[i].getAttribute("type").startsWith("acc")) {
					relPos = gtfObj[i].getStart() + offset - accShift;
				} else {
					relPos = gtfObj[i].getEnd() + symShift;
					delta = 2 * symShift
							+ (gtfObj[i].getEnd() - gtfObj[i].getStart());
				}
			}
	
			if (delta < 0) { // swap start/end
				relPos = relPos + delta;
				delta = -delta;
			}
			--relPos; // convert to 0-based
	
			// int absPos= getAbsPos(relPos, ali[x]);
			// if (absPos< 0) {
			// System.err.println(gtfObj[i].getStart()+ " outside of "+
			// gtfObj[i].getSeqname());
			// pr.println(gtfObj[i].getStart()+ " outside of "+
			// gtfObj[i].getSeqname());
			// continue;
			// }
	
			// int maxDelta= 0; // get max absolute Delta
			// for (int j = 0; j < ali.length; j++) {
			// String subseq= ali[j].substring(absPos);
			// int absDelta= getAbsPos(newDelta, subseq);
			// maxDelta= Math.max(maxDelta, absDelta);
			// }
			int maxDelta = delta; // get max absolute Delta
	
			// get subseqs and compare
			pr.println(gtfObj[i].getAttribute("type") + " "
					+ gtfObj[i].getStrand() + " "
					+ gtfObj[i].getAttribute("region") + " "
					+ gtfObj[i].getAttribute("splicing") + " "
					+ gtfObj[i].getSeqname() + " " + gtfObj[i].getStart() + " "
					+ gtfObj[i].getEnd());
			String[] seqBlock = new String[names.length];
			for (int j = 0; j < seqBlock.length; j++)
				seqBlock[j] = getSubstring(names[j], relPos, relPos + maxDelta);
	
			// for (int j = 0; j < seqBlock.length; j++) {
			// int mis= 0;
			// if (j!= x) {
			// for (int k = 0; k < seqBlock[x].length(); k++)
			// if (seqBlock[x].charAt(k)!= seqBlock[j].charAt(k))
			// mis++;
			// gtfObj[i].addAttribute(names[j], Integer.toString(mis));
			// gtfObj[i].addAttribute("ss_length", Integer.toString(delta));
			// }
			// pr.println(names[j]+"\t"+seqBlock[j]+"\t"+delta+"\t"+mis);
			// }
			// pr.flush(); pr.close();
	
			if (!gtfObj[i].isStrand())
				seqBlock[x] = gphase.tools.Arrays.reverseComplement(seqBlock[x]);
			boolean yes = false;
			String chk;
			if (gtfObj[i].getAttribute("type").startsWith("don"))
				chk = seqBlock[x].substring(-donShift, -donShift + 2);
			else
				chk = seqBlock[x].substring(delta - accShift - 2, delta
						- accShift);
			if ((gtfObj[i].getAttribute("type").startsWith("don") && seqBlock[x]
					.substring(-donShift, -donShift + 2).equalsIgnoreCase("GT"))
					|| ((gtfObj[i].getAttribute("type").startsWith("acc") && seqBlock[x]
							.substring(delta - accShift - 2, delta - accShift)
							.equalsIgnoreCase("AG"))))
				yes = true;
			if (!yes)
				System.currentTimeMillis();
			for (int j = 0; yes && j < seqBlock.length; j++) {
				int[] cons = new int[seqBlock[j].length()];
				for (int k = 0; k < cons.length; k++)
					cons[k] = 1;
				if (j != x) {
					if (!gtfObj[i].isStrand()) { // watch out! i is in
														// human!!
						// Arrays.reverse(cons);
						seqBlock[j] = gphase.tools.Arrays.reverseComplement(seqBlock[j]);
					}
					for (int k = 0; k < seqBlock[x].length(); k++)
						if (Character.toUpperCase(seqBlock[x].charAt(k)) == Character
								.toUpperCase(seqBlock[j].charAt(k)))
							cons[k] = 1;
						else
							cons[k] = 0;
					gtfObj[i].addAttribute((names[j] + "_cons_profile"), Arrays
							.toString(cons));
					gtfObj[i]
							.addAttribute("ss_length", Integer.toString(delta));
				}
				pr.println(names[j] + "\t" + seqBlock[j] + "\t" + delta + "\t"
						+ Arrays.toString(cons));
			}
			pr.flush();
			pr.close();
		}
	
		return gtfObj;
	}

	public Vector[][] calcConservation(Graph g) {
		
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ssNr.length; i++) {	// nr SSs per gene
			ssNr[i]= new int[ge.length][];
			ssColNr[i]= new int[ge.length][][][];
			for (int j = 0; j < ssNr[i].length; j++) {
				ssNr[i][j]= new int[4];
				ssColNr[i][j]= new int[5][][];
				for (int k = 0; k < ssNr[i][j].length; k++) 
					ssNr[i][j][k]= 0;
			}
		}
	
		Vector[] v5= new Vector[5], vCDS= new Vector[5];
		for (int i = 0; i < v5.length; i++) { 
			v5[i]= new Vector();
			vCDS[i]= new Vector();
		}
		
		for (int i = 0; i < ge.length; i++) {
			
			ASVariation[] vars= ge[i].getASVariations(ASMultiVariation.FILTER_CODING_REDUNDANT);
			
			SpliceSite[] s0= ge[i].getSpliceSites();
			if (s0== null|| s0.length< 1) {
				ssColNr[0][i]= new int[0][][];
				ssColNr[1][i]= new int[0][][];
				continue;
			}
			DirectedRegion utr5= ge[i].getReal5UTR();
			if (utr5!= null) {
				ASVariation[] ev5UTR= (ASVariation[]) gphase.tools.Arrays.filter(vars, "is5UTRFunc");
				SpliceSite[] ss= (SpliceSite[]) DirectedRegion.contained(utr5, s0);
				
				double[] r= calcConservationAD_AA_CD_CA_All(utr5, ev5UTR, ss, ssNr[0][i], ssColNr[0][i]);
				for (int j = 0; j < r.length; j++) 
					if (r[j]>= 0)
						v5[j].add(new Double(r[j]));
			} 
			
			DirectedRegion cds= ge[i].getMaxCDS();		// getRealCDS()
			if (cds!= null) {
				ASVariation[] evCDS= (ASVariation[]) gphase.tools.Arrays.filter(vars, "isCodingFunc");
				SpliceSite[] ss= (SpliceSite[]) DirectedRegion.contained(cds, s0);
				double[] r= calcConservationAD_AA_CD_CA_All(cds, evCDS, ss, ssNr[1][i], ssColNr[1][i]);
				for (int j = 0; j < r.length; j++) 
					if (r[j]>= 0)
						vCDS[j].add(new Double(r[j]));
			} 
		}
		
		return new Vector[][] {v5, vCDS};
	}

	private int[][] countEvents(Graph g) {
		
		Gene[] ge= g.getGenes();
		int[][] result= new int[2][];	// UTR, CDS
		for (int j = 0; j < result.length; j++) {
			result[j]= new int[2];		// AD, AA
			for (int i = 0; i < result[j].length; i++) 
				result[j][i]= 0;
		}
		
		for (int i = 0; i < ge.length; i++) {
			
			ASVariation[] vars= ge[i].getASVariations(ASMultiVariation.FILTER_CODING_REDUNDANT);
			
			SpliceSite[] s0= ge[i].getSpliceSites();
			DirectedRegion utr5= ge[i].getReal5UTR();
			ASVariation[] ev5UTR= (ASVariation[]) gphase.tools.Arrays.filter(vars, "is5UTRFunc");
			if (utr5!= null) {				
				result[0]= countAD_AA_CD_CA_All(utr5, ev5UTR, result[0]);
				SpliceSite[] ss= (SpliceSite[]) DirectedRegion.contained(utr5, s0);
			} 
			
			DirectedRegion cds= ge[i].getMaxCDS();		// getRealCDS()
			ASVariation[] evCDS= (ASVariation[]) gphase.tools.Arrays.filter(vars, "isCodingFunc");
			if (cds!= null) {
				result[1]= countAD_AA_CD_CA_All(cds, evCDS, result[1]);
				SpliceSite[] ss= (SpliceSite[]) DirectedRegion.contained(cds, s0);
			} 
		}
		
		return result;
	}
	
	private int[] countAD_AA(ASVariation[] evts, int[] r) {
		int ctrAA= 0;
		int ctrAD= 0;
		for (int j = 0; evts!= null&& j < evts.length; j++) {
			if (evts[j].toString().equals(ASVariation.ID_PURE_AA))
					++ctrAA;
			if (evts[j].toString().equals(ASVariation.ID_PURE_AD))
				++ctrAD;
		}
		return new int[] {ctrAD, ctrAA};
	}
	
	private double[] calcConservationAD_AA_CD_CA_All(DirectedRegion reg, ASVariation[] vars, SpliceSite[] ss, int[] ssCnt, int[][][] ssColCnt) {
		int don5F= 2;	// 5'flank of donor
		int don3F= 4;	// 6
		int acc5F= 13;	// 15
		int acc3F= 2;
		
		double csvAD= -1d, csvAA= -1d, csvCD= -1d, csvCA= -1d, csvAll= -1d;
		
		ASVariation[] evAD= ASAnalyzer.getVariation("(1^ // 2^)", vars);
		if (evAD!= null) {
			for (int i = 0; i < evAD.length; i++) {
				//SpliceSite[] ad= extractSS(evAD);
				SpliceSite[] ad= evAD[i].getSpliceUniverse();
				if (ad[0].getPos()< ad[1].getPos())
					ad= new SpliceSite[] {ad[0]};
				else
					ad= new SpliceSite[] {ad[1]};
				ssCnt[0]= ad.length;
				ssColCnt[0]= new int[ad.length][];
				if (evAD[i].isNotProteinCoding())
					csvAD= calcConserv(ad, don5F, don3F, reg, ssColCnt[0], evAD[i]);
				else
					csvAD= calcConserv(ad, don5F, don3F, reg, ssColCnt[0], null);
			}
		} else {
			ssColCnt[0]= new int[1][];
		}
		
			// for ad_prox/dist analysis
		if (evAD!= null) {
//			SpliceSite[] ad= extractSS(evAD);
			for (int i = 0; i < evAD.length; i++) {
				//SpliceSite[] ad= extractSS(evAD);
				SpliceSite[] ad= evAD[i].getSpliceUniverse();
				if (ad[0].getPos()> ad[1].getPos())
					ad= new SpliceSite[] {ad[0]};
				else 
					ad= new SpliceSite[] {ad[1]};
				ssCnt[1]= ad.length;
				ssColCnt[1]= new int[ad.length][];
				if (evAD[i].isNotProteinCoding())
					csvAA= calcConserv(ad, don5F, don3F, reg, ssColCnt[1], evAD[i]);
				else
					csvAA= calcConserv(ad, don5F, don3F, reg, ssColCnt[1], null);
			}
		} else {
			ssColCnt[1]= new int[1][];
		}
		
//		ASVariation[] evAA= ASAnalyzer.getVariation("(1= // 2=)", vars);
//		if (evAA!= null) {
//			SpliceSite[] aa= extractSS(evAA);
//			ssCnt[1]= aa.length;
//			ssColCnt[1]= new int[aa.length][];
//			csvAA= calcConserv(aa, acc5F, acc3F, reg, ssColCnt[1], null);
//		} else {
//			ssColCnt[1]= new int[1][];
//		}
	
		SpliceSite[] cs= (SpliceSite[]) gphase.tools.Arrays.filter(ss, "isConstitutive");
		SpliceSite[] cd= (SpliceSite[]) gphase.tools.Arrays.filter(cs, "isDonor");
		if (cd!= null) {
			ssCnt[2]= cd.length;
			ssColCnt[2]= new int[cd.length][];
			csvCD= calcConserv(cd, don5F, don3F, reg, ssColCnt[2], null);
		} else {
			ssColCnt[2]= new int[1][];
		}
		
		SpliceSite[] ca= (SpliceSite[]) gphase.tools.Arrays.filter(cs, "isAcceptor");
		if (ca!= null) {
			ssCnt[3]= ca.length;
			ssColCnt[3]= new int[ca.length][];
			csvCA= calcConserv(ca, acc5F, acc3F, reg, ssColCnt[3], null);
		} else {
			ssColCnt[3]= new int[1][];
		}
		
		ssColCnt[4]= new int[1][];
		if (reg.getStart()< reg.getEnd()) { 	// only valid regions
			ssColCnt[4][0]= new int[reg.getEnd()- reg.getStart()+ 1];
			csvAll= calcConserv(reg, ssColCnt[4][0]);
		} else {
			ssColCnt[4]= new int[0][];
		}
		
		return new double[] {csvAD, csvAA, csvCD, csvCA, csvAll};
	}

	private int[] countAD_AA_CD_CA_All(DirectedRegion reg, ASVariation[] vars, int[] ssCnt) {
		int don5F= 2;	// 5'flank of donor
		int don3F= 6;
		int acc5F= 15;
		int acc3F= 2;
		
		int csvAD= 0, csvAA= 0, csvCD= 0, csvCA= 0, csvAll= 0;
		
		ASVariation[] evAD= ASAnalyzer.getVariation("(1^ // 2^)", vars);
		if (evAD!= null) {
			SpliceSite[] ad= extractSS(evAD);
			ssCnt[0]+= count(ad, don5F, don3F, reg);
		}
		
		ASVariation[] evAA= ASAnalyzer.getVariation("(1= // 2=)", vars);
		if (evAA!= null) {
			SpliceSite[] aa= extractSS(evAA);
			ssCnt[1]+= count(aa, acc5F, acc3F, reg);
		} 

		return ssCnt;
	}
	
	double calcConserv(SpliceSite[] ss, int f5, int f3, DirectedRegion reg, int[][] ssColCnt, ASVariation var) {
	
			// SSs, extract regions
		double result= 0d;
		if (ss!= null) {
			for (int i = 0; i < ss.length; i++) {
				
//				if (!reg.contains(ss[i]))
//					continue;
				
				int p1, p2;
				if (ss[i].isDonor()) {
					p1= ss[i].getPos()+ 1- f5;	// exon border -> ss dinucleotide
					p2= ss[i].getPos()+ 2+ f3;
				} else {
					p1= ss[i].getPos()- 2- f5;
					p2= ss[i].getPos()- 1+ f3;
				}
				p1= Math.max(reg.get5PrimeEdge(), p1);
				p2= Math.min(reg.get3PrimeEdge(), p2);
				
				DirectedRegion reg2= new DirectedRegion(p1, p2, reg.getStrand());
				
				ssColCnt[i]= new int[Math.max(p1,p2)- Math.min(p1,p2)+ 1];	
				for (int j = 0; j < ssColCnt[i].length; j++) 
					ssColCnt[i][j]= 0;
				double delta= calcConserv(reg2, ssColCnt[i]);
				result+= delta;
				if (var!= null&& ss[i].isDonor()&& !ss[i].isConstitutive()) {
					int dist= 0;
					if (var.getTranscript1().containsSS(ss[i]))
						dist= ss[i].getPos()- var.getTranscript1().get5PrimeEdge();
					else
						dist= ss[i].getPos()- var.getTranscript2().get5PrimeEdge();
					System.out.println(dist+"\t"+delta);
				}
			}
			result= (result/ ss.length);
		}
	
		return result;
	}

	int count(SpliceSite[] ss, int f5, int f3, DirectedRegion reg) {

			// SSs, extract regions
		int success= 0;
		if (ss!= null) {
			for (int i = 0; i < ss.length; i++) {
				int p1, p2;
				if (ss[i].isDonor()) {
					p1= ss[i].getPos()+ 1- f5;	// exon border -> ss dinucleotide
					p2= ss[i].getPos()+ 2+ f3;
				} else {
					p1= ss[i].getPos()- 2- f5;
					p2= ss[i].getPos()- 1+ f3;
				}
//				p1= Math.max(reg.get5PrimeEdge(), p1);
//				p2= Math.min(reg.get3PrimeEdge(), p2);
				
				DirectedRegion reg2= new DirectedRegion(p1, p2, reg.getStrand());
				reg2.setChromosome(reg.getChromosome());
				
				if (reg.contains(reg2))
					++success;
				else {
					reg2= ENCODE.convertToChromosomalCoord(reg2);
					System.currentTimeMillis();
				}
			}
		}

		return success;
	}
	
	double calcConserv(DirectedRegion reg, int[] colProfile) {
	
			// prepare coordinates
		int p1= reg.getStart();
		int p2= reg.getEnd();
		if (Math.abs(p1)> Math.abs(p2)) {	// change start/end on reverse
			int h= p1;
			p1= p2;
			p2= h;
		}
		p1= Math.abs(p1);
		p2= Math.abs(p2);
		
			// extract sequences (always convert to forward for SS visual control)
		String[] seqBlock = new String[names.length];
		for (int j = 0; j < seqBlock.length; j++) {
			seqBlock[j] = getSubstring(names[j], p1, p2+ 1);	// end excluded
			if (!reg.isLeadingStrand())
				seqBlock[j] = gphase.tools.Arrays.reverseComplement(seqBlock[j]);
		}
		
			// calculate conservation
		int x = 0;
		for (x = 0; x < names.length; x++)
			if (names[x].equalsIgnoreCase("human"))
				break;
		double[] result= new double[seqBlock[x].length()];
		for (int i = 0; i < seqBlock[x].length(); i++) {
			if (isQuestionable(seqBlock[x].charAt(i)))
				continue;
			char c= Character.toUpperCase(seqBlock[x].charAt(i));
			int base= 0;
			int csvd= 0;
			for (int j = 0; j < seqBlock.length; j++) {
				if (j== x)
					continue;
				if (!isQuestionable(seqBlock[j].charAt(i))) {
					++base;
					if (Character.toUpperCase(seqBlock[j].charAt(i))== c)
						++csvd;
				}
			}
			if (base> 0) {
				result[i]= csvd/ ((double) base);
				colProfile[i]= base;
			}
			else
				result[i]= 0d;
		}
		
			// calculate average
		double cnsv= 0d;
		int col= 0;
		for (int j = 0; j < result.length; j++) {
			if (result[j]!= 0d) {
				cnsv+= result[j];
				++col;
			}
		}
	
		return (cnsv/ col);
	}

	private static boolean isQuestionable(char c) {
		if (Character.toUpperCase(c)== 'N'
			|| !Character.isLetter(c))
			//|| Character.isLowerCase(c))
			return true;
		return false;
	}
	
	public double[][] getPWDistances(SpliceSite[] ss) {

		double[][] cons = new double[names.length][];
		for (int j = 0; j < cons.length; j++)
			cons[j] = new double[names.length];

		int x = 0;
		for (x = 0; x < names.length; x++)
			if (names[x].equalsIgnoreCase("human"))
				break;
		int a= 0;
		for (int i = 0; i < ss.length; i++) {
			
			DefaultRegion reg= ENCODE.getEncodeRegion(Math.abs(ss[i].getPos()));
			if (reg== null) {
				System.err.println("outsite encode "+ ss[i].getGene().getChromosome()+ " "+ Math.abs(ss[i].getPos()));
				continue;
			} 
			//if (encRegion== null)
				encRegion= reg.getID();
//>			else if (!reg.getID().equalsIgnoreCase(encRegion))
//				continue; 
			
			PrintStream pr = null;
			try {
				pr = new PrintStream(new BufferedOutputStream(
						new FileOutputStream(encRegion, true)));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}

			// don: 12 + 2 + 13 = 27
			// acc: 15 + 2 + 10 = 27
			int relPos = -1;
			int delta = 27;
			int offset = 0;
			int donShift = -12;	//-12;
			int accShift = 10;	//10;
			if (ss[i].getGene().isForward()) {
				if (ss[i].isDonor()) {
					relPos = ENCODE.convertToEncodeCoord(ss[i]).getStart() + offset + donShift;
				} else {
					relPos =  ENCODE.convertToEncodeCoord(ss[i]).getEnd() + 1 - offset + accShift;
					delta = -delta;
				}
			} else {
				if (ss[i].isDonor()) {
					relPos = ENCODE.convertToEncodeCoord(ss[i]).getEnd() + 1 - offset - donShift;
					delta = -delta;
				} else {
					relPos = ENCODE.convertToEncodeCoord(ss[i]).getStart() + offset - accShift+ 2;
				}
			}

			if (delta < 0) { // swap start/end
				relPos = relPos + delta;
				delta = -delta;
			}
			--relPos; // convert to 0-based

			// get subseqs and compare
//			pr.println(encRegion + " " + ss[i].getGene().getStrand() + " "
//					+ ss[i].getGene().getChromosome() + " " + Math.abs(ss[i].getPos())
//					+ (ss[i].isDonor() ? "don" : "acc"));
			String[] seqBlock = new String[names.length];
			for (int j = 0; j < seqBlock.length; j++) {
				seqBlock[j] = getSubstring(names[j], relPos, relPos + delta);
				if (!ss[i].getGene().isForward())
					seqBlock[j] = gphase.tools.Arrays.reverseComplement(seqBlock[j]);
			}

			boolean yes = false;
			String chk;
			int chkPos= -1;
			if (ss[i].isDonor()) {
				if (ss[i].getGene().isForward()) 
					chkPos= -donShift;
				else
					chkPos= delta+donShift-3;
			} else {
				if (ss[i].getGene().isForward()) 
					chkPos= delta - accShift - 3;
				else
					chkPos= delta - accShift - 2;
			}
			chk = seqBlock[x].substring(chkPos, chkPos + 2);
			
			if ((ss[i].isDonor() && seqBlock[x].substring(chkPos,
					chkPos + 2).equalsIgnoreCase("GT"))
					|| (!ss[i].isDonor() && seqBlock[x].substring(
							chkPos, chkPos+ 2)
							.equalsIgnoreCase("AG")))
				yes = true;
			
			for (int j = 0; j < seqBlock.length; j++) {
				if (j== x)
					continue;
				if ((ss[i].isDonor() && !seqBlock[j].substring(chkPos,
						chkPos + 2).equalsIgnoreCase("GT"))
						|| (!ss[i].isDonor() && !seqBlock[j].substring(
								chkPos, chkPos+ 2)
								.equalsIgnoreCase("AG"))) {
					a++;
					break;
				}
			}
			System.currentTimeMillis();
			for (int j = 0; yes && j < seqBlock.length; j++) {
				for (int k = j + 1; k < cons.length; k++) {
					int dd = 0;
					for (int m = 0; m < seqBlock[k].length(); m++)
						if (Character.toUpperCase(seqBlock[k].charAt(m)) != Character
								.toUpperCase(seqBlock[j].charAt(m)))
							++dd;
					cons[j][k] = (double) dd / (seqBlock[k].length() - 2);
					cons[k][j] = cons[j][k];
				}
			}
		}

		System.out.println(a+" / "+ss.length);
		return cons;
	}

	public int[][][] getSsNr() {
		return ssNr;
	}

	public int[][][][][] getSsColNr() {
		return ssColNr;
	}

}
