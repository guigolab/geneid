/*
 * Created on Oct 18, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Vector;



import qalign.algo.dialign.Constants;
import qalign.tools.MSFWrapper;
import repeat.data.AlignedTupel;
import repeat.data.Alignment;
import repeat.data.Repeat;

/**
 * 
 * 
 * @author micha
 */
public class PSComparator {
	
	String[] seqNamesTest= null;
	String[] seqNamesRef= null;
	double avgHitPercent= -1d;	
	
	/**
	 * alignRelation[seqA][seqB][RepA]= Rep# of main alignment partner of RepA in B
	 */
	int[][][] repAligned= null;
	
	/**
	 * repProfile[SeqA][SeqB][RepA]= repeats of SeqB aligned with RepA
	 */
	int[][][][] repProfile= null;
	/**
	 * refTupels[seqA][seqB][repA][repB][]
	 */
	AlignedTupel[][][][][] refTupels= null;
	/**
	 * refCoreBlocks[seqA][seqB][][]
	 */
	int[][][][] refCoreTupelCount= null;
	
	/**
	 * testTupels[seqA][seqB][]
	 */ 
	AlignedTupel[][][][][] testTupels= null;
	
	/**
	 * hits[seqA][seqB][repA][repB]
	 */
	int[][][][] hitCountAll= null;

	/**
	 * hits[seqA][seqB][repA][repB]
	 */
	int[][][][] hitCountCore= null;	
	int[][] seq2HitCountAll= null;
	
	int[][] seq2HitCountCore= null;
	int totalHitCountAll= -1;
	int totalHitCountCore= -1;
	
	String fileBase= null;
	static final String fileExt= null;
	BAliManager baliRef= null;
	Alignment aliTest= null;
	
	static BufferedWriter statsWriter= null;
	
	public static void start(Class baseClass, String ext) {
		
		try {									// create log writer
			statsWriter= new BufferedWriter(new FileWriter(ext+".tab"));
		} catch (IOException e) {
			; // :)
		}
		
		Class[] pars= {String.class};			// get redirection method
		Method m= null;
		Constructor c= null;
		try {
			c= baseClass.getConstructor(new Class[] {String.class, String.class});
			m= PSComparator.class.getMethod("writeStats",null);
		} catch (NoSuchMethodException e) {
			; // :)
		}
		
		BaliConqueror conqueror= new BaliConqueror(c,m,ext);	// iterate BaliBase
		conqueror.conquer();
		
		try {									// close writer
			statsWriter.flush();
			statsWriter.close();
		} catch (IOException e) {
			; // :)
		}
	}
	
	int getRepAlignmentPartner(AlignedTupel[][] tst, boolean tst2) {
		
		return 0;
	}
	
	public int[][][] getRepAligned(boolean triggerCore) {
		
		if (repAligned == null) {
			
				// init
			getRefTupels();
			getTestTupels();
			
				// determine alignment partner
			repAligned= new int[testTupels.length][][];
			for (int i= 0; i < testTupels.length; i++) {
				repAligned[i]= new int[testTupels[i].length][];
				for (int j= 0; j < repAligned[i].length; j++) {
					repAligned[i][j]= new int[testTupels[i][j].length];
					for (int k= 0; k < repAligned[i][j].length; k++) 
						repAligned[i][j][k]= getRepAlignmentPartner(testTupels[i][j][k], triggerCore);
				}
			}
			
		}

		return repAligned;
	}
	
	public void writeStats() throws IllegalArgumentException {

			// check
		if (aliTest== null|| baliRef== null)
			throw new IllegalArgumentException();
			
			// stat string
		//compareTestRef();
		String baliFile=
			fileBase.substring(
				fileBase.substring(0,fileBase.lastIndexOf(File.separator)).lastIndexOf(File.separator)+ 1,
				fileBase.length()
			);
		float accCore= ((float) getTotalHitCountCore())/((float) (getTotalHitCountCore()+ getXCountCore()));
		float accAll= ((float) getTotalHitCountAll())/((float) (getTotalHitCountAll()+ getXCountAll()));
		String str= 
			baliFile+ "\t"+
			getTotalHitCountCore()+ "\t"+getTotalHitCountAll()+"\t"
			+getXCountCore()+"\t"+getXCountAll()+"\t"
			+accCore+"\t"+accAll+"\n";
//			((float) getAverageIdentity(true))+"\t"+getAverageIdentity(false)+"\n";
//			getAvgHitPercent(true)+ " "+getAvgHitPercent(false)+"\t"+
//			getMinMaxHitPercent(true,true)+ " "+getMinMaxHitPercent(false,true)+"\t"+
//			getMinMaxHitPercent(true,false)+ " "+getMinMaxHitPercent(false,false)+"\n";
		
			// write to writer
		System.out.println(str);
		try {
			statsWriter.write(str);
			statsWriter.flush();
		} catch (IOException e) {
		}
	}
	
	public int getXCountCore() {

		if (xCountCore< 0) {
			countX();
		}
		return xCountCore;
	}
	
	public int getXCountAll() {
		
		if (xCountAll< 0) {
			countX();
		}
		return xCountAll;
	}
	
	int xCountCore= -1, xCountAll= -1;
	/**
	 * @deprecated xCountCore not working!
	 * @param core
	 * @return
	 */
	void countX() {
		
		xCountCore= 0;
		xCountAll= 0;
		initRepeatArray();
		
			// derive gapfree repeats (only for core)
		baliRef.getAlis();
		int ll= 0;		
		for (int i= 0; i < baliRef.alis.length; i++) 
			ll+= baliRef.alis[i].getSequences().length;
		String[] gfr= new String[ll];
		int[] gfrStart= new int[ll];
		String[] gfrNames= new String[ll];
		ll= 0;
		for (int i= 0; i < baliRef.alis.length; i++) 
			 for (int j= 0; j < baliRef.alis[i].getSequences().length; j++) {
				StringBuffer sb= new StringBuffer(baliRef.alis[i].getSequences()[j]);
				for (int k= 0; k < sb.length(); k++) // eliminate repeats
					if (sb.charAt(k)== '.')
						sb.deleteCharAt(k--);
				gfr[ll]= sb.toString();
				gfrNames[ll]= baliRef.alis[i].getSeqIDs()[j];
				gfrStart[ll++]= baliRef.alis[i].getStartPos()[j];
			}
		

		for (int i= 0; i < refTupels.length; i++) {		// sequences
			for (int j= 0; j < refTupels[i].length; j++) {			// mirror !		
				if (i==j)
					continue;
				
				AlignedTupel[] tTupel= nativeTupels[i][j]; 
				if (tTupel== null)
					continue;
				
				boolean noHit= true;
				for (int index= 0; noHit&& index < tTupel.length; index++) {
					AlignedTupel t= tTupel[index];
					for (int k= 0; noHit&& k < refTupels[i][j].length; k++) {	// iterate over repeats
						for (int l= 0; noHit&& l < refTupels[i][j][k].length; l++) { 
							AlignedTupel[] rTupel= refTupels[i][j][k][l];
							if (rTupel== null)
								continue;
							int x= Arrays.binarySearch(
								rTupel, t, new AlignedTupel.NaturalOrderComparator());	// t sorted?!
							if (x>= 0) 
								noHit= false;	// tupel exists
						}
					}
						
					if (!noHit)
						continue;
						
					//AlignedTupel r= rTupel[x];	// else..
					if (Character.isLowerCase(t.getCharA())||			// unaligned! Dialign truca! 
							Character.isLowerCase(t.getCharB()))
						continue;  

					int tA= t.getPositionA();
					int tB= t.getPositionB();
					boolean rA= isRepeatPosition(i,tA);
					boolean rB= isRepeatPosition(j,tB);
					if (rA^ rB) {
						xCountAll++;	
						
							// HERE! xCountCore
						int pos= rA?tA:tB;
						char c= rA?t.getCharA():t.getCharB();
						String seqName= rA?t.getSeqNameA():t.getSeqNameB();
						int k= 0;
						while (k< gfrNames.length) {
							
							for (; k < gfrNames.length; k++)		// find next seq# 
								if (gfrNames[k].startsWith(seqName))
									break; 
							if (k>= gfrNames.length)
								break;			// no more found
							
							if (pos< gfrStart[k]|| pos>= gfrStart[k]+ gfr[k].length()) {
								++k;
								continue;		// out of range, look for next match
							}
								
							if (Character.toLowerCase(gfr[k].charAt(pos- gfrStart[k]))!=	// repeat block found
									Character.toLowerCase(c))
								System.err.println("wrong repeat");	// check
							if (Character.isUpperCase(gfr[k].charAt(pos- gfrStart[k])))
								xCountCore++;
							break;	// only one match possible
							
						}
						if (k>= gfrNames.length)
							System.err.println("seqName not found "+seqName);
						
						
							
					}
				}
						
			}
		}
		
	}
	
	Repeat[][] repeatIdx= null;
	void initRepeatArray() {
		
			// init array
		Repeat[] tmpRepeats= baliRef.getRepeats();
		Vector[] sortRepeats= new Vector[seqNamesRef.length];
		for (int i= 0; i < sortRepeats.length; i++) 
			sortRepeats[i]= new Vector();
		
			// sort ascending repeats
		Repeat currRepeat;
		Repeat cmpRepeat;
		for (int i= 0; i < tmpRepeats.length; i++) {
			currRepeat= tmpRepeats[i];
			for (int j= 0; j < seqNamesRef.length; j++) 
				if (seqNamesRef[j].equalsIgnoreCase(currRepeat.getSeqName())) {
					currRepeat.setSeqNb(j);
					break;
				}
			if (currRepeat.getSeqNb()< 0)
				System.err.println("mapping error!");
			int j;
			for (j= 0; j < sortRepeats[currRepeat.getSeqNb()].size(); j++) {
				cmpRepeat= (Repeat) sortRepeats[currRepeat.getSeqNb()].elementAt(j);
				if (cmpRepeat.getStart()> currRepeat.getStart())
					break;
			}
			sortRepeats[currRepeat.getSeqNb()].insertElementAt(currRepeat, j);
		}
		
		repeatIdx= new Repeat[sortRepeats.length][];
		for (int i= 0; i < sortRepeats.length; i++) {
			repeatIdx[i]= new Repeat[sortRepeats[i].size()];
			for (int j= 0; j < sortRepeats[i].size(); j++) 
				repeatIdx[i][j]= (Repeat) sortRepeats[i].elementAt(j);
		}
	}
	
	boolean isRepeatPosition(int seqNb, int pos) {
		
			// iterate repeats
		Repeat inRepeat= null;
		for (int i= 0; i < repeatIdx[seqNb].length; i++) {
			if (pos< repeatIdx[seqNb][i].getStart())
				break;		// no more repeat pos
			if (pos>= repeatIdx[seqNb][i].getStart()+ repeatIdx[seqNb][i].getLength())
				continue;	// not in this repeat
			inRepeat= repeatIdx[seqNb][i];	// else: is in range of this repea
			break;
		}
		
//		if (!(pos>= inRepeat.getStart()&& pos< inRepeat.getStart()+ inRepeat.getLength()))
//			System.err.println("error");
		if (inRepeat!= null)
			return true;
		return false;
	}
	
	public double getAverageIdentity(boolean core) {

		refTupels= getRefTupels();
		refCoreTupelCount= getRefCoreTupelCount();
		testTupels= getTestTupels();
		
		int totalSim= 0;
		for (int i= 0; i < refTupels.length; i++) {
			for (int j= 0; j < refTupels[i].length; j++) {			// mirror !		
				if (i==j)
					continue;
				for (int k= 0; k < refTupels[i][j].length; k++) {
					for (int l= 0; l < refTupels[i][j][k].length; l++) {
						totalSim+= getSimilarity(testTupels[i][j][k][l], refTupels[i][j][k][l], core);
					}
				}
			}
		}
		
		double divisor= 0d;
		if (core) 
			divisor= (double) getTotalHitCountCore();
		else
			divisor= (double) getTotalHitCountAll();
		
		divisor= countAllTestTupels();
		double res= ((double) totalSim)/ divisor;
		
		return res;
	}
	
	int countAllTestTupels() {
		getTestTupels();
		int tup= 0;
		for (int i= 0; i < testTupels.length; i++) {
			for (int j= 0; j < testTupels[i].length; j++) {
				if (i== j)
					continue;
				for (int k= 0; testTupels[i][j]!= null&& k < testTupels[i][j].length; k++) {
					for (int l= 0; testTupels[i][j][k]!= null&& l < testTupels[i][j][k].length; l++) {
						tup+= testTupels[i][j][k][l].length;
					}
				}
			}
		}
		
		return tup;
	}
	
	int getSimilarity(AlignedTupel[] tTupel, AlignedTupel[] rTupel, boolean core) {
		
		if (tTupel== null|| rTupel== null)
			return 0;
		
			// assume sorted arrays! 
		int sim= 0;
		int[] mapping= Constants.getMapBChars();
		for (int i= 0; i < tTupel.length; i++) {
			AlignedTupel t= tTupel[i];
			int x= Arrays.binarySearch(
				rTupel, t, new AlignedTupel.NaturalOrderComparator());
			if (x< 0) 
				continue;
			AlignedTupel r= rTupel[x];
			
			char r_a= Character.toUpperCase(r.getCharA());
			char t_a= Character.toUpperCase(t.getCharA());
			char r_b= Character.toUpperCase(r.getCharB());
			char t_b= Character.toUpperCase(t.getCharB());
			 
			if ((r_a!= t_a)|| (r_b!= t_b))		// test
				System.err.println("Characters dont match:\n\t"+t+"\n\t"+r);
			if(!core|| (core&& r.isCore()))
				try {
					sim+= Constants.getBlosum()
						[mapping[r_a- 65]]
						[mapping[r_b- 65]];
				} catch (ArrayIndexOutOfBoundsException e) {
					; // not supported char
				}
		}
		
		return sim;
		
	}	
	
	public int[][][][] getRepeatProfile() {
		if (repProfile == null) {
			
			testTupels= getTestTupels();
			repProfile= new int[testTupels.length][][][];
			for (int i= 0; i < testTupels.length; i++) {
				repProfile[i]= new int[testTupels[i].length][][];
				for (int j= 0; j < testTupels[i].length; j++) {
					if (i== j)
						continue;
					repProfile[i][j]= new int[testTupels[i][j].length][];
					for (int k= 0; k < testTupels[i][j].length; k++) {
						Vector alignedRepeats= new Vector();
						for (int l= 0; l < testTupels[i][j][k].length; l++) {
							if (testTupels[i][j][k][l].length> 0)
								alignedRepeats.addElement(new Integer(l));
						}
						repProfile[i][j][k]= new int[alignedRepeats.size()];
						for (int l= 0; l < repProfile[i][j][k].length; l++) 
							repProfile[i][j][k][l]= ((Integer) alignedRepeats.elementAt(l)).intValue();
					}
				}
			}
		}

		return repProfile;
	}
		
	
	public PSComparator(String newFileBase, String newFileExt) {
		this.fileBase= newFileBase;
		baliRef= new BAliManager(fileBase);
		try {
			initAlignment(newFileExt);
		} catch (Exception e) {
			aliTest= null;			
		}
	}
	

	
	public void initAlignment(String newFileExt) throws Exception {
		 this.aliTest= readAlignment(newFileExt);
	}
	
	public Alignment readAlignment(String newFileExt)  throws Exception {
		String readAli= this.fileBase+newFileExt;
		MSFWrapper msf= new MSFWrapper(readAli);
		try {
			msf.read();
		} catch (Exception e) {
			System.err.println("Could not read test alignment "+readAli+"!");
			throw new Exception();
		}
		Alignment ali= new Alignment(msf.getSequences());
		ali.setSeqIDs(msf.getSeqNames());
		
		return ali;		
	}
	
	public int getTotalHitCountAll() {
		
		if (totalHitCountAll< 0) {
			
			seq2HitCountAll= getSeq2HitCountAll();		// init pw hit counts
			
			totalHitCountAll= 0;						// add up total hits
			for (int i= 0; i < seq2HitCountAll.length; i++) 
				for (int j= (i+1); j < seq2HitCountAll[i].length; j++)	
					totalHitCountAll+=  seq2HitCountAll[i][j];	// ONLY one triangle (otherwise double hits)
		}

		return totalHitCountAll;
	}
	
	public int getTotalHitCountCore() {
		
		if (totalHitCountCore< 0) {
			
			seq2HitCountCore= getSeq2HitCountCore();		// init pw hit counts
			
			totalHitCountCore= 0;						// add up total hits
			for (int i= 0; i < seq2HitCountCore.length; i++) 
				for (int j= (i+1); j < seq2HitCountCore[i].length; j++)	
					totalHitCountCore+=  seq2HitCountCore[i][j];	// ONLY one triangle (otherwise double hits)
		}

		return totalHitCountCore;
	}	
	
	public int[][] getSeq2HitCountAll() {

		if (seq2HitCountAll == null) {
			
			hitCountAll= getHitCountAll();			// init hit count
			
			seq2HitCountAll= new int[hitCountAll.length][hitCountAll.length];
			for (int i= 0; i < hitCountAll.length; i++) {
				for (int j= 0; j < hitCountAll[i].length; j++) {
					
					seq2HitCountAll[i][j]= 0;
					if (i== j) 					// equal seqs: 0
						continue;
					
					for (int k= 0; k < hitCountAll[i][j].length; k++) 	// iterate over repeats
						for (int l= 0; l < hitCountAll[i][j][k].length; l++) 
							seq2HitCountAll[i][j]+= hitCountAll[i][j][k][l];
				}
			}
			
		}

		return seq2HitCountAll;
	}
	public int[][] getSeq2HitCountCore() {

			if (seq2HitCountCore == null) {
			
				hitCountCore= getHitCountCore();			// init hit count
			
				seq2HitCountCore= new int[hitCountCore.length][hitCountCore.length];
				for (int i= 0; i < hitCountCore.length; i++) {
					for (int j= 0; j < hitCountCore[i].length; j++) {
					
						seq2HitCountCore[i][j]= 0;
						if (i== j) 					// equal seqs: 0
							continue;
					
						for (int k= 0; k < hitCountCore[i][j].length; k++) 	// iterate over repeats
							for (int l= 0; l < hitCountCore[i][j][k].length; l++) 
								seq2HitCountCore[i][j]+= hitCountCore[i][j][k][l];
					}
				}
			
			}

			return seq2HitCountCore;
		}	
	
	public int[][][][] getHitCountAll() {

		if (hitCountAll == null) {
			
			refTupels= getRefTupels();
			refCoreTupelCount= getRefCoreTupelCount();
			testTupels= getTestTupels();
			hitCountAll= new int[refTupels.length][refTupels.length][][];
			countHits(hitCountAll, false);
			
		}

		return hitCountAll;
	}
	
	public int[][][][] getHitCountCore() {

		if (hitCountCore == null) {
			
			refTupels= getRefTupels();
			refCoreTupelCount= getRefCoreTupelCount();
			testTupels= getTestTupels();
			hitCountCore= new int[refTupels.length][refTupels.length][][];
			countHits(hitCountCore, true);		// only look for core matches
			
		}

		return hitCountCore;
	}	
	
	void countHits(int[][][][] hitCount, boolean core) {

		for (int i= 0; i < refTupels.length; i++) {
			for (int j= 0; j < refTupels[i].length; j++) {			// mirror !		
				if (i==j)
					continue;
				hitCount[i][j]= new int[refTupels[i][j].length][];
				for (int k= 0; k < refTupels[i][j].length; k++) {
					hitCount[i][j][k]= new int[refTupels[i][j][k].length];
					for (int l= 0; l < refTupels[i][j][k].length; l++) {
						hitCount[i][j][k][l]= countHits(testTupels[i][j][k][l], refTupels[i][j][k][l], core);
					}
				}
			}
		}
	}
	
	public double getAvgHitPercent(boolean core) {
		
		repProfile= getRepeatProfile();


		int hits= 0;
		int refs= 0;
		for (int i= 0; i < repProfile.length; i++) {				// iterate seq pairs
			for (int j= (i+1); j < repProfile[i].length; j++) {		// only one triangle
				for (int k= 0; k < repProfile[i][j].length; k++) {	// all repeats in A
					for (int l= 0; l < repProfile[i][j][k].length; l++) { // rep nbs in B
						if (core) {
							hits+= hitCountCore[i][j][k][repProfile[i][j][k][l]];
							refs+= refCoreTupelCount[i][j][k][repProfile[i][j][k][l]];
						} else {
							hits+= hitCountAll[i][j][k][repProfile[i][j][k][l]];
							refs+= refTupels[i][j][k][repProfile[i][j][k][l]].length;
						}
						if (hits!= refs)
							System.currentTimeMillis();
					}
				}
			}
		}

		double avgHitPercent= (double) hits/ (double) refs;
		return avgHitPercent;
	}
	
	public double getMinMaxHitPercent(boolean core, boolean min) {
		
		repProfile= getRepeatProfile();

		double ratio= 0d; 	// error extrapolation, I know, better to save hits and refs..doit
		for (int i= 0; i < repProfile.length; i++) {				// iterate seq pairs
			for (int j= (i+1); j < repProfile[i].length; j++) {		// only one triangle
				for (int k= 0; k < repProfile[i][j].length; k++) {	// all repeats in A
					double bestRatio= 0d;
					if (min)
						bestRatio= Double.MAX_VALUE;
					for (int l= 0; l < repProfile[i][j][k].length; l++) { // rep nbs in B
						if (core) 
							ratio= (double) hitCountCore[i][j][k][repProfile[i][j][k][l]] /
									(double) refCoreTupelCount[i][j][k][repProfile[i][j][k][l]];
						else 
							ratio= (double) hitCountAll[i][j][k][repProfile[i][j][k][l]] /
									(double) refTupels[i][j][k][repProfile[i][j][k][l]].length;
						if (min)
							bestRatio= Math.min(bestRatio, ratio);
						else
							bestRatio= Math.max(bestRatio, ratio);
					}
					ratio+= bestRatio;
				}
			}
		}

		return ratio;
	}	
	
	int countHits(AlignedTupel[] tTupel, AlignedTupel[] rTupel, boolean core) {
		
		if (tTupel== null|| rTupel== null)
			return 0;
		
			// assume sorted arrays! 
		int hits= 0;
		for (int i= 0; i < tTupel.length; i++) {
			AlignedTupel t= tTupel[i];
			int x= Arrays.binarySearch(
				rTupel, t, new AlignedTupel.NaturalOrderComparator());
			if (x< 0) 
				continue;
			AlignedTupel r= rTupel[x];
			
			if (Character.toLowerCase(r.getCharA())!= Character.toLowerCase(t.getCharA())|| 
					Character.toLowerCase(r.getCharB())!= Character.toLowerCase(t.getCharB()))		// test
				System.err.println("Characters dont match:\n\t"+t+"\n\t"+r);
			if(!core|| (core&& r.isCore()))
				++hits;
		}
		
		return hits;
	}
	
	
	public void compareTestRef() {
		
		getRefTupels();
		getTestTupels();
		
			// get test
		Vector allAll= new Vector();
		for (int i= 0; i < testTupels.length; i++) {
			for (int j= (i+1); j < testTupels[i].length; j++) {
				for (int k= 0; k < testTupels[i][j].length; k++) {
					for (int l= (k+1); l < testTupels[i][j][k].length; l++) {
						for (int m= 0; m < testTupels[i][j][k][l].length; m++) {
							allAll.add(testTupels[i][j][k][l][m]);
						}
					}
				}
			}
		}
		AlignedTupel[] tTupels= new AlignedTupel[allAll.size()];
		for (int i= 0; i < allAll.size(); i++) 
			tTupels[i]= (AlignedTupel) allAll.elementAt(i);
		Arrays.sort(tTupels, new AlignedTupel.NaturalOrderComparator());
		
			// get reference
		allAll= new Vector();
		for (int i= 0; i < refTupels.length; i++) {
			for (int j= (i+1); j < refTupels[i].length; j++) {
				for (int k= 0; k < refTupels[i][j].length; k++) {
					for (int l= (k+1); l < refTupels[i][j][k].length; l++) {
						for (int m= 0; m < refTupels[i][j][k][l].length; m++) {
							allAll.add(refTupels[i][j][k][l][m]);
						}
					}
				}
			}
		}
		AlignedTupel[] rTupels= new AlignedTupel[allAll.size()];
		for (int i= 0; i < allAll.size(); i++) 
			rTupels[i]= (AlignedTupel) allAll.elementAt(i);
		Arrays.sort(rTupels, new AlignedTupel.NaturalOrderComparator());
		
		
			// compare
		for (int i= 0; i < rTupels.length; i++) {
			
			int idx= Arrays.binarySearch(tTupels, rTupels[i], new AlignedTupel.NaturalOrderComparator());
			if (idx< 0)
				System.err.println("Not found in test tupels: "+ rTupels[i]);
		}
	}
	
	
	AlignedTupel[][][] nativeTupels;
	public AlignedTupel[][][][][] getTestTupels() {
		
		if (testTupels == null) {
						
				// get tupels
			nativeTupels= 
				aliTest.getAlignedTupels();	
			
				// enforce alphbet. order of sequences
			String[] seqNames= aliTest.getSeqIDs();
			seqNamesTest= new String[seqNames.length];
			for (int i= 0; i < seqNamesTest.length; i++) 
				seqNamesTest[i]= seqNames[i].toString();	// no clone?!
			Arrays.sort(seqNamesTest);
			int[] reorder= new int[seqNamesTest.length];
			for (int i= 0; i < seqNames.length; i++)			// reorder[source][dest] 
				for (int j= 0; j < seqNamesTest.length; j++) 
					if (seqNames[i].equals(seqNamesTest[j])) {
						reorder[i]= j;
						break;
					}
			for (int i= 0; i < reorder.length; i++) {
				
				if (reorder[i]== i)
					continue;
				
				AlignedTupel[][] tmp= nativeTupels[reorder[i]];	// swap primary dimension
				nativeTupels[reorder[i]]= nativeTupels[i];
				nativeTupels[i]= tmp;
				
				for (int j= 0; j < nativeTupels.length; j++) {	// swap secondary dimensions
					AlignedTupel[] temp= nativeTupels[j][i];
					nativeTupels[j][i]= nativeTupels[j][reorder[i]];
					nativeTupels[j][reorder[i]]= temp;
				}

				int t= reorder[reorder[i]];		// swap in reorder[]
				reorder[reorder[i]]= reorder[i];
				reorder[i]= t;
				if (reorder[i]!= i)		// still not in final position
					i--;				// swap again
			}					
			for (int i= 0; i < nativeTupels.length; ++i){		// test & chk same seq names 
				String nameTest= seqNamesTest[i];
				String nameRef= seqNamesRef[i];
				if (!nameTest.toLowerCase().startsWith(nameRef.toLowerCase())
					&& !nameRef.toLowerCase().startsWith(nameTest.toLowerCase()))
					System.err.println("Seq mismatch in read files : "+
						nameTest+" (test) <> "+nameRef+" (ref)");
			}
				
			
				// order according to repeats
			testTupels= new AlignedTupel[nativeTupels.length][nativeTupels.length][][][];
			for (int i= 0; i < nativeTupels.length; i++) {
				for (int j= 0; j < nativeTupels[i].length; j++) {
					if (i== j)
						continue;
					testTupels[i][j]= getTestTupels(
						nativeTupels[i][j],
						seqNamesRef[i],
						seqNamesRef[j]
					); 
									// mirror possible
				}
			}
		}
			
		return testTupels;
	}
	
	AlignedTupel[][][] getTestTupels(AlignedTupel[] nativeTupels, String seqA, String seqB) {
		
			// retrieve repeats for the sequences
		Repeat[] repeats= baliRef.getRepeats();
	
			// get repeats for the two sequences A and B
		Vector seqARepVec= new Vector();
		Vector seqBRepVec= new Vector();
		for (int i= 0; i < repeats.length; i++) {
			if (repeats[i].getSeqName().toLowerCase().startsWith(seqA.toLowerCase()))
				seqARepVec.add(repeats[i]);
			if (repeats[i].getSeqName().toLowerCase().startsWith(seqB.toLowerCase()))
				seqBRepVec.add(repeats[i]);
		}
		Repeat[] seqARep= new Repeat[seqARepVec.size()];
		for (int i= 0; i < seqARepVec.size(); i++) 
			seqARep[i]= (Repeat) seqARepVec.elementAt(i);
		Repeat[] seqBRep= new Repeat[seqBRepVec.size()];
		for (int i= 0; i < seqBRepVec.size(); i++) 
			seqBRep[i]= (Repeat) seqBRepVec.elementAt(i);
			
			// for all repeat pairs get aligned tupels
		AlignedTupel[][][] result= new AlignedTupel[seqARep.length][seqBRep.length][];
		for (int i= 0; i < seqARep.length; i++) {
			for (int j= 0; j < seqBRep.length; j++) {
				result[i][j]= getTestTupels(seqARep[i], seqBRep[j], nativeTupels);
			}
		}
		
		return result;			 						
	}
	
	
	AlignedTupel[] getTestTupels(Repeat repA, Repeat repB, AlignedTupel[] nativeTupels) {
		
		Vector resultVec= new Vector();
		for (int i= 0; i < nativeTupels.length; i++) 
			if ((repA.contains(nativeTupels[i].getPositionA()))
				&& (repB.contains(nativeTupels[i].getPositionB())))	// lets assume sequence order is ok..
				resultVec.add(nativeTupels[i]);
		
		AlignedTupel[] res= new AlignedTupel[resultVec.size()];
		for (int i= 0; i < res.length; i++) {
			res[i]= (AlignedTupel) resultVec.elementAt(i);
		}
		Arrays.sort(res, new AlignedTupel.NaturalOrderComparator());

		return res;
	}
		
	/**
	 * @return
	 */
	public int[][][][] getRefCoreTupelCount() {

		if (refCoreTupelCount == null) 
			refCoreTupelCount= baliRef.getCoreTupelCount();
		
		return refCoreTupelCount;
	}

	/**
	 * @return
	 */
	public AlignedTupel[][][][][] getRefTupels() {

		if (refTupels == null) {
			refTupels= baliRef.getAlignedTupels();
			seqNamesRef= baliRef.getSeqNames();
		}

		return refTupels;
	}
	
}
