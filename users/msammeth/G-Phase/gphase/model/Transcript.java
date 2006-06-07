/*
 * Created on Mar 3, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.Constants;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * 
 * 
 * @author micha
 */
public class Transcript extends DirectedRegion {

	static final long serialVersionUID = 2863324463934791891L;
	int type = Constants.NOINIT;
	int confidence = Constants.NOINIT;
		
	
	/**
	 * 
	 * @param newType type descriptor
	 * @return <code>false</code> if the type is already set or the type
	 * descriptor is not valid, <code>true</code> otherwise.  
	 */
	public boolean setType(String newType) {
		
		if (newType== null|| type!= Constants.NOINIT) 
			return false;
		
		int i= Constants.findIgnoreCase(newType, Constants.TYPES);
		if (i< 0) {
			System.err.println("Unknown type "+ newType);
			return false;
		}
		type= i;
		return true;
	}
	
	public String toString() {
		return getTranscriptID();
	}

	/**
	 * 
	 * @param newConfidence confidence descriptor
	 * @return <code>false</code> if confidence is already set or the confidence
	 * descriptor is not valid, <code>true</code> otherwise.  
	 */
	public boolean setConfidence(String newConfd) {
		
		if (newConfd== null|| type!= Constants.NOINIT)
			return false;
		
		int i= Constants.findIgnoreCase(newConfd, Constants.CONFIDENCES);
	
		if (i< 0) {
			System.err.println("Unknown confidence "+ newConfd);
			return false;
		}
	
		confidence= i;
		return true;
	}
	Translation[] translations= null;
	SpliceSite[] spliceChain= null;	// sorted!!
	public final static Transcript[] toTranscriptArray(Vector v) {
		if (v== null)
			return null;
		Transcript[] result= new Transcript[v.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Transcript) v.elementAt(i);
		return result;
	}	
	
	/**
	 * 
	 * @return the first translation registered for <code>this</code> transcript.
	 */
	public Translation getDefaultTranslation() {
		if (translations== null) {
			translations= new Translation[]{new Translation(this)};
		}
		return translations[0];
	}
	
	public SpliceSite[] getSpliceChain(){
		
		spliceChain= null;
		if (spliceChain== null) {
			if (exons== null)
				return null;
			
			spliceChain= new SpliceSite[exons.length* 2- 2];
			int pos= 0;
			for (int i = 0; i < exons.length; i++) {
				if(i> 0)
					spliceChain[pos++]= exons[i].getAcceptor();
				if(i< exons.length- 1)
					spliceChain[pos++]= exons[i].getDonor();
			}
			
			Arrays.sort(spliceChain, new AbstractSite.PositionComparator());
		}
		return spliceChain;
	}
	
	public SpliceSite getPredSpliceSite(SpliceSite s) {
		
		SpliceSite[] ss= getSpliceChain();
		int p= Arrays.binarySearch(ss, s, new SpliceSite.PositionComparator());
		if ((p< 0)|| (p== 0)) 
			return null;
		return ss[p-1];
	}
	
	public SpliceSite getSuccSpliceSite(SpliceSite s) {
		
		SpliceSite[] ss= getSpliceChain();
		int p= Arrays.binarySearch(ss, s, new SpliceSite.PositionComparator());
		if ((p< 0)|| (p> (ss.length-2))) 
			return null;
		return ss[p+1];
	}
	
	public SpliceSite getSuccSpliceSite(int pos) {
		SpliceSite[] ss= getSpliceChain();
		return getSuccSpliceSite(ss, pos);
	}
	public static SpliceSite getSuccSpliceSite(SpliceSite[] ss, int pos) {
		int p= Arrays.binarySearch(ss, new Integer(pos+1), new AbstractSite.PositionToSpliceSiteComparator());	// abstract site for only comparing pos
		if(p>= 0) {	// found ?!
			return ss[p];
		}
		p= -(p+ 1);	// insertionpoint= point of successor
		if (p< ss.length&& ss[p].getPos()> pos)
			return ss[p];
		else
			return null;
	}
	
	public static int getSuccPos(int[] ss, int pos) {
		int p= Arrays.binarySearch(ss, pos+1);	// abstract site for only comparing pos
		if(p>= 0) {
			return ss[p];
		}
		p= -(p+ 1);	// insertionpoint= point of successor
		if (p< ss.length&& ss[p]> pos)
			return ss[p];
		else
			return 0;
	}
	
	public SpliceSite getPredSpliceSite(int pos) {
		SpliceSite[] ss= getSpliceChain();
		return getPredSpliceSite(ss, pos);
	}
	
	public static SpliceSite getPredSpliceSite(SpliceSite[] ss, int pos) {
		int p= Arrays.binarySearch(ss, new Integer(pos-1), new AbstractSite.PositionToSpliceSiteComparator());	// abstract site for only comparing pos
		if(p>= 0) {	// found ?!
			return ss[p];
		}
		p= -(p+ 1)- 1;	// insertion point- 1= point of predecessor
		if (p>= 0&& p< ss.length&& ss[p].getPos()< pos)
			return ss[p];
		else
			return null;
	}
	
	public static int getPredPos(int[] ss, int pos) {
		int p= Arrays.binarySearch(ss, pos-1);	// abstract site for only comparing pos
		if(p>= 0) {
			return ss[p];
		}
		p= -(p+ 1)- 1;	// insertion point- 1= point of predecessor
		if (p>= 0&& p< ss.length&& ss[p]< pos)
			return ss[p];
		else
			return 0;
	}
	
	public static SpliceSite getSpliceSiteByPos(SpliceSite[] ss, int pos) {
		int p= Arrays.binarySearch(ss, new Integer(pos), new AbstractSite.PositionToSpliceSiteComparator());	// abstract site for only comparing pos
		if(p>= 0) 
			return ss[p];
		return null;
	}
	
	public static int getPos(int[] ss, int pos) {
		int p= Arrays.binarySearch(ss, pos);	// abstract site for only comparing pos
		if(p>= 0) 
			return ss[p];
		return -1;
	}
	
	public boolean addTranslation(Translation newTrans) {
		
			// new transcipt array
		if (translations== null) {
			translations= new Translation[] {newTrans};
			return true;
		}
			
			// search transcript for same translation, not necessary
			
			// add translation
		Translation[] nTranslations= new Translation[translations.length+ 1];
		for (int i= 0; i < translations.length; i++) 
			nTranslations[i]= translations[i];
		nTranslations[nTranslations.length- 1]= newTrans;
		translations= nTranslations;
		return true;
	}
	
	public Translation[] getTranslation(int start, int end) {
		if (translations== null)
			return null;
		Vector result= new Vector();
		for (int i = 0; i < translations.length; i++) {
			if ((translations[i].getStart()<= start)&&
				(translations[i].getEnd()>= end))
				result.add(translations[i]);
		}
		
		Object o= gphase.tools.Arrays.toField(result);
		if (o== null|| result.size()< 1)
			return null;
		return (Translation[]) o;
	}
	
	public Translation[] getTranslations() {
		return translations;
	}
	
	public boolean isNonCoding() {
		return (translations== null|| translations.length== 0);
	}
	
	public boolean is5UTR(int pos) {
		
		int i;
		for (i = 0; i < translations.length; i++) {
			if (getGene().isForward()&& translations[i].getStart()<= pos)
				continue;
			if (!getGene().isForward()&& translations[i].getEnd()> pos)		// reverse strand: rev coordinates
				continue;
			break;
		}
		if (i< translations.length)
			return true;
		return false;
	}

	/**
	 * when 
	 * @param pos
	 * @return
	 */
	public boolean isCDS(int pos) {
		
		if (translations== null)
			return false;
		
		int i;
		for (i = 0; i < translations.length; i++) 
			if (pos>= translations[i].getStart()&& pos<= translations[i].getEnd())
				break;
		
		if (i< translations.length)
			return true;
		return false;
	}
	
	public boolean isForward() {
		if (gene!= null)
			return gene.isForward();
		return super.isForward();
	}
	
	Gene gene= null;

	String transcriptID= null;
	Exon[] exons= null;	// sorted !!
	public Transcript(Gene newGene, String stableTranscriptID) {

		this.strand= newGene.getStrand();
		this.gene= newGene;
		this.transcriptID= stableTranscriptID;
	}
	
	public Transcript(String newID) {
		this.transcriptID= newID;
	}
	
	/**
	 * @return
	 */
	public Exon[] getExons() {
		return exons;
	}
	
	/**
	 * @return
	 */
	public Gene getGene() {
		return gene;
	}

	/**
	 * @return
	 */
	public String getTranscriptID() {
		return transcriptID;
	}
	
	/**
	 * @return
	 */
	public String getStableID() {
		
		return transcriptID;
	}	

	/**
	 * @param exons
	 */
	public void setExons(Exon[] exons) {
		this.exons= exons;
	}

	/**
	 * @deprecated mandatory in constructor now
	 * @param gene
	 */
	public void setGene(Gene gene) {
		this.gene= gene;
	}

	/**
	 * @deprecated not valid for HOX-genes, AY... genes, ..
	 * @param i
	 */
	public void setTranscriptID(String i) {
		transcriptID= i;
	}

	public String getChromosome() {
		return getGene().getChromosome();
	}
	
	public Species getSpecies() {
		return getGene().getSpecies();
	}
	
	/**
	 * @param b
	 */
	public boolean checkStrand(boolean b) {
		
		return (b== getGene().isForward());
	}
	
	public boolean checkStrand(String newStrand) {
		
		String nStrand= newStrand.trim();
		if (nStrand.equals("1")) 	// || nStrand.equals("forward")
			return checkStrand(true); 
		else if (nStrand.equals("-1"))	// || nStrand.equals("reverse")
			return checkStrand(false);
		
		return false; // error			
	}


	public void addCDS(int start, int end) {
		if (!isForward()) {
			start= -start;
			end= -end;
		}
			
		if (translations== null) {
			translations= new Translation[] {new Translation(this)};
			translations[0].setStart(start);
			translations[0].setEnd(end);
			return;
		}
			// else
		if (start< translations[0].getStart())
			translations[0].setStart(start);
		if (end> translations[0].getEnd())
			translations[0].setEnd(end);
	}
	
	/**
		 * Inserts the exons in an array sorted according to ascending order
		 * of their start/stop position. <b>IMPORTANT</b>: add exons AFTER adding 
		 * transcripts to ensure the correct init of AS types.
		 * 
		 * @param newExon
		 * @return the exon already contained or <code>newExon</code> case of the exon was added successfully
		 */
		public boolean addExon_new(Exon newExon) {
	
				// new exon array
			if (exons== null) 
				exons= new Exon[] {newExon};
			else {
				
					// search for identical exon (HERE necessary)
				int p= Arrays.binarySearch(
						exons,
						newExon,
						new AbstractRegion.PositionComparator()	
					);
		
				if (p>= 0) 
					return false;	// already contained, not added
				
					// new Exon: search for overlapping exons 
					// and accordingly construct SuperExon
	//			Exon[] exs= getGene().getExons();		// in ALL transcripts
	//			for (int i = 0; i < exons.length; i++) {
	//				if ((exs[i].getStart()>= newExon.getStart()&& exs[i].getStart()< newExon.getEnd())||
	//						newExon.getStart()>= exs[i].getStart()&& newExon.getStart()< exs[i].getEnd()) { // intersecting
	//					if (exs[i].getSuperExon()!= null)
	//						if (newExon.getSuperExon()!= null)
	//							newExon.getSuperExon().merge(exs[i].getSuperExon()); 	// merge two super-exs
	//						else
	//							exs[i].getSuperExon().add(newExon);	// add to super-exon of the other
	//					else {
	//						SuperExon superEx= newExon.getSuperExon();
	//						if (superEx== null) {
	//							superEx= new SuperExon();
	//							superEx.add(newExon);
	//						}
	//						superEx.add(exs[i]);
	//					}
	//				}
	//			}
				
					// add exon
				exons= (Exon[]) gphase.tools.Arrays.insert(this.exons, newExon, p);
			}
	
				// init splice sites 
			int p= Arrays.binarySearch(
					this.exons,
					newExon,
					new AbstractRegion.PositionComparator()	
			);
			
			if ((p== 0)&& exons.length>1) {	// ex-first exon now has an acceptor
				Exon ex= this.exons[1];
				int posAcceptor= isForward()?ex.getStart()-2:ex.getEnd()+1;
				SpliceSite acceptor= new SpliceSite(getGene(), posAcceptor, false);
				SpliceSite ss= getGene().checkSpliceSite(acceptor);
				if (ss== null) 
					getGene().addSpliceSite(acceptor);
				else
					acceptor= ss;
				ex.setAcceptor(acceptor);
			}
			
			if (p> 0) {										// has acceptor
				int posAcceptor= isForward()?newExon.getStart()-2:newExon.getEnd()+1;
				SpliceSite acceptor= new SpliceSite(getGene(), posAcceptor, false);
				SpliceSite ss= getGene().checkSpliceSite(acceptor);
				if (ss== null) 
					getGene().addSpliceSite(acceptor);
				else
					acceptor= ss;
				exons[p].setAcceptor(acceptor);
			}
			
			if (p< this.exons.length- 1) {					// has donor
				int posDonor= isForward()?newExon.getEnd()+1:newExon.getStart()-2;
				SpliceSite donor= new SpliceSite(getGene(), posDonor, true);
				SpliceSite ss= getGene().checkSpliceSite(donor);
				if (ss== null) 
					getGene().addSpliceSite(donor);
				else
					donor= ss;
				exons[p].setDonor(donor);
			}
	
			if ((p== this.exons.length- 1)&& exons.length>1) {	// ex-last exon now has an donor
				Exon ex= this.exons[exons.length- 2];
				int posDonor= isForward()?ex.getEnd()+1:ex.getStart()-2;
				SpliceSite donor= new SpliceSite(getGene(), posDonor, true);
				SpliceSite ss= getGene().checkSpliceSite(donor);
				if (ss== null) 
					getGene().addSpliceSite(donor);
				else
					donor= ss;
				ex.setDonor(donor);
			}
			
			return true;
		}

		/**
		 * Inserts the exons in an array sorted according to ascending order
		 * of their start/stop position. <b>IMPORTANT</b>: add exons AFTER adding 
		 * transcripts to ensure the correct init of AS types.
		 * 
		 * @param newExon
		 * @return the exon already contained or <code>newExon</code> case of the exon was added successfully
		 */
		public boolean addExon(Exon newExon) {
	
				// new exon array
			if (exons== null) 
				exons= new Exon[] {newExon};
			else {
				
					// search for identical exon (HERE necessary)
				int p= Arrays.binarySearch(
						exons,
						newExon,
						new AbstractRegion.PositionComparator()	// has to be directed, end< start for two ident. regions	
					);
		
				if (p>= 0) 
					return false;	// already contained, not added
				
					// new Exon: search for overlapping exons 
					// and accordingly construct SuperExon
	//			Exon[] exs= getGene().getExons();		// in ALL transcripts
	//			for (int i = 0; i < exons.length; i++) {
	//				if ((exs[i].getStart()>= newExon.getStart()&& exs[i].getStart()< newExon.getEnd())||
	//						newExon.getStart()>= exs[i].getStart()&& newExon.getStart()< exs[i].getEnd()) { // intersecting
	//					if (exs[i].getSuperExon()!= null)
	//						if (newExon.getSuperExon()!= null)
	//							newExon.getSuperExon().merge(exs[i].getSuperExon()); 	// merge two super-exs
	//						else
	//							exs[i].getSuperExon().add(newExon);	// add to super-exon of the other
	//					else {
	//						SuperExon superEx= newExon.getSuperExon();
	//						if (superEx== null) {
	//							superEx= new SuperExon();
	//							superEx.add(newExon);
	//						}
	//						superEx.add(exs[i]);
	//					}
	//				}
	//			}
				
					// add exon
				exons= (Exon[]) gphase.tools.Arrays.insert(this.exons, newExon, p);
			}
	
				// init splice sites 
			int p= Arrays.binarySearch(
					this.exons,
					newExon,
					new AbstractRegion.PositionComparator()		// here also abstractregion possible	
			);
			
			if (p== 0) {
				int tss= isForward()?newExon.getStart():newExon.getEnd();
				setStart(tss);
				if (exons.length>1) {	// ex-first exon now has an acceptor
					Exon ex= this.exons[1];
					int posAcceptor= isForward()?ex.getStart():ex.getEnd();
					SpliceSite acceptor= new SpliceSite(getGene(), posAcceptor, false, ex);
					SpliceSite ss= getGene().checkSpliceSite(acceptor);
					if (ss== null) 
						getGene().addSpliceSite(acceptor);
					else {
						acceptor= ss;
						ss.addExon(ex);
					}
					ex.setAcceptor(acceptor);
				}
			}
			
			if (p> 0) {										// has acceptor
				int posAcceptor= isForward()?newExon.getStart():newExon.getEnd();
				SpliceSite acceptor= new SpliceSite(getGene(), posAcceptor, false, newExon);
				SpliceSite ss= getGene().checkSpliceSite(acceptor);
				if (ss== null) 
					getGene().addSpliceSite(acceptor);
				else {
					acceptor= ss;
					ss.addExon(newExon);
				}
				exons[p].setAcceptor(acceptor);
			}
			
			if (p< this.exons.length- 1) {					// has donor
				int posDonor= isForward()?newExon.getEnd():newExon.getStart();
				SpliceSite donor= new SpliceSite(getGene(), posDonor, true, newExon);
				SpliceSite ss= getGene().checkSpliceSite(donor);
				if (ss== null) 
					getGene().addSpliceSite(donor);
				else {
					donor= ss;
					ss.addExon(newExon);
				}
				exons[p].setDonor(donor);
			}
	
			if (p== this.exons.length- 1) {
				int tes= isForward()?newExon.getEnd():newExon.getStart();
				setEnd(tes);
				if (exons.length>1) {	// ex-last exon now has an donor
					Exon ex= this.exons[exons.length- 2];
					int posDonor= isForward()?ex.getEnd():ex.getStart();
					SpliceSite donor= new SpliceSite(getGene(), posDonor, true, ex);
					SpliceSite ss= getGene().checkSpliceSite(donor);
					if (ss== null) 
						getGene().addSpliceSite(donor);
					else {
						donor= ss;
						ss.addExon(ex);
					}
					ex.setDonor(donor);
				}
			}
			
			return true;
		}

		/**
		 * Inserts the exons in an array sorted according to ascending order
		 * of their start/stop position. <b>IMPORTANT</b>: add exons AFTER adding 
		 * transcripts to ensure the correct init of AS types.
		 * 
		 * @param newExon
		 * @return the exon already contained or <code>newExon</code> case of the exon was added successfully
		 */
		public void setBoundaries(Exon newExon) {
	
			if (getStart()== 0|| newExon.getStart()< getStart())
				setStart(newExon.getStart());
			if (getEnd()== 0|| newExon.getEnd()> getEnd())
				setEnd(newExon.getEnd());
		}

	/**
	 * Finds the first exon containing the corresponding position
	 * 
	 * @deprecated inconsistent for overlapping exons		 
	 * @param absPos
	 * @return
	 */
	public Exon getExon(int absPos) {
		
		for (int i = 0; i < exons.length; i++) 
			if (exons[i].contains(absPos))
				return exons[i];
		
		return null;
	}

	public int getOtherSideOfExon(int pos) {
		for (int i = 0; i < exons.length; i++) {
			if (exons[i].getStart()== pos)
				return exons[i].getEnd();
			else if (exons[i].getEnd()== pos)
				return exons[i].getStart();
		}
		
		return -1;
	}
	
	public int getTSSPos() {
		return exons[0].getStart();
	}

	public AbstractSite getTSS() {
		return new TSSite(this, getTSSPos(), true);
	}
	
	public AbstractSite getTES() {
		return new TSSite(this, getTESPos(), false);
	}
	
	public int getTESPos() {
		return exons[exons.length- 1].getEnd();
	}
	
	public Exon getExon(String stableID) {
		
		if (exons== null)
			return null;
		
		for (int i = 0; i < exons.length; i++) 
			if (exons[i].getExonID().equals(stableID))
				return exons[i];
	
		return null;
	}

	public Exon getLastExon() {
		
		if (exons== null|| exons.length< 1)
			return null;
		return exons[exons.length- 1];
	}

	public SpliceSite getSpliceSite(int pos) {
		
		Comparator compi= new SpliceSite.PositionComparator();
		SpliceSite ss= new SpliceSite(null, pos, true);
		int p= Arrays.binarySearch(spliceChain, ss, compi);
		if (p>= 0)
			return spliceChain[p];
		return null;
	}



}
