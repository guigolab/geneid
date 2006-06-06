/*
 * Created on Nov 29, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import java.util.Date;


/**
 * 
 * 
 * @author msammeth
 */
public class EncodeFragment extends DefaultRegion {

	EncodeRegion encodeRegion= null;	// regEnc, "ENm007"
	int chromLength= -1;	    // "4753", redundant but info from flat file
	boolean positiveStrand= false;	// strand, "+"
	String accessionVersion;    // "super_182599"
	int numBases= -1;			// normally same as chromLength, "4753"      
	int numN= -1;    			// "1573"
	int fragmentID= -1;  	// thisContigNum, "23"
	String comment= null;		// "."
	
	public EncodeFragment(EncodeRegion newEncodeRegion) {
		this.encodeRegion= newEncodeRegion;
	}

	
	public String getAccessionVersion() {
		return accessionVersion;
	}
	public void setAccessionVersion(String accessionVersion) {
		this.accessionVersion = accessionVersion;
	}
	public int getChromLength() {
		return chromLength;
	}
	public void setChromLength(int chromLength) {
		this.chromLength = chromLength;
	}
	public String getComment() {
		return comment;
	}
	public void setComment(String comment) {
		this.comment = comment;
	}
	public int getNumBases() {
		return numBases;
	}
	public void setNumBases(int numBases) {
		this.numBases = numBases;
	}
	public int getNumN() {
		return numN;
	}
	
	public Species getSpecies() {
		return getEncodeRegion().getSpecies();
	}
	
	public void getSpecies(Species spec) {
		getEncodeRegion().setSpecies(spec);
	}	
	public void setNumN(int numN) {
		this.numN = numN;
	}
	public EncodeRegion getEncodeRegion() {
		return encodeRegion;
	}
	public void setEncodeRegion(EncodeRegion newEncodeRegion) {
		this.encodeRegion = newEncodeRegion;
	}
	public int getFragmentID() {
		return fragmentID;
	}
	public void setFragmentID(int newFragmentID) {
		this.fragmentID = newFragmentID;
	}
	public boolean isPositiveStrand() {
		return positiveStrand;
	}
	public void setPositiveStrand(boolean positiveStrand) {
		this.positiveStrand = positiveStrand;
	}
	
	public String toString() {
		return getEncodeRegion().getName()+ ","+ getFragmentID()+ ":"+ super.toString();
	}
}
