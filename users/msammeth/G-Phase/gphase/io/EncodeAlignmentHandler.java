/*
 * Created on Nov 29, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.io;

import gphase.Constants;
import gphase.db.EnsemblDBAdaptor;
import gphase.model.EncodeFragment;
import gphase.model.EncodeRegion;
import gphase.model.Species;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Date;
import java.util.StringTokenizer;

import com.sun.corba.se.spi.legacy.connection.GetEndPointInfoAgainException;

/**
 * 
 * 
 * @author msammeth
 */
public class EncodeAlignmentHandler {

	public static final String ENCODE_BASE_DIR= Constants.HOME_DIR+ File.separator+ "encode"; 	//"/projects/encode";
	public static final String ENCODE_MSA_SUBDIR= "data/msa";
	public static final String ENCODE_MSA_SEQ_SUBDIR= "sequences";
	final static String[] months= {
			"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
	};

	public static void main(String[] args) {
		
		EncodeAlignmentHandler testHandler= new EncodeAlignmentHandler();
		Species[] specs= testHandler.getEncodeRegions(Species.createByBinomialName(EnsemblDBAdaptor.SPECIES_ENCODE));
		
		for (int k = 0; k < specs.length; k++) {
			EncodeRegion[] regions= specs[k].getEncodeRegions();
			for (int i = 0; i < regions.length; i++) {
				System.out.println(regions[i]);
				for (int j = 0; j < regions[i].getFragments().length; j++) 
					System.out.println("\t"+ regions[i].getFragments()[j]);
			}
		}
	}
	
	/**
	 * parses American dates of the form "MAY-2005" or "02-may-2004".
	 * @param dateStr
	 * @return
	 */
	public static final Date parseDate(String dateStr) {

		char sepChar= '-';
		int fst= dateStr.indexOf(sepChar);
		int lst= dateStr.lastIndexOf(sepChar);
		
		int day;
		String s_mon;
		if (fst!= lst) {
			day= Integer.parseInt(dateStr.substring(0,fst));
			s_mon= dateStr.substring(fst+1, lst);
		} else
			s_mon= dateStr.substring(0, fst);
		String year= dateStr.substring(lst+1, dateStr.length());
		int mon= 0;
		for(;mon< months.length;mon++)
			if (s_mon.equalsIgnoreCase(months[mon]))
				break;
		
		return new Date(Integer.parseInt(year), mon+1, 1);		
	}
	
	/**
	 * 
	 * @return the canonical path of the newest encode sequences directory
	 */
	public static final String getEncodeSequencesDir(){
		
		File[] dir= new File(ENCODE_BASE_DIR+ File.separator+
				ENCODE_MSA_SUBDIR+ File.separator+ ENCODE_MSA_SEQ_SUBDIR).listFiles();
		String newest= null;
		Date refDate= null;
		for (int i = 0; i < dir.length; i++) {
			if (!dir[i].isDirectory())
				continue;
			
			Date newDate= parseDate(dir[i].getName()); 
			if (refDate== null|| newDate.after(refDate)) {	// find newest
				refDate= newDate;
				try {
					newest= dir[i].getCanonicalPath();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		return newest;

	}

	/**
	 * 
	 * @param species
	 * @return
	 */
	public Species[] getEncodeRegions(Species[] species) {
	
		File[] encodeRegions= new File(getEncodeSequencesDir()).listFiles();
		for (int i = 0; i < encodeRegions.length; i++) {
			
			if ((!encodeRegions[i].getName().startsWith("EN"))
					|| (!encodeRegions[i].isDirectory()))
				continue;
			
			File[] fastaFiles= encodeRegions[i].listFiles();
			for (int j = 0; j < species.length; j++) {
				
				String name= species[j].getCommonName().toLowerCase();		// sometimes they use common names..
				int pos;
				for (pos= 0; pos< fastaFiles.length; ++pos) 
					if (fastaFiles[pos].getName().startsWith(name))
						break;
				if (pos>= fastaFiles.length) {								// .. and sometimes they use first names
					name= species[j].getNameFirst();
					for (pos= 0; pos< fastaFiles.length; ++pos) 
						if (fastaFiles[pos].getName().startsWith(name))
							break;
				}

				if (pos>= fastaFiles.length) {								// some species do not appear in a region
					System.err.println("No alignment region given for "+species[j].getCommonName()+
							" in region "+encodeRegions[i].getName()+ "!");
					continue;
				}
				
				EncodeRegion region= getRegion(species[j], fastaFiles[pos]);	// but sometimes we find something
				species[j].addEncodeRegion(region);
			}
		}
		
		return species;
	}
	
	EncodeRegion getRegion(Species spe, File fastaFile) {
		
		EncodeRegion currRegion= null;
		try {
			RandomAccessFile raf= new RandomAccessFile(fastaFile, "r");
			String line= "";
			int offset= 0;
			while(true) {
				while (!line.startsWith(">"))
					line= raf.readLine().trim();
				EncodeRegion[] transporter= new EncodeRegion[] {currRegion};	// array container for returning value
				int seqlength= parseFastaHeader(line.substring(1), transporter);
				currRegion= transporter[0];
				if (seqlength< 1)
					break;
				
				int lineLength= 50;
				offset+= line.length()+ 1+ seqlength+ (seqlength/ lineLength);
				raf.seek(offset);
				line="";
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NullPointerException e) {
			e.printStackTrace();
		}
		
		return currRegion;
	}
	
	/**
	 * Construtcts parses fasta header and creates new <code>EncodeFragment</code>
	 * or <code>EncodeRegion</code> (iff <code>null</code>) instances.
	 * 
	 * @param line
	 * @param currRegion
	 * @return number of positions to skip to next sequence in multi-fasta file or
	 * (<code>-1</code>) if last sequence has been spotted
	 */
	// niscOrganism    encodeRegion    freezeDate      ncbiTaxonId     assemblyProvider        assemblyDate    assemblyId      
	// chromosome      chromStart      chromEnd        chromLength     strand
	// accessionVersion        numBases        numN    thisContigNum   totalNumContigs comment
	// >chicken|ENm001|SEP-2005|9031|CGSC Feb. 2004|30-MAR-2005|galGal2|chr1|21449967|22183480|188239860|+|NT_107450.4|733514|7422|1|1|.
	int parseFastaHeader(String line, EncodeRegion[] currRegion) {
		
		StringTokenizer toki= new StringTokenizer(line, "|");
		String tmp;
		if (currRegion[0]== null) {
			currRegion[0]= new EncodeRegion(new Species(Species.getBinomialForSomeName(toki.nextToken().trim())));	// mandatory
			currRegion[0].setSilentAndNumber(toki.nextToken().trim());		// mandatory
			try {
				currRegion[0].setFreezeDate(parseDate(toki.nextToken().trim()));		// not madatory
			} catch (StringIndexOutOfBoundsException e) {;}
			try {
				currRegion[0].setNcbiTaxonId(Integer.parseInt(toki.nextToken().trim()));		// not mandatory
			} catch (NumberFormatException e) {;}
			currRegion[0].setAssemblyProvider(toki.nextToken().trim());
			try {
				currRegion[0].setAssemblyDate(parseDate(toki.nextToken().trim()));		// not mandatory
			} catch (StringIndexOutOfBoundsException e) {;}
			
			tmp= toki.nextToken().trim();		// not madatory
			if (!(tmp.equals("")|| tmp.equals(".")))
				currRegion[0].setAssemblyId(tmp);
		} else  // end region init
			for (int i = 0; i < 7; i++) 
				toki.nextToken();
				
		EncodeFragment currFrag= new EncodeFragment(currRegion[0]);
		currRegion[0].addFragment(currFrag);
		tmp= toki.nextToken().trim();		// not madatory
		if (!(tmp.equals("")|| tmp.equals(".")))
			currFrag.setChromosome(tmp);
		try {
			currFrag.setStart(Integer.parseInt(toki.nextToken().trim()));		// not mandatory
		} catch (NumberFormatException e) {;}
		try {
			currFrag.setEnd(Integer.parseInt(toki.nextToken().trim()));		// not mandatory
		} catch (NumberFormatException e) {;}
		try {
			currFrag.setChromLength(Integer.parseInt(toki.nextToken().trim()));		// not mandatory
		} catch (NumberFormatException e) {;}
		currFrag.setPositiveStrand(toki.nextToken().trim().equals("+")?true:false); // always "+", not inited?
		tmp= toki.nextToken().trim();											// not madatory ?
		if (!(tmp.equals("")|| tmp.equals(".")))
			currFrag.setAccessionVersion(tmp);
		try {
			currFrag.setNumBases(Integer.parseInt(toki.nextToken().trim()));		// not mandatory?
		} catch (NumberFormatException e) {;}
		try {
			currFrag.setNumN(Integer.parseInt(toki.nextToken().trim()));		// not mandatory?
		} catch (NumberFormatException e) {;}
		try {
			currFrag.setFragmentID(Integer.parseInt(toki.nextToken().trim()));		// not mandatory?
		} catch (NumberFormatException e) {;}
		
		int maxNb= Integer.parseInt(toki.nextToken().trim());
		if (maxNb> 1) {																// mandatorytotal number of fragments 
			currRegion[0].setStart(-1);											// null region coords		
			currRegion[0].setEnd(-1);
			currRegion[0].setChromosome(null);
		}
		if (maxNb== currFrag.getFragmentID())									// last sequence in multi-fasta file
			return (-1);
		
		if (toki.hasMoreTokens()) {
			tmp= toki.nextToken().trim();											// not madatory 
			if (!(tmp.equals("")|| tmp.equals(".")))
				currFrag.setComment(tmp);
		}
		
		return currFrag.getNumBases();
	}
	
	/**
	 * cross-checks the the sequences used for the Encode regions 
	 * with the substrings extracted from the genomes
	 * @return <code>true</code> if everything's ok
	 */
	boolean checkEncodeRegionSequences() {
		 
		
	}
	
	
	

}
