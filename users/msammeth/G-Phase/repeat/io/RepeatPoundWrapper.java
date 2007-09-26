/*
 * Created on Oct 14, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import repeat.data.Repeat;

/**
 * 
 * 
 * @author micha
 */
public class RepeatPoundWrapper {
	
	public static final String MARKER_POUND= "#";
	public static final String MARKER_TOR= "TYPE_OF_REPEATS";
	public static final String MARKER_EOTOR= "END_OF_REPEAT_TYPE";
	public static final String MARKER_SUBTYPES= "SUBTYPES";
	public static final String MARKER_LIST= "LIST_OF_REPEATS";
	public static final String MARKER_ALI= "ALIGNMENT";
	
	protected String fileBase= null;
	protected Repeat[] repeats= null;
	
	public static void main(String[] args) {
	}
	
	public RepeatPoundWrapper(String newFileBase) {
		this.fileBase= newFileBase;		
	}
	
	public RepeatPoundWrapper(String newFileBase, Repeat[] newRepeats) {
		this(newFileBase);
		this.repeats= newRepeats;		
	}
	
	protected Repeat[] read() {
		
			// create reader
		BufferedReader buffy= null;
		try {
			buffy= new BufferedReader(new FileReader(fileBase));
		} catch (FileNotFoundException e) {System.err.println(e);;}	// :)
		
			// read in
		Vector repVec= new Vector();
		try {
			while (buffy.ready()) {
				Vector subVec= readRepeatType(buffy);
				for (int i= 0; (subVec!= null)&& i < subVec.size(); i++) 
					repVec.add(subVec.elementAt(i));
			}
		} catch (IOException e) {;}	// :)

		
			// convert
		Repeat[] rep= new Repeat[repVec.size()];
		for (int i= 0; i < rep.length; i++) 
			rep[i]= (Repeat) repVec.elementAt(i);
			
		return rep;
	}
	
	protected Vector readRepeatType(BufferedReader buffy) throws IOException {
		
			// tor
		String marker= MARKER_POUND+" "+MARKER_TOR;
		String line= buffy.readLine();
		if (line== null)
			return null;
		if (line.trim().length()== 0) {
			line= buffy.readLine();
			if (line== null)
				return null; 
		}
		while(line.indexOf(marker)< 0)
			line= buffy.readLine();
		String type= buffy.readLine().trim();
		
			// list 
		marker= MARKER_POUND+" "+MARKER_LIST;
		line= buffy.readLine();
		while(line.indexOf(marker)< 0)
			line= buffy.readLine();
		line= buffy.readLine();		// header
		line= buffy.readLine();
		Pattern patty= Pattern.compile("^(\\S+)\\s+.(\\d+).\\s+(\\d+)\\s+(\\d+)\\s+(\\w+).*");	
		Matcher matty= patty.matcher(line);
		Vector repVec= new Vector();
		while (matty.matches()) {
			
			Repeat tmpRep= new Repeat();
			tmpRep.setSeqName(matty.group(1));
			tmpRep.setNb(Integer.parseInt(matty.group(2)));
			tmpRep.setStart(Integer.parseInt(matty.group(3)));
			tmpRep.setLength(Integer.parseInt(matty.group(4)));
			tmpRep.setType(matty.group(5));
			repVec.add(tmpRep);

			line= buffy.readLine();			// next line
			matty= patty.matcher(line);
		}
		
			// check for alignment
		while (buffy.ready()) {
			if (line.indexOf(MARKER_EOTOR)> 0)		// end-of-type
				break;
			if (line.indexOf(MARKER_ALI)> 0) {		// alignment
				line= buffy.readLine();
				patty= Pattern.compile("^(\\S+)\\s+.(\\d+).\\s+(.+)\\w*");
				matty= patty.matcher(line);
				while(!matty.matches()) {
					line= buffy.readLine();
					matty= patty.matcher(line);
				}
				while(matty.matches()) {
					String name= matty.group(1);
					int nb= Integer.parseInt(matty.group(2));
					for (int i= 0; i < repVec.size(); i++) {
						Repeat tmpRep= (Repeat) repVec.elementAt(i);
						if (name.equals(tmpRep.getSeqName())
							&& nb == tmpRep.getNb()) {
								tmpRep.setAlignedSeq(matty.group(3));
								tmpRep.setAlignmentID(type);
							}
					}
					line= buffy.readLine();
					matty= patty.matcher(line);
				}
				continue;
			}
			line= buffy.readLine();
		}
		
		return repVec;
	}
	
	public Repeat[] getRepeats() {
		if (repeats == null) {
			repeats= read();
		}
		return repeats;
	}
	
	public void write() {
		
			// open writer
		BufferedWriter buffy= null;
		try {
			buffy= new BufferedWriter(new FileWriter(fileBase));
		} catch (IOException e) {
			; // :)
		}
		
			// count types
		Vector collectedTypes= new Vector();
		for (int i= 0; i < repeats.length; i++) {
			String id= repeats[i].getAlignmentID();
			int j;
			for (j= 0; j < collectedTypes.size(); j++) {
				if (id.equals(collectedTypes.elementAt(j)))
					break;
			}
			if (j>= collectedTypes.size())
				collectedTypes.insertElementAt(id, 0);		// insert at front for efficient comparison
		}
		
			// group and submit for writing
		for (int i= collectedTypes.size()- 1; i>= 0; --i) {
			Vector filteredRepeats= new Vector();			// filter repeats
			for (int j= 0; j < repeats.length; j++) 
				if (repeats[j].getAlignmentID().equals((String) collectedTypes.elementAt(i)))
					filteredRepeats.add(repeats[j]);
			
			Repeat[] reps= new Repeat[filteredRepeats.size()];	// convert & sort
			for (int j= 0; j < reps.length; j++) 
				reps[j]= (Repeat) filteredRepeats.elementAt(j);
			Arrays.sort(reps, new Comparator() {
				public boolean equals(Object obj) {
					
					if (compare(this, obj)== 0)
						return true;
					else
						return false;
				}

				public int compare(Object o1, Object o2) {
					
					Repeat r1= (Repeat) o1;
					Repeat r2= (Repeat) o2;
					
						// first: seqName
					if (!r1.getSeqName().equals(r2.getSeqName()))
						return r1.getSeqName().compareTo(r1.getSeqName());
					
						// else.. start
					if (r1.getStart()!= r2.getStart())
						if (r1.getStart()< r2.getStart())
							return (-1);
						else
							return 1;
					
						// (else.. end)
					if (r1.getLength()!= r2.getLength())
						if (r1.getLength()< r2.getLength())
							return (-1);
						else
							return 1;
					
					return 0;// (else subtype?)
				}
			});
			
				// submit to writing
			try {
				writeRepeatType(buffy, (String) collectedTypes.elementAt(i), reps);
			} catch (IOException e) {
				; // :)
			}
		}
		
			// close writer
		try {
			buffy.flush();
			buffy.close();
		} catch (IOException e) {
			; // :)
		}
	}
	
	protected void writeRepeatType(BufferedWriter buffy, String typeID, Repeat[] rep) throws IOException {
		
		buffy.write(MARKER_POUND+" "+MARKER_TOR+"\n"+typeID+"\n\n");	// header
		
			//count and write subtypes
		int longestName= 0;
		Vector subtypeIDs= new Vector();
		for (int i= 0; i < rep.length; i++) {
			if (rep[i].getSeqName().length()> longestName)
				longestName= rep[i].getSeqName().length();	// remember longest Name
			Object id= rep[i].getType();
			int j;
			for (j= 0; j < subtypeIDs.size(); j++) {
				if (id.equals(subtypeIDs.elementAt(j)))
					break;
			}
			if (j>= subtypeIDs.size())
				subtypeIDs.insertElementAt(id, 0);		// insert at front for efficient comparison
		}
		buffy.write(MARKER_POUND+" "+MARKER_SUBTYPES+"\n");
		for (int i= subtypeIDs.size()- 1; i >= 0; --i) {
			buffy.write(subtypeIDs.elementAt(i)+"\n");
		}
		buffy.write("\n");
		
			// write list header
		buffy.write(MARKER_POUND+" "+MARKER_LIST+"\n");
		String stretcher="";
		int insertTabs= (longestName)/ 8;			// "Start" = 5 chars, (xy)= 5 chars: 0 total
		for (int j= 0; j < insertTabs; j++) {
			stretcher+= "\t";
		}
		buffy.write("Name"+stretcher+"Start\tLength\tSubtype\t[Score]\tSeqID_Nb\n");
		Hashtable hash= new Hashtable();
		for (int i= 0; i < rep.length; i++) {
			
			int ctr= 0;		// counter for repeats in same seq
			if (hash.get(rep[i].getSeqName())== null) 
				ctr= 1;
			else
				ctr= ((Integer) hash.get(rep[i].getSeqName())).intValue();
			hash.put(rep[i].getSeqName(), new Integer(ctr+1));
			
			stretcher="";
			insertTabs= (longestName+5)/8- 
				(rep[i].getSeqName().length()+3+Integer.toString(ctr).length())/8;	
			for (int j= 0; j < insertTabs; j++) {
				stretcher+= "\t";
			}
			buffy.write(
				rep[i].getSeqName()+ " ("+ctr+")\t"+
				stretcher+
				rep[i].getStart()+ "\t"+
				rep[i].getLength()+ "\t"+
				rep[i].getType()+ "\t"+		// Subtype
				"-\n"						// Score
			);
		}
		buffy.write("\n");

			// write alignment		
		buffy.write(MARKER_POUND+" "+MARKER_ALI+"\n");
		hash= new Hashtable();
		for (int i= 0; i < rep.length; i++) {
			
			int ctr= 0;		// counter for repeats in same seq
			if (hash.get(rep[i].getSeqName())== null) 
				ctr= 1;
			else
				ctr= ((Integer) hash.get(rep[i].getSeqName())).intValue();
			hash.put(rep[i].getSeqName(), new Integer(ctr+1));
			
			stretcher="";			
			insertTabs= (longestName+5)/8- 
				(rep[i].getSeqName().length()+3+Integer.toString(ctr).length())/8;	
			for (int j= 0; j < insertTabs; j++) {
				stretcher+= "\t";
			}
			buffy.write(
				rep[i].getSeqName()+ " ("+ctr+")\t"+
				stretcher+
				rep[i].getAlignedSeq()+ "\n"
			);
		}
		buffy.write("\n");
		buffy.write(MARKER_POUND+" "+MARKER_EOTOR+"\n\n");
		
	}
	/**
	 * @param repeats
	 */
	public void setRepeats(Repeat[] repeats) {
		this.repeats= repeats;
	}

}
