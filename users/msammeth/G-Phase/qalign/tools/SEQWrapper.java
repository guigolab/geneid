package qalign.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.sql.Date;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class SEQWrapper extends DefaultIOWrapper {

	protected String[] sequences= null;
	protected String[] seqNames= null;
	protected String[] seqComments= null;

	public SEQWrapper(String newFName, String newFPath) {
		
		super(newFName, newFPath);
	}
	
	public SEQWrapper(String absFName) {
		
		super(absFName);
	}
	
	/**
	 * @throws all exceptions during file handling to keep <code>this</code> independent from UI
	 */
	public void read() throws Exception {
		BufferedReader fptr = null;
//		try {
			fptr = new BufferedReader(new FileReader(fPath+ File.separator+ fName));
			
			String read= "";
			sequences= new String[1];
			seqComments= new String[1];
			seqNames= new String[1];
			seqNames[0]= "Unknown Sequence";
			sequences[0]= "";
			seqComments[0]= "";

			while ((fptr.ready())&& (!read.startsWith("^^"))) {
				seqComments[0]+= read;
				read= fptr.readLine().trim();
			}
			
			while (fptr.ready())
				sequences[0]+= fptr.readLine().trim();

/*		} catch (Exception e) {
			e.printStackTrace();
			// do nothing
		}
*/		
	}
	
	/**
	 * @throws all exceptions during file handling to keep <code>this</code> independent from UI
	 */
	public void write() throws Exception {
		
		for (int i= 0; i< sequences.length; ++i) {
			
				// create file
			String modName= fPath+ File.separator;
			if (i> 0)
				modName+= i+ fName;
			else
				modName+= fName;
			FileWriter writer= new FileWriter(modName);
			
				// write name
			writer.write(seqNames[i]);
			
				// write comment			
			if ((seqComments!= null)|| (seqComments.length> 0))
				writer.write(seqComments[i]+"\n");
			else
				writer.write(new Date(System.currentTimeMillis()).toString()+"\n");
	
				// write separator
			writer.write("^^\n");
		
				// write sequence
			writer.write(sequences[i]+"\n");		
			writer.close();
		}
			
	}
	
	public SequenceWrapper[] getWrappedSequences() {
		
		SequenceWrapper[] result= new SequenceWrapper[sequences.length];
		for (int i= 0; i< result.length; ++i) {
			result[i]= new SequenceWrapper(
				seqNames[i], 
				seqComments[i],
				fName,
				getFDate()
			);
			result[i].setSequence(sequences[i]);
			result[i].setWrapper(this);
		}
		
		return result;
	}
	
	public void setWrappedSequences(SequenceWrapper[] newData) {
		
		seqComments= new String[newData.length];
		seqNames= new String[newData.length];
		sequences= new String[newData.length];
		for(int i= 0; i< newData.length; ++i) {
			
			seqComments[i]= newData[i].getComment();
			seqNames[i]= newData[i].getName();
			sequences[i]= newData[i].getSequence();
		}
		
	}	
}
