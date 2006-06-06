package qalign.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class TimeMSFWrapper extends MSFWrapper {

	protected float hyperTime= -1f, consistentTime= -1f;
	protected boolean consistentTimeFound= false;

	public TimeMSFWrapper(String absFName) {
		
		super(absFName);
	}
	
	public TimeMSFWrapper(String newFName, String newFPath) {
		
		super(newFName, newFPath);
	}	

	public void initMSF() throws Exception {

		BufferedReader fptr = 
			new BufferedReader(new FileReader(fPath+ File.separator+ fName));

		readHeader(fptr);
		readTimes(fptr);	// this is new
		readDescriptor(fptr);
		readSequences(fptr);
	}
	
	public void readTimes(BufferedReader fptr) throws Exception {

		String readStr = "";
		while ((fptr.ready()) && (readStr.indexOf("Times:") == (-1))) {
			if (readStr.trim().length() > 0)
				throw new Exception("Undefined Description Found:\n" + readStr);
			if (readStr.indexOf("Name:") != (-1))
				throw new Exception("No times found:\n" + readStr);
			readStr= fptr.readLine();
		}
		
			// Time line found
		StringTokenizer tokenizer= new StringTokenizer(readStr.trim());
		tokenizer.nextToken();	// Times:
		hyperTime= Float.parseFloat(tokenizer.nextToken());
		if (readStr.indexOf(" + ")!= (-1)) {
			consistentTimeFound= true;
			tokenizer.nextToken(); // (sec)
			tokenizer.nextToken(); // (DCA)
			tokenizer.nextToken(); // ( + )
			consistentTime= Float.parseFloat(tokenizer.nextToken());
		}
		
		
	}	
	/**
	 * Returns the consistentTime.
	 * @return float
	 */
	public float getConsistentTime() {
		return consistentTime;
	}

	/**
	 * Returns the consistentTimeFound.
	 * @return boolean
	 */
	public boolean isConsistentTimeFound() {
		return consistentTimeFound;
	}

	/**
	 * Returns the hyperTime.
	 * @return float
	 */
	public float getHyperTime() {
		return hyperTime;
	}

}
