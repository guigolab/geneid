package gphase.io;

import java.io.File;
import java.util.Date;



/**
 *  Default implementation of the IOWrapper Interface
 * 
 * @author micha
 */
public abstract class DefaultIOWrapper implements IOWrapper {

	/**
	 * File name and extension.
	 */
	protected String fName= null;
	
	/**
	 * File path without seperator.
	 */
	protected String fPath= null;
	
	/**
	 * File last modification
	 */
	protected long fLastModified= 0L;
	
	public DefaultIOWrapper(String newFName, String newFPath) {	
		this.fName= newFName;
		this.fPath= newFPath;
	}
	
	public DefaultIOWrapper(String absFilePath) {		
		int p= absFilePath.length();
		while ((--p>= 0)&& (absFilePath.charAt(p)!= File.separatorChar))
			; // count down
		this.fPath= absFilePath.substring(0,p);
		this.fName= absFilePath.substring((p+1), absFilePath.length());		
	}
	
	
	/**
	 * @see qalign.tools.IOWrapper#getFileName()
	 */
	public String getFileName() {
		return fName;
	}
	
	public String getAbsFileName() {		
		if (fPath!= null)
			return getFilePath()+ File.separator+ getFileName();			
		return getFileName();
	}
	
	public String getFilePath() {		
		if (fPath!= null)
			return fPath;			
		return null;
	}
	
	
	public long getfLastModified() {		
		if (fLastModified == 0L) {
			File file= new File(fPath+ File.separator+ fName);
			if (file!= null)
				fLastModified = file.lastModified();
		}
		return fLastModified;
	}

	public Date getFDate() {		
		if (getfLastModified()== 0L)
			return null;		
		Date date= new Date(getfLastModified());
		return date;
	}
	
	public static int[] extentIntArray(int[] old) {		
		int[] extended= new int[old.length+ 1];
		for (int i= 0; i< old.length; ++i)
			extended[i]= old[i];
		return extended;
	}
	
	public static float[] extentFloatArray(float[] old) {		
		float[] extended= new float[old.length+ 1];
		for (int i= 0; i< old.length; ++i)
			extended[i]= old[i];
		return extended;
	}	
	
	public static String[] extentStringArray(String[] old) {
		if (old== null) {
			String[] extended= new String[0];
			old= extended;
			return extended;
		} else {
			String[] extended= new String[old.length+ 1];
			for (int i= 0; i< old.length; ++i)
				extended[i]= old[i];
			old= extended;
			return extended;
		}
	}	
}
