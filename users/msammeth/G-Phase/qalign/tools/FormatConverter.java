package qalign.tools;

import java.io.File;
import java.util.StringTokenizer;

/**
 * @author sammeth
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class FormatConverter {
	
	public static String[] extensions= {
		"MSF",
		"MS"
	};
	protected String fileName= null;
	protected int inFormat= 0;
	protected int outFormat= 0;

	public static void main(String[] args) {
	}
	
	public FormatConverter(String fname, String outExtension) {
		
		setFileName(fname);
		setInputFormat(fname);
		setOutputFormat(outExtension);
	}
	
	public void setInputFormat(String fname) {
		
		String ext= getExtension(fname);
		ext= ext.toUpperCase();
		
		String[] validExts= null;
			// check for MSF
		
		if (ext.equals("")
	}
	
	public String getExtension(String fname) {
		
		StringTokenizer st= new StringTokenizer(fname, ".");
		
		String token= null;
		while (st.hasMoreElements())	// get last Token
			token= st.nextToken();
		
		return token;
	}
	
}
