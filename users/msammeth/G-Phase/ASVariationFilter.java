import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

/*
 * Created on Mar 31, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

/**
 * 
 * 
 * @author msammeth
 */
public class ASVariationFilter {

	String fileAbsPath= null;
	String outFileAbsPath= null;
	
	public static void main(String[] args) {
		ASVariationFilter myFilter= new ASVariationFilter("results\\encode_beta.txt");
		myFilter.filter();
	}
	
	public ASVariationFilter(String newFileAbsPath) {
		this.fileAbsPath= newFileAbsPath;
		this.outFileAbsPath= newFileAbsPath+ "_filtered";
	}

	public void filter() {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fileAbsPath));
			BufferedWriter buffy2= new BufferedWriter(new FileWriter(outFileAbsPath));
			
			while (buffy.ready()) {
				String line= buffy.readLine();
				line= filterLength(line);
				if (line!= null)
					buffy2.write(line+"\n");
			}
			
			buffy2.flush();
			buffy2.close();
		} catch (Exception e) {
			System.err.println(e);
		}
	}
	
	String filterLength(String line) {
		if (line.length()> 50)
			return null;
		return line;
	}
}
