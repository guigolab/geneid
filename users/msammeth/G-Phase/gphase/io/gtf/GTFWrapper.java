/*
 * Created on Mar 2, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package gphase.io.gtf;

import gphase.io.DefaultIOWrapper;
import gphase.tools.Arrays;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * see http://genes.cs.wustl.edu/GTF2.html for description
 * @author micha
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class GTFWrapper extends DefaultIOWrapper {
	
	GTFObject[] gtfObj= null;
	
	public static void main(String[] args) {
		GTFWrapper myWrapper= new GTFWrapper(new File("encode/44regions_genes_CHR_coord.gtf").getAbsolutePath()); // testGTF.gtf
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		for (int i = 0; i < myWrapper.getGtfObj().length; i++) 	
			System.out.println(myWrapper.getGtfObj()[i]);
		
	}
	public GTFWrapper(String absFName) {
		super(absFName);
	}
	
	/* (non-Javadoc)
	 * @see gphase.io.IOWrapper#isApplicable()
	 */
	public boolean isApplicable() throws Exception {
		// TODO Auto-generated method stub
		return false;
	}
	
	/**
	 * @deprecated
	 */
	public void write() throws Exception {
		// TODO Auto-generated method stub

	}
	

	GTFObject createGTFObject(){
		return new GTFObject();
	}
	
	
	/* (non-Javadoc)
	 * @see gphase.io.IOWrapper#read()
	 */
	public void read() throws Exception {
		
		BufferedReader buffy= new BufferedReader(new FileReader(fPath+ File.separator+ fName));
		String line;
		int lineCtr= 0;
		Vector gtfVec= new Vector();
		while (buffy.ready()) {
			lineCtr++;
			line= buffy.readLine();
			StringTokenizer toki= new StringTokenizer(line, "\t");	// must be tab, see specification
			if (toki.countTokens()< 8)
				System.err.println("line "+ lineCtr+ ": skipped (<8 token)!\n\t"+ line);
			// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
			GTFObject newObj= createGTFObject();
			try {				
				newObj.seqname= toki.nextToken();
				newObj.source= toki.nextToken();
				newObj.setFeature(toki.nextToken());
				newObj.start= Integer.parseInt(toki.nextToken());
				newObj.end= Integer.parseInt(toki.nextToken());
				newObj.setScore(toki.nextToken());
				newObj.setLeadingStrand(toki.nextToken());
				newObj.setFrame(toki.nextToken());
			} catch (Exception e) {
				System.err.println("*line "+ lineCtr+": "+ e);
				e.printStackTrace();
				//continue;
			}
			
				// optional attributes
			int smc= line.indexOf(';');		// GTF2
			if (smc>= 0) {
				int x= 0;
				toki= new StringTokenizer(line, "\t");	// must be tab, see specification
				for (int i = 0; i < 8; i++) 
					x+= toki.nextToken().length()+ 1;
//				String h= line.substring(0, smc);			// last ';'
//				h= line.substring(0, h.lastIndexOf(' '));	// two ' ' tokens before
//				h= line.substring(0, h.lastIndexOf(' '));
				String h= line.substring(x, line.length()).trim();	// skip that part
				
				toki= new StringTokenizer(h, ";");		// attributes
				while (toki.hasMoreTokens()) {
					h= toki.nextToken().trim();
					int sep= h.indexOf(' ');
					if (sep < 0) {						// comments
						String s= h;
						while (toki.hasMoreTokens())
							s+= " "+ toki.nextToken();
						newObj.setComments(s);
					}
					String id= h.substring(0, sep);
					String val= h.substring(sep+ 1, h.length());
					newObj.addAttribute(id, val);
				}
			}
				
			gtfVec.add(newObj);
		}
		
		gtfObj= (GTFObject[]) Arrays.toField(gtfVec);
	}
	/**
	 * @return Returns the gtfObj.
	 */
	public GTFObject[] getGtfObj() {
		return gtfObj;
	}
}
