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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Collection;
import java.util.Iterator;
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
	String[] sortAttributes= null;
	
	
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
	
	public GTFWrapper() {
	}
	
	/* (non-Javadoc)
	 * @see gphase.io.IOWrapper#isApplicable()
	 */
	public boolean isApplicable() throws Exception {
		// TODO Auto-generated method stub
		return false;
	}
	
	public void addFileSuffix(String newSfx) {
		if (fName== null)
			return;
		this.fName+= newSfx;
	}
	
	/**
	 */
	public void write() throws Exception {
		writeGTF();
	}
	
	public void writeGTF() throws Exception {
			String outName= getAbsFileName();
	//		if (new File(outName).exists())
	//			outName+= "_out";
			BufferedWriter buffy= new BufferedWriter(new FileWriter(outName));
			for (int i = 0; i < gtfObj.length; i++) {
				// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
				buffy.write(
						gtfObj[i].seqname+ "\t"+	//!! getter removes "chr"
						gtfObj[i].getSource()+ "\t"+
						gtfObj[i].getFeature()+ "\t"+
						gtfObj[i].getStart()+ "\t"+
						gtfObj[i].getEnd()+ "\t"+
						gtfObj[i].getScore()+ "\t"+
						gtfObj[i].getStrand()+ "\t"+
						gtfObj[i].getFrame()+ "\t"
				);
				Collection c= gtfObj[i].getAttributes().keySet();
				Iterator iter= c.iterator();
				String[] keys= new String[c.size()];
				int cc= 0;
				while (iter.hasNext())
					keys[cc++]= (String) iter.next();
				
				if (sortAttributes!= null) 
					for (int j = 0; j < sortAttributes.length; j++) 
						for (int k = (j+1); k < keys.length; k++) 
							if (sortAttributes[j].equals(keys[k])) {
								String h= keys[k];
								keys[k]= keys[j];
								keys[j]= h;
								break;
							}
				
				for (int j = 0; j < keys.length; j++) {
					buffy.write(keys[j]+ " \""+ gtfObj[i].getAttribute(keys[j])+ "\"; ");
				}
				buffy.write("\n");
			}
			buffy.flush(); buffy.close();
		}
	public void writeGFF() throws Exception {
		String outName= getAbsFileName();
		BufferedWriter buffy= new BufferedWriter(new FileWriter(outName));
		
		for (int i = 0; gtfObj!= null&& i < gtfObj.length; i++) {
           buffy.write(
					gtfObj[i].seqname+ "\t"+	//!! getter removes "chr"
					gtfObj[i].getSource()+ "\t"+
					gtfObj[i].getFeature()+ "\t"+
					gtfObj[i].getStart()+ "\t"+
					gtfObj[i].getEnd()+ "\t"+
					gtfObj[i].getScore()+ "\t"+
					gtfObj[i].getStrand()+ "\t"+
					gtfObj[i].getFrame()+ "\t"
			);
			Collection c= gtfObj[i].getAttributes().keySet();
			Iterator iter= c.iterator();
			String[] keys= new String[c.size()];
			int cc= 0;
			while (iter.hasNext())
				keys[cc++]= (String) iter.next();
			
			java.util.Arrays.sort(keys);
//			if (sortAttributes!= null) 
//				for (int j = 0; j < sortAttributes.length; j++) 
//					for (int k = (j+1); k < keys.length; k++) 
//						if (sortAttributes[j].equals(keys[k])) {
//							String h= keys[k];
//							keys[k]= keys[j];
//							keys[j]= h;
//							break;
//						}
			
			for (int j = 0; j < keys.length; j++) {
				buffy.write(gtfObj[i].getAttribute(keys[j])+ "\t");
			}
			buffy.write("\n");
		}
		buffy.flush(); buffy.close();
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
			StringTokenizer toki= new StringTokenizer(line, " \t");	// must be tab, see specification
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
				newObj.setStrand(toki.nextToken());
				newObj.setFrame(toki.nextToken());
			} catch (Exception e) {
				System.err.println("*line "+ lineCtr+": "+ e);
				e.printStackTrace();
				//continue;
			}
			
				// optional attributes
			int smc= line.indexOf(';');		// GTF2
			if (smc>= 0) {
				String ss= toki.nextToken();
//				toki= new StringTokenizer(line, " \t");	// must be tab, see specification
//				for (int i = 0; i < 8; i++) 
//					ss= toki.nextToken();
//				String h= line.substring(0, smc);			// last ';'
//				h= line.substring(0, h.lastIndexOf(' '));	// two ' ' tokens before
//				h= line.substring(0, h.lastIndexOf(' '));
				String h= line.substring(line.indexOf(ss), line.length()).trim();	// skip that part
				
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
	public void setGtfObj(GTFObject[] gtfObj) {
		this.gtfObj = gtfObj;
	}
	public void setSortAttributes(String[] sortAttributes) {
		this.sortAttributes = sortAttributes;
	}
}
