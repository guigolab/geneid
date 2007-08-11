package gphase.io.bed;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;

import gphase.io.DefaultIOWrapper;
import gphase.tools.Arrays;
import gphase.tools.File;

public class BEDwrapper extends DefaultIOWrapper {

	BEDobject[] beds= null;
	
	public BEDwrapper(String newFilePath) {
		super(newFilePath);
	}
	
	public boolean isApplicable() throws Exception {
		// TODO Auto-generated method stub
		return false;
	}

	public void read() throws Exception {

		Vector objV= new Vector();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(this.fPath+File.separator+this.fName));
			while (buffy.ready()) {
				String line= buffy.readLine();
				if (line== null)
					break;
				String[] tokens= line.split("\\s");
				if (tokens.length< 3) {
					System.out.println("WARNING: skipped incomplete line with "+tokens.length+" token");
					continue;
				}
				
				BEDobject bed= new BEDobject();
				bed.setChrom(tokens[0]);
				bed.setStart(Integer.parseInt(tokens[1]));
				bed.setEnd(Integer.parseInt(tokens[2]));
				if (tokens.length> 3) 
					bed.setName(tokens[3]);
				if (tokens.length> 4) 
					bed.setScore(Integer.parseInt(tokens[4]));
				
				if (tokens.length> 5) 
					bed.setStrand(tokens[5]);
				if (tokens.length> 6) 
					bed.setThickStart(Integer.parseInt(tokens[6]));
				if (tokens.length> 7) 
					bed.setThickEnd(Integer.parseInt(tokens[7]));
					
				if (tokens.length> 8) 
					bed.setCol(tokens[8]);

				if (tokens.length> 9) 
					bed.setBlockCount(Integer.parseInt(tokens[9]));
				if (tokens.length> 10) 
					bed.setBlockSizes(tokens[10]);
				if (tokens.length> 11)  
					bed.setBlockStarts(tokens[11]);

				objV.add(bed);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		

		beds= (BEDobject[]) Arrays.toField(objV);
	}

	public void write() throws Exception {
		write(false);
	}
	
	public void write(boolean append) {
		try {
			BufferedWriter buffy= new BufferedWriter(new FileWriter(this.fPath+File.separator+this.fName, append));
			for (int i = 0; beds!= null&& i < beds.length; i++) {
				buffy.write(beds[i].toString()+"\n");
			}
			buffy.flush();
			buffy.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	public BEDobject[] getBeds() {
		return beds;
	}

	public void setBeds(BEDobject[] beds) {
		this.beds = beds;
	}

}
