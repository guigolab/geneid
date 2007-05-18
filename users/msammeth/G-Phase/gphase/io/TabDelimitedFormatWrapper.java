package gphase.io;

import gphase.tools.Arrays;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.util.Vector;

public class TabDelimitedFormatWrapper extends DefaultIOWrapper {

	String[][] table= null;

	public TabDelimitedFormatWrapper(String absFName) {
		super(absFName);
	}
	
	public void write() throws Exception {
		System.err.println("implement..");
	}

	public boolean isApplicable() throws Exception {
		System.err.println("implement..");
		return false;
	}
	
	public void read() throws Exception {
		Vector v= new Vector();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fPath+ File.separator+ fName));
			int fCnt= -1;
			while (buffy.ready()) {
				String line= buffy.readLine();
				if (line.trim().length()< 1)
					continue;
				StringTokenizer toki= new StringTokenizer(line, "\t");
				if (toki.countTokens()> fCnt)
					fCnt= toki.countTokens();
//				else
//					if (fCnt!= toki.countTokens()) {
//						System.err.println("Invalid line, expected "+fCnt+" tokens, but found "
//								+toki.countTokens()+"\n"+line);
//						continue;
//					}
				String[] ln= new String[fCnt];
				int len= toki.countTokens();
				for (int i = 0; i < len; i++) 
					ln[i]= toki.nextToken();
				for (int i = len; i < fCnt; i++) 
					ln[i]= "MISS";
				v.add(ln);
			}
			buffy.close();			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		table= (String[][]) Arrays.toField(v);
	}

	public String[][] getTable() {
		return table;
	}
	
	public String[] getColumn(int colNr) {
		if (table== null)
			return null;
		String[] col= new String[table.length];
		for (int i = 0; i < col.length; i++) 
			col[i]= table[i][colNr];
		return col;
	}

	public void setTable(String[][] table) {
		this.table = table;
	}
}
