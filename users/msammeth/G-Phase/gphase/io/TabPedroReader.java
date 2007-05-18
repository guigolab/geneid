package gphase.io;

import gphase.tools.Arrays;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.Vector;

public class TabPedroReader extends TabDelimitedFormatWrapper {

	long bytesRead= 0l;
	
	public TabPedroReader(String absFName) {
		super(absFName);
	}

	@Override
	public void read() throws Exception {
		Vector v= new Vector();
		try {
			File f= new File(fPath+ File.separator+ fName);
			String lf= FileTools.determineLineFeed(f);
			InputStream stream= new FileInputStream(f);
			stream.skip(bytesRead);
			BufferedReader buffy= new BufferedReader(new InputStreamReader(stream));
			int fCnt= -1;
			String chrTag= null;
			while (buffy.ready()) {
				String line= buffy.readLine();
				if (line.trim().length()< 1)
					continue;				
				String[] tokens= line.split("\t");
				if (chrTag== null)
					chrTag= tokens[0];
				else if (!chrTag.equals(tokens[0])) 
					break;
				
				bytesRead+= line.length()+ lf.length();
				if (fCnt< 0)
					fCnt= tokens.length;
				else
					if (fCnt!= tokens.length) {
						System.err.println("Invalid line, expected "+fCnt+" tokens, but found "
								+tokens.length+"\n"+line);
						continue;
					}
				
				v.add(tokens);
			}
			buffy.close();			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		table= (String[][]) Arrays.toField(v);
	}

	public void sweepToChromosome(String chrom) {
		try {
			File f= new File(fPath+ File.separator+ fName);
			String lf= FileTools.determineLineFeed(f);
			InputStream stream= new FileInputStream(f);
			stream.skip(bytesRead);
			BufferedReader buffy= new BufferedReader(new InputStreamReader(stream));
			int fCnt= -1;
			boolean startOver= true;
			if (bytesRead== 0l)
				startOver= false;
			while (true) {
				String line= buffy.readLine();
				if (line== null) {
					if (startOver) {
						buffy= new BufferedReader(new InputStreamReader(stream));
						bytesRead= 0l;
						continue;
					} else
						break;
				}
				String[] tokens= line.split("\t");
				if (tokens== null|| tokens.length< 1)
					continue;
				if (chrom.equals(tokens[0])) 
					break;
				
				bytesRead+= line.length()+ lf.length();
			}
			buffy.close();			
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
}
