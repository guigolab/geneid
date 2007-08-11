import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

import javax.annotation.processing.Filer;

public class FASTAsplit {
	static final String usage= "Hi!\n\nI will split a multi-FASTA file if you tell me:\nFASTAsplit inputFile [outputDir]\n\nHave fun!";

	String fName= null, outDir= null;

	public static void main(String[] args) {
		if (args.length< 1)  { 
			System.err.println(usage);
			System.exit(-1);
		}

		File f= new File(args[0]);
		if (!f.exists()) {
			System.err.println("File not found "+f.getAbsolutePath());
			System.exit(-1);
		}		
		String s= f.getAbsolutePath();
		FASTAsplit splitter= new FASTAsplit(s);
		int p= s.lastIndexOf(File.separator);
		if (p>= 0)
			s= s.substring(0, p);
		splitter.setOutDir(s);
		
		if (args.length> 1) {
			f= new File(args[1]);
			if ((!f.exists())|| (!f.isDirectory())) {
				System.err.println("Invalid output dir "+f.getAbsolutePath());
				System.exit(-1);
			}
			splitter.setOutDir(f.getAbsolutePath());
		}
		
		splitter.run();
	}
	
	public FASTAsplit(String newFName) {
		this.fName= newFName;
	}

	public String getOutDir() {
		return outDir;
	}

	public void setOutDir(String outDir) {
		this.outDir = outDir;
	}
	
	public void run() {
		try {
			BufferedReader reader= new BufferedReader(new FileReader(fName));
			BufferedWriter writer= null;
			while (reader.ready()) {
				String line= reader.readLine();
				if (line.charAt(0)== '>') {
					if (writer!= null) {
						writer.flush();
						writer.close();
					}
					String[] tokens= line.substring(1).split("\\s");					
					writer= new BufferedWriter(new FileWriter(outDir+File.separator+tokens[0]+".fa"));					
				}
				writer.write(line+"\n");
			}
			if (writer!= null) {
				writer.flush();
				writer.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
