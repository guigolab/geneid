package gphase.tools;

import java.io.IOException;
import java.io.PrintStream;

public class File extends java.io.File {
	
	public static boolean checkForOverwrite(PrintStream p, File f) {
		if (!f.exists())
			return true;
		p.println("Confirm overwriting file "+f+" (y/n)");
		int b= 'n';
		try {
			b= System.in.read();
		} catch (IOException e) {
			e.printStackTrace();
		}
		if (b== 'y'|| b== 'Y') {
			f.delete();
			return true;
		}
		return false;
	}
	
	public File(String name) {
		super(name);
	}
	
	public String getPathOnly() {
		int pos= getAbsolutePath().lastIndexOf(File.separator);
		if (pos< 0)
			return getAbsolutePath();
		return getAbsolutePath().substring(0,pos);		
	}
	public String getFileNameOnly() {
		int pos= getAbsolutePath().lastIndexOf(File.separator);
		if (pos< 0)
			return getAbsolutePath();
		return getAbsolutePath().substring(pos+1);
	}

	public String getExtension() {
		int pos= getFileNameOnly().lastIndexOf('.');
		if (pos< 0)
			return null;
		return getFileNameOnly().substring(pos+1);
	}
	public String getFileNameWithoutExtension() {
		int pos= getFileNameOnly().lastIndexOf('.');
		if (pos< 0)
			return null;
		return getFileNameOnly().substring(0, pos);
	}

}
