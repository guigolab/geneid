package gphase.tools;

import com.sun.java_cup.internal.internal_error;
import com.sun.org.apache.bcel.internal.generic.FNEG;

public class File extends java.io.File {
	public File(String name) {
		super(name);
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
