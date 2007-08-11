package gphase.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

public class FileTools {
	/**
	 * determines the line feed chars of a file by reading 2 lines from the file
	 * assumes that file has 2 lines and line is not longer than 10000 
	 * @return
	 */
	public static String determineLineFeed(File f) {
		try {
			BufferedReader reader= new BufferedReader(new FileReader(f));
			String line= reader.readLine();
			char[] buf= new char[line.length()+ 2];
			reader.close();
			reader= new BufferedReader(new FileReader(f));
			reader.read(buf, 0, line.length()+ 2);
			String result= "";
			if (Character.isWhitespace(buf[buf.length- 2]))
				result+= buf[buf.length- 2];
			if (Character.isWhitespace(buf[buf.length- 1]))
				result+= buf[buf.length- 1];
			reader.close();
			return result;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
}
