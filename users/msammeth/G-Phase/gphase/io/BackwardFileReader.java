package gphase.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;

public class BackwardFileReader extends RandomAccessFile {

	public BackwardFileReader(String arg0, String arg1)
			throws FileNotFoundException {
		super(arg0, arg1);
		try {
			seek(new File(arg0).length());	// set pointer to the end of file
		} catch (IOException e) {
			e.printStackTrace();
		}	
		
	}

	public BackwardFileReader(File arg0, String arg1)
			throws FileNotFoundException {
		super(arg0, arg1);
		try {
			seek(arg0.length());	// set pointer to the end of file
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}

	public byte readByteBackward() throws IOException {
		seek(getFilePointer()- 1);
		byte b= readByte();
		seek(getFilePointer()- 1);
		return b;
	}
	
	public char readCharBackward() throws IOException {
		return (char) readByteBackward();
	}
	
	public String readLineBackward() throws IOException {
		StringBuffer sb= new StringBuffer();
		char c= readCharBackward();
		while (c!= '\n') {
			sb.insert(0, c);
			c= readCharBackward();
		}
		return sb.toString();
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	
}
