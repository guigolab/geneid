/*
 * Created on Jul 22, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.PrintStream;
import java.util.Date;

/**
 * 
 * 
 * @author msammeth
 */
public class Logger {

	static boolean dateStamp= true;
	static PrintStream ps= null;
	static {
		ps= System.out;
//		try {
//			ps= new PrintStream(new FileOutputStream(System.getProperty("user.dir")+ File.separator+ "run.log"));
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
	}
	public static final void println(String message, boolean error) {
		
		if (dateStamp)
			message= new Date(System.currentTimeMillis())+"\t"+ message;
		if (error)
			message= "!!! "+ message;
		
		ps.println(message);
		ps.flush();
	}
	
	public static final void print(String message) {
		
		ps.print(message);
		ps.flush();
	}

}
