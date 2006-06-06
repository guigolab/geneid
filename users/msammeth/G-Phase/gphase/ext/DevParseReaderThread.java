package gphase.ext;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Vector;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class DevParseReaderThread extends DevPipeReaderThread {

	protected Object origin= null;	// needed for refective invocation
	protected Vector commands= null, methods= null;
	/**
	 * Constructor for DevParseReaderThread.
	 * @param in
	 * @param out
	 */
	public DevParseReaderThread(InputStream in, OutputStream out) {
		super(in, out);
	}
	
	/**
	 * Constructor for DevParseReaderThread.
	 * @param in
	 * @param out
	 */
	public DevParseReaderThread(InputStream in, Object newOrigin) {

		this.in= in;
		this.origin= newOrigin;
	}	
	
	/**  
	 * The main method.  
	 * Reads continously from the input stream(s) as long as there is something to read and
	 * interrupt() is not called
	 */
	public void run() {

		
		BufferedReader buffy= new BufferedReader(
			new InputStreamReader(in));
		finished= false;
		String read= "";
		
		try {
		do {
			
			try {
				sleep(10);
			} catch (InterruptedException e) {
				; //nothing
			}
				
			// read line from stream
			try {				
				
//				read= null;
//				try {
					read = buffy.readLine();
//				} catch (ThreadDeath e) {
//					; // nothing, if stopped
//				}
				if (read== null) {
					finished= true;
					break;
				}
				parse(read);

			} catch (IOException e) {
				finished= true;
			}

		} while (!finished);

		} catch (Throwable e) {	// ThreadDeath
			finished= true;
		}
		
//		System.out.println("thread died: "+getName());
	}
	
	/**
	 * @see gphase.ext.DevPipeReaderThread#setFinished(boolean)
	 */
	public void setFinished(boolean finished) {
		
		super.setFinished(finished);
		commands= null;
		methods= null;
		
//		System.out.println("finished: "+getName());
	}

	
	protected void parse(String read) {
		
		if ((commands== null)|| (read== null))
			return;

		String readLC= read.toLowerCase(); 	// lower case comparison
		Object[] args= new Object[1];
//		System.out.println("*"+readLC+"*");
		for (int i= 0; (!finished)&& (i< commands.size()); ++i) {
			if (readLC.indexOf((String) commands.elementAt(i))!= -1) {
				args[0]= readLC;
				try {
					((Method) methods.elementAt(i)).invoke(origin, args);
				} catch (InvocationTargetException e) {
					; // nothing
				} catch (IllegalAccessException e) {
					; // nothing
				}
				break;
			} 
//			else System.out.println(readLC);
			
		}
//		if (readLC.indexOf("%]")!= -1)
//			System.out.println(readLC);
		
	}
	
	public void addCommand(String newCommand, Method newMethod) {
		
		if (commands== null) {
			commands= new Vector();
			methods= new Vector();
		}
		
		commands.add(newCommand.toLowerCase());
		methods.add(newMethod);
	}	

	/**
	 * Constructor for DevParseReaderThread.
	 * @param in
	 * @param in2
	 * @param out
	 */
	public DevParseReaderThread(
		InputStream in,
		InputStream in2,
		OutputStream out) {
		super(in, in2, out);
	}

}
