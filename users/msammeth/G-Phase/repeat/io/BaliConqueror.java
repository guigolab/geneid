/*
 * Created on Oct 19, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;

/**
 * 
 * 
 * @author micha
 */
public class BaliConqueror {
	
	static String pathToG6= "D:\\Eigene Dateien\\repeats\\data\\ref6"; 
	Method redirectMethod= null;
	Constructor constructor= null;
	String fileExt= null;
	
	public static void main(String[] args) {
	}
	
	public void conquer() {

			// retrieve sub-dirs		
		File dir= new File(pathToG6);
		String[] subDirs= dir.list();
		Arrays.sort(subDirs);

			// process subdirs
		boolean singleFile= false;
		boolean jumpDir= false;		
		boolean jump= false;
		for (int i= 0; i < subDirs.length; i++) {
			if (jumpDir&& !subDirs[i].equals("test"))
				continue;
			//jumpDir= false;			
			if (!singleFile)
				jumpDir= false;				
			System.out.println("\n== ["+subDirs[i]+"] ==");
			File sdir= new File(pathToG6+File.separator+subDirs[i]);
			String[] files= sdir.list();
			for (int j= 0; j < files.length; j++) {
				int idx= files[j].indexOf("_schema");			// unique marker "_schema" ?!?
				if (idx< 0)
					continue;
				String procFile= files[j].substring(0,idx);
				if (jump&& !procFile.equals("dead_ref6"))	
					continue;
				if (!singleFile)
					jump= false;				
				String absFile= pathToG6+File.separator+subDirs[i]+File.separator+procFile;
				System.out.print("processing "+procFile+"..");
				Object[] args= {absFile, fileExt};
				try {
					if (constructor== null) 
						redirectMethod.invoke(null, args);
					else {
						Object instance= constructor.newInstance(args);
						redirectMethod.invoke(instance, null);
					}
				} catch (InstantiationException e) {
					; // :)
				} catch (IllegalAccessException e) {
					; // :)
				} catch (InvocationTargetException e) {
					; // here if not read an alignment correctly
				}
				System.out.println("done.");
			}
		}
	
	}
	
	public BaliConqueror(Constructor newConstructor, Method newTarget, String newFileExt) {
		this.constructor= newConstructor;
		this.redirectMethod= newTarget;
		this.fileExt= newFileExt;
	}
}
