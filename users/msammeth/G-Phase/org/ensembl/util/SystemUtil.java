/*
 Copyright (C) 2003 EBI, GRL

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
package org.ensembl.util;

/**
 * Utility methods for looking at system settings and command line program.
 * 
 * <p>
 * Can be used from the command line to determine where resources, e.g. classes
 * will be loaded from. Delegates parameter to runtimeResourceLocation(String).
 * e.g.
 * <code>java org.ensembl.util.SystemUtil "org.w3c.dom.Node" </code>
 * </p>
 * @see #runtimeResourceLocation(String) for more information on command line parameters
 */

import java.io.File;
import java.net.URL;
import java.util.StringTokenizer;

public class SystemUtil {

	/**
	 * Converts each element in the array to a string via element.toString() and
	 * creates a formated string containing all the strings.
	 * 
	 * Format:
	 * 
	 * <pre>
	 * 
	 *  0 = item0
	 *  1 = item2
	 *  3 = item3
	 *  ...
	 *  
	 * </pre>
	 */
	public static String toString(Object[] array) {
		StringBuffer buf = new StringBuffer();
		for (int i = 0; i < array.length; ++i) {
			buf.append(i).append(" = ").append(array[i]).append("\n");
		}
		return buf.toString();
	}

	/**
	 * Converts array into a string representation.
	 * 
	 * String format: <i>index = item </i>
	 * 
	 * <pre>
	 * 
	 *  0 = item0
	 *  1 = item2
	 *  3 = item3
	 *  ...
	 *  
	 * </pre>
	 */
	public static String toString(int[] array) {
		StringBuffer buf = new StringBuffer();
		for (int i = 0; i < array.length; ++i) {
			buf.append(i).append(" = ").append(array[i]).append("\n");
		}
		return buf.toString();
	}

	public static String getClasspath() {
		StringBuffer buf = new StringBuffer();

		StringTokenizer paths = new StringTokenizer(System.getProperty(
				"java.class.path", "."), ":");
		while (paths.hasMoreTokens())
			buf.append(paths.nextToken()).append("\n");

		return buf.toString();
	}

	/**
	 * Memory usage for JVM, useful for debugging.
	 * 
	 * @return string containing total and free memory usage for JVM.
	 */
	public static final String memoryStatus() {
		Runtime rt = Runtime.getRuntime();
		long total = rt.totalMemory();
		long free = rt.freeMemory();
		long used = total - free;
		return "memory: total = " + total + ", free = " + free + ", used = "
				+ used;

	}

	public static String environmentDump() {
		return "Ensj classes loaded from (sample class shown): " + runtimeResourceLocation("org.ensembl.Example")
		+ "\nMySQL classes loaded from (sample class shown): "  + runtimeResourceLocation("com.mysql.jdbc.Driver")
		 + "\n\nCLASSPATH = " + SystemUtil.getClasspath().toString();
	}

	
	/**
	 * Attempts to resolve the src to a URL.
	 * 
	 * This is useful for checking where a class or resource
	 * file is being loaded from by the classloadeder. 
	 * <ul>e.g. src = ...
	 * <li>org.ensembl.Example</li>
	 * <li>org/ensembl/Example.class</li>
	 * <li>com.mysql.jdbc.Driver</li>
	 * </ul>
	 * 
	 * Performs some magic to try to resolve the name if necessary,
	 * e.g. org.ensembl.Example will resolve to org/ensembl/Example.class
	 * @param src src file path or qualified class name.
	 * @return resolved URL or null if none found.
	 */
	public static URL runtimeResourceLocation(String src) {
		
		URL url = null;
		
		// try the following path variations for src 
		String[] variants = { src, src + ".class",
				src.replaceAll("\\.", File.separator),
				src.replaceAll("\\.", File.separator) + ".class" };
		for (int i = 0; url==null && i < variants.length; i++) 
			url = SystemUtil.class.getClassLoader().getResource(variants[i]);
		
		return url;
	}

	public static void main(String[] args) {
		

		if (args.length > 0) {
			System.out.println(args[0] + " is loaded from path: " + runtimeResourceLocation(args[0]));
		}
		else {
			System.out.println(environmentDump());
		}
	}
}
