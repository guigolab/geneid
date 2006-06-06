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

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Properties;
import java.util.StringTokenizer;


/**
 * Convenience class offering methods for easily accessing property files stored as files available via file system or URL. Makes it
 * easy to extract single property values or all of them as a property object.
 */
public class PropertiesUtil {

    /**
     * Creates property object from propertyFiles. Firstly attempts to open each propertyFile as a file, secondly attempts to open
     * it as a resource (available via CLASSPATH).
     * 
     * @return property object corresponding to contents of files (empty if can't open files or files are empty).
     */
    public static Properties createProperties(String[] propertyFiles) {

        Properties p = new Properties();
        for (int i = 0; i < propertyFiles.length; ++i) {
            String fileName = propertyFiles[i];
            if (fileName != null) {
                Properties p1 = createProperties(fileName);
                if (p1 != null) p.putAll(p1);
            } else
                System.err.println("WARNING: PropertiesUtil.createProperties(String[]): ignoring null filename");
        }
        return p;

    }

    /**
     * Creates property object from propertyFile(s). Firstly attempts to open property as a file, secondly attempts to open it as a
     * resource (available via CLASSPATH).
     * 
     * @param propertyFile list of 1 or more files separated by commas.
     * @return property object corresponding to file, or null if can't open file.
     */
    public static Properties createProperties(String propertyFile) {

        if (propertyFile == null) return null;

        if (propertyFile.indexOf(",") > 1) {

            StringTokenizer tokens = new StringTokenizer(propertyFile, ",");
            final int nTokens = tokens.countTokens();
            String[] files = new String[nTokens];
            for (int t = 0; t < nTokens; ++t)
                files[t] = tokens.nextToken();

            return createProperties(files);
        }

        URL url = null;

        File f = new File(propertyFile);
        if (f.exists()) {
            try {
                url = f.toURL();
            } catch (MalformedURLException e) {
                e.printStackTrace();
                url = null;
            }
        }

        if (url == null) {
            url = PropertiesUtil.class.getClassLoader().getResource(propertyFile);
        }

        if (url == null) {
            System.err.println("WARNING: PropertiesUtil.createProperties(String): Failed to load properties from file:"
                    + propertyFile);
            return null;
        } else {
            return createProperties(url);
        }
    }

    /**
     * Creates property object from url.
     * 
     * @return property object corresponding to file, or null if can't open resource.
     */
    public static Properties createProperties(URL url) {

        Properties p = new Properties();
        try {
        		InputStream is = url.openStream(); 
            p.load(is);
            is.close();
            return p;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }

    /**
     * @return value corresponding to key from file, null if can't open file or key not present.
     */
    public static String getProperty(String propertyFile, String key) {

        Properties p = createProperties(propertyFile);
        if (p == null) return null;
        return p.getProperty(key);
    }

    /**
     * @return value corresponding to key from url, null if can't open file or key not present.
     */
    public static String getProperty(URL url, String key) {

        Properties p = createProperties(url);
        if (p == null) return null;
        return p.getProperty(key);
    }

    /**
     * @return true if property value is is "true" (case insensitive), false if value is not "true", _default_ if proerty not set.
     */
    public static boolean booleanValue(Properties properties, String propertyName, boolean defaultValue) {

        String v = properties.getProperty(propertyName);
        if (v == null) return defaultValue;
        if (v.equalsIgnoreCase("true"))
            return true;
        else
            return false;
    }

    /**
     * Sorts pairs by key name and prints one name = value pair per line.
     */
    public static String toString(Properties p) {

        StringBuffer buf = new StringBuffer();
        ArrayList keys = new ArrayList(p.keySet());
        Collections.sort(keys);

        for (int i = 0; i < keys.size(); ++i) {
            String key = (String) keys.get(i);
            buf.append(key).append(" = ").append(p.getProperty(key)).append("\n");
        }

        return buf.toString();
    }

    /**
     * Filters the input set of properties only allowing those with keys that begin with _prefix_.
     * 
     * @param prefix key prefix to filter on.
     * @param in input properties
     * @return new properties object containing a subset of _in_ matching the condition.
     */
    public static Properties filterOnPrefix(String prefix, Properties in) {

        Properties out = new Properties();

        Enumeration keys = in.keys();
        while (keys.hasMoreElements()) {
            String key = (String) keys.nextElement();
            if (key.startsWith(prefix)) out.put(key, in.getProperty(key));
        }

        return out;
    }

    /**
     * Removes the string _prefix_ from the beggining pf any key that bigins with it.
     * 
     * @param prefix prefix to be removed from the beginning of any keys which begin with it.
     * @param in input properties
     * @return new properties object containing the same key to value values as _in_ except that all keys beginning with the
     *         _prefix_ will have been modified.
     */
    public static Properties removeKeyPrefix(String prefix, Properties in) {

        Properties out = new Properties();

        final int prefixLen = prefix.length();
        Enumeration keys = in.keys();
        while (keys.hasMoreElements()) {
            String key = (String) keys.nextElement();
            if (key.startsWith(prefix)) out.put(key.substring(prefixLen), in.getProperty(key));
        }

        return out;
    }

    /**
     * Write some properties to file. Optionally filter on property name.
     * 
     * @param props The properties to write.
     * @param fileName The file to write to.
     * @param pat If not null, only properties whose names match this string will be written.
     */
    public static void writeToFile(Properties props, String fileName, String pat) {

        try {

            OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(fileName));
            Enumeration en = props.propertyNames();
            
            while (en.hasMoreElements()) {
                String name = (String) en.nextElement();
                String value = (String) props.get(name);
                if (pat == null || name.indexOf(pat) >= 0) {
                    writer.write(name + "=" + value + "\n");
                }
            }

            writer.close();

        } catch (IOException e) {

            e.printStackTrace();
        }

    }

}// PropertyUtil
