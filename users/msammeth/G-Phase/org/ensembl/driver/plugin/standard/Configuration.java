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

package org.ensembl.driver.plugin.standard;

import java.util.Enumeration;
import java.util.Properties;

/**
 * Properties object augmented with convenience methods to support both 
 * key=value and prefix.key=value pairs.
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 */
public class Configuration extends Properties {

	private static final long serialVersionUID = 1L;

	public Configuration(Properties properties) {
    putAll(properties);
  }

  public boolean hasKeyWithPrefix(String prefix) {
    for (Enumeration keysEnum = keys(); keysEnum.hasMoreElements();) {
      String key = (String) keysEnum.nextElement();
      if (key.startsWith(prefix))
        return true;
    }
    return false;
  }

  public void putProperty(String prefix, String rawKey, String value) {
    put(prefix+"."+rawKey,value);
  }

  public String getProperty(String prefix, String rawKey) {
    String value = getProperty(prefix + "." + rawKey);
    if (value == null)
      value = getProperty(rawKey);
    return value;
  }
}
