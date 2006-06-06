/*
Copyright (C) 2001 EBI, GRL

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

package org.ensembl.datamodel.variation;

import java.util.HashMap;
import java.util.Map;


/**
 * TranscriptVariation consequence types.
 * 
 * Usage: 
 * 
 * Either use one of the public final static instances directly
 * or use the factory method ConsequenceType.createConsequenceType(String s) where s is one
 * of "INTRONIC", "UPSTREAM", "DOWNSTREAM", "SYNONYMOUS_CODING",
 * "NON_SYNONYMOUS_CODING", "FRAMESHIFT_CODING", "5PRIME_UTR",
 * "3PRIME_UTR" 
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 */
public class ConsequenceType {
	
	private static Map types = new HashMap();
	
	public final static ConsequenceType INTRONIC = new ConsequenceType("INTRONIC");
  public final static ConsequenceType UPSTREAM = new ConsequenceType("UPSTREAM");
  public final static ConsequenceType DOWNSTREAM = new ConsequenceType("DOWNSTREAM");
  public final static ConsequenceType SYNONYMOUS_CODING = new ConsequenceType("SYNONYMOUS_CODING");
  public final static ConsequenceType NON_SYNONYMOUS_CODING = new ConsequenceType("NON_SYNONYMOUS_CODING");
  public final static ConsequenceType FRAMESHIFT_CODING = new ConsequenceType("FRAMESHIFT_CODING");
  public final static ConsequenceType FIVE_PRIME_UTR = new ConsequenceType("5PRIME_UTR");
  public final static ConsequenceType THREE_PRIME_UTR = new ConsequenceType("3PRIME_UTR");
  
  
	/**
	 * String representation of this type.
	 */
	public final String string;
	
	private ConsequenceType(String string) {
		this.string = string;
		types.put(string, this);
	}
	
	/**
	 * Retrieves one of the final static instances matching the string,
	 * otherwise returns null.
	 * @param string consequence type.
	 * @return one of the static instances or null.
	 */
	public static ConsequenceType createConsequenceType(String string) {
		Object o = types.get(string);
		return (o!=null) ? (ConsequenceType) o : null;
	} 
	
	/**
	 * @return the string representing this consequence type.
	 */
	public String toString() {
		return string;
	}
	
  
}
