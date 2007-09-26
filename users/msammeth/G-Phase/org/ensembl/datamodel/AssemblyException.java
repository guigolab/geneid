/*
 * Created on 10-Nov-2003
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package org.ensembl.datamodel;

import org.ensembl.datamodel.impl.BaseFeatureImpl;

/**
 * @author arne
 *
 * To change the template for this generated type comment go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
public class AssemblyException extends BaseFeatureImpl {

  /**
   * Used by the (de)serialization system to determine if the data 
   * in a serialized instance is compatible with this class.
   *
   * It's presence allows for compatible serialized objects to be loaded when
   * the class is compatible with the serialized instance, even if:
   *
   * <ul>
   * <li> the compiler used to compile the "serializing" version of the class
   * differs from the one used to compile the "deserialising" version of the
   * class.</li>
   *
   * <li> the methods of the class changes but the attributes remain the same.</li>
   * </ul>
   *
   * Maintainers must change this value if and only if the new version of
   * this class is not compatible with old versions. e.g. attributes
   * change. See Sun docs for <a
   * href="http://java.sun.com/j2se/1.4.2/docs/guide/serialization/">
   * details. </a>
   *
   */
  private static final long serialVersionUID = 1L;



	public Location originalLocation;
	public Location linkedLocation;
	public String type;
	
	public static void main(String[] args) {
	}

	public AssemblyException(Location org, Location linked, String type) {
		this.originalLocation = org;
		this.linkedLocation = linked;
		this.type = type;
	}
	
	public AssemblyException() {
	}
}
