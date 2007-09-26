/*
 * Created on 12-Nov-2003
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package org.ensembl.driver;

import org.ensembl.datamodel.AssemblyException;
import org.ensembl.datamodel.Location;

/**
 * @author arne
 *
 */
public interface AssemblyExceptionAdaptor extends Adaptor {

  final static String TYPE = "assembly_exception";

  /**
  * Retrieve AssemblyExceptions that are relevant for given location.
  * 
  * <p>They have to be on the same sequence region with their linked Location.
  * </p>
  * @param loc
  * @return AssemblyExceptions that are relevant for given location
  * @throws AdaptorException
  */
  AssemblyException[] fetchLinks(Location loc) throws AdaptorException;

}
