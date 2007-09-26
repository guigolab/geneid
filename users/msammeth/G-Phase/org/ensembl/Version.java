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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package org.ensembl;

/**
 * Command line utility that prints version information about ensj.jar
 * library. Does not work when run outside ensj.jar.  
 *
 * <p><code>prompt&gt;java org.ensembl.Version</code></p>
 *
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 * @version $Revision: 1.1 $
 */
public class Version{

  
  /**
   * Prints Specification-Title, Specification-Version and
   * Specification-Title to standard out.  
   * @param args ignored
   */
  public static void main(String[] args) {

    // We need to create an instance of a class from the org.ensembl package
    // before we can see the manifest information about the package.

    Package p = new Version().getClass().getPackage();
    
    System.out.println("Specification-Title: "
                       + p.getSpecificationTitle());

    System.out.println("Specification-Version: "
                       + p.getSpecificationVersion());

    System.out.println("Implementation-Version: "
                       + p.getImplementationVersion());
    
  }
} // Version
