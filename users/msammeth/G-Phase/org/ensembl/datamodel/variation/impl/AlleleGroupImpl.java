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

package org.ensembl.datamodel.variation.impl;

import java.util.ArrayList;
import java.util.List;

import org.ensembl.datamodel.impl.PersistentImpl;
import org.ensembl.datamodel.variation.AlleleGroup;
import org.ensembl.datamodel.variation.Population;
import org.ensembl.datamodel.variation.Variation;
import org.ensembl.datamodel.variation.VariationGroup;
import org.ensembl.driver.variation.VariationDriver;

/**
 * Allele Group Implementation.
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 */
public class AlleleGroupImpl extends PersistentImpl implements AlleleGroup {

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



  private List variations = new ArrayList();
  
  private List alleles = new ArrayList();
  
  private String name;
  
  private double frequency;
  
  private Population population;
  
  private String source;
  
  private VariationGroup variatioGroup;

private VariationDriver vdriver;
  

	/**
	 * @param vdriver
	 * @param name
	 * @param frequency
	 * @param population
	 * @param source
	 * @param variatioGroup
	 */
	public AlleleGroupImpl(VariationDriver vdriver, String name, double frequency,
			Population population, String source, VariationGroup variatioGroup) {
		this.vdriver = vdriver;
		this.name = name;
		this.frequency = frequency;
		this.population = population;
		this.source = source;
		this.variatioGroup = variatioGroup;
	}
  /**
   * @see org.ensembl.datamodel.variation.AlleleGroup#getName()
   */
  public String getName() {
    return name;
  }

  /**
   * @see org.ensembl.datamodel.variation.AlleleGroup#getVariationGroup()
   */
  public VariationGroup getVariationGroup() {
    return variatioGroup;
  }

  /**
   * @see org.ensembl.datamodel.variation.AlleleGroup#getPopulation()
   */
  public Population getPopulation() {
    return population;
  }

  /**
   * @see org.ensembl.datamodel.variation.AlleleGroup#getSource()
   */
  public String getSource() {
    return source;
  }

  /**
   * @see org.ensembl.datamodel.variation.AlleleGroup#getFrequency()
   */
  public double getFrequency() {
    return frequency;
  }

  /**
   * @see org.ensembl.datamodel.variation.AlleleGroup#addVariation(org.ensembl.datamodel.variation.Variation, java.lang.String)
   */
  public void addVariation(Variation variation, String allele) {
    variations.add(variation);
    alleles.add(allele);
  }

  /**
   * @see org.ensembl.datamodel.variation.AlleleGroup#getVariations()
   */
  public List getVariations() {
    return variations;
  }

  /**
   * @see org.ensembl.datamodel.variation.AlleleGroup#getAlleles()
   */
  public List getAlleles() {
    return alleles;
  }

}
