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

import java.util.List;

import org.ensembl.datamodel.impl.PersistentImpl;
import org.ensembl.datamodel.variation.Individual;
import org.ensembl.datamodel.variation.Population;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.RuntimeAdaptorException;
import org.ensembl.driver.variation.VariationDriver;

/**
 * The point of this class is....
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 */
public class IndividualImpl extends PersistentImpl implements Individual {

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



	private VariationDriver vdriver;

	private String name;

	private String description;

	private String gender;

	private Population population;

	private Individual father;

	private Individual mother;

	private List children;

	private long populationID;

	private long fatherID;

	private long motherID;

	/**
	 */
	public IndividualImpl(VariationDriver vdriver, String name,
			String description, String gender, long populationID,
			long fatherID, long motherID) {
		this.vdriver = vdriver;
		this.name = name;
		this.description = description;
		this.gender = gender;
		this.populationID = populationID;
		this.fatherID = fatherID;
		this.motherID = motherID;
	}

	/**
	 * @see org.ensembl.datamodel.variation.Individual#getName()
	 */
	public String getName() {
		return name;
	}

	/**
	 * @see org.ensembl.datamodel.variation.Individual#getDescription()
	 */
	public String getDescription() {
		return description;
	}

	/**
	 * @see org.ensembl.datamodel.variation.Individual#getGender()
	 */
	public String getGender() {
		return gender;
	}

	/**
	 * @see org.ensembl.datamodel.variation.Individual#getPopulation()
	 */
	public Population getPopulation() {
		if (population == null && populationID > 0 && vdriver!=null)
			try {
				population = vdriver.getPopulationAdaptor().fetch(populationID);
			} catch (AdaptorException e) {
				throw new RuntimeAdaptorException(
						"Failed to lazy load population of Individual "
								+ internalID, e);
			}
		return population;
	}

	/**
	 * Gets father, lazy loads if possible.
	 * 
	 * @see org.ensembl.datamodel.variation.Individual#getFather()
	 */
	public Individual getFather() {
		if (father == null && fatherID > 0)
			try {
				father = vdriver.getIndividualAdaptor().fetch(fatherID);
			} catch (AdaptorException e) {
				throw new RuntimeAdaptorException(
						"Failed to lazy load father of Individual "
								+ internalID, e);
			}
		return father;
	}

	/**
	 * @see org.ensembl.datamodel.variation.Individual#getMother()
	 */
	public Individual getMother() {
		if (mother == null && motherID > 0  && vdriver!=null)
			try {
				mother = vdriver.getIndividualAdaptor().fetch(motherID);
			} catch (AdaptorException e) {
				throw new RuntimeAdaptorException(
						"Failed to lazy load mother of Individual "
								+ internalID, e);
			}
		return mother;
	}

	/**
	 * @see org.ensembl.datamodel.variation.Individual#getChildren()
	 */
	public List getChildren() {
		if (children == null  && vdriver!=null)
			try {
				children = vdriver.getIndividualAdaptor().fetch(this);
			} catch (AdaptorException e) {
				throw new RuntimeAdaptorException(
						"Failed to lazy load children of Individual "
								+ internalID, e);
			}
		return children;
	}

	public String toString() {
		StringBuffer buf = new StringBuffer();

		buf.append("[");
		buf.append(super.toString());
		buf.append(", name=").append(name);
		buf.append(", description=").append(description);
		buf.append(", gender=").append(gender);
		buf.append(", population=").append(
				population != null ? population.getName() : "unset");
		buf.append(", mother=").append(
				mother != null ? mother.getName() : "unset");
		buf.append(", father=").append(
				father != null ? father.getName() : "unset");
		buf.append("]");

		return buf.toString();
	}

}
