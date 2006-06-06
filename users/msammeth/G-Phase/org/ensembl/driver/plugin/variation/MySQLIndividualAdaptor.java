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
package org.ensembl.driver.plugin.variation;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;

import org.ensembl.datamodel.variation.Individual;
import org.ensembl.datamodel.variation.Population;
import org.ensembl.datamodel.variation.impl.IndividualImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.variation.IndividualAdaptor;
import org.ensembl.driver.variation.VariationDriver;

/**
 * Implementation of IndividualAdaptor that fetches Individuals from an ensembl
 * database.
 * 
 * Uses a cache to improve performance.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 *  
 */
public class MySQLIndividualAdaptor extends BasePersistentAdaptor implements
		IndividualAdaptor {

	private final static String BASE_QUERY = "SELECT i.individual_id, i.name, i.description, i.population_id, i.gender, i.father_individual_id, i.mother_individual_id FROM   individual i";

	/**
	 *  
	 */
	public MySQLIndividualAdaptor(VariationDriver vdriver) {
		super(vdriver,1000);
	}

	/**
	 * @see org.ensembl.driver.variation.IndividualAdaptor#fetch(long)
	 */
	public Individual fetch(long internalID) throws AdaptorException {

		Object o = cache.get(internalID);
		if (o != null)
			return (Individual) o;

		String sql = BASE_QUERY + " WHERE  i.individual_id = " + internalID;
		return (Individual) fetchByQuery(sql);
	}

	/**
	 * @see org.ensembl.driver.variation.IndividualAdaptor#fetch(java.lang.String)
	 */
	public List fetch(String name) throws AdaptorException {
		// can't use cache because name is potentially 1 to many
		// and even if cache has >=1 etry matching name there might
		// be others in the db that are not in the cache.
		String sql = BASE_QUERY + " WHERE  i.name = \"" + name + "\"";
		return fetchListByQuery(sql);
	}

	/**
	 * @see org.ensembl.driver.variation.IndividualAdaptor#fetch(org.ensembl.datamodel.variation.Population)
	 */
	public List fetch(Population population) throws AdaptorException {
		String sql = BASE_QUERY + " WHERE  i.population_id = "
				+ population.getInternalID() + "";
		return fetchListByQuery(sql);
	}

	/**
	 * @see org.ensembl.driver.variation.IndividualAdaptor#TYPE
	 * @see org.ensembl.driver.Adaptor#getType()
	 */
	public String getType() throws AdaptorException {
		return TYPE;
	}

	protected Object createObject(ResultSet rs) throws SQLException,
			AdaptorException {
		if (rs.isAfterLast())
			return null;

		final long internalID = rs.getLong("i.individual_id");
		Individual individual = new IndividualImpl(vdriver, rs
				.getString("i.name"), rs.getString("i.description"), rs
				.getString("i.gender"), rs.getLong("i.population_id"), rs
				.getLong("i.father_individual_id"), rs
				.getLong("i.mother_individual_id"));
		individual.setInternalID(internalID);
		
		cache.put(individual, individual.getInternalID());

		rs.next();
		
		return individual;

	}

	/**
	 * Fetches children of this individual.
	 * 
	 * If gender is set then look for children with this individual as father or
	 * mother as appropriatte. If gender is unset then look for children where
	 * this individual is first considered as a mother. If there no children are
	 * found the repeat the query but time treating this individual as a father.
	 * 
	 * @see org.ensembl.driver.variation.IndividualAdaptor#fetch(org.ensembl.datamodel.variation.Individual)
	 */
	public List fetch(Individual parent) throws AdaptorException {

		if (parent.getGender() != null)
			return fetchByParent(parent, parent.getGender());

		List r = fetchByParent(parent, "Female");
		if (r.size() == 0)
			r = fetchByParent(parent, "Male");
		return r;

	}

	private List fetchByParent(Individual parent, String gender)
			throws AdaptorException {
		String sql = BASE_QUERY
				+ " WHERE "
				+ (gender.equals("Male") ? " i.father_individual_id = "
						: " i.mother_individual_id = ")
				+ parent.getInternalID();
		return fetchListByQuery(sql);
	}
}
