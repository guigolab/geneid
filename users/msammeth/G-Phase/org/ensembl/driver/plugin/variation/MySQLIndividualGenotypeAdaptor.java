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

package org.ensembl.driver.plugin.variation;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;

import org.ensembl.datamodel.variation.Individual;
import org.ensembl.datamodel.variation.IndividualGenotype;
import org.ensembl.datamodel.variation.impl.IndividualGenotypeImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.variation.IndividualGenotypeAdaptor;

/**
 * Implementation of the IndividualGenotypeAdaptor that works
 * with ensembl databases.
 *
 * Uses a cache.
 *
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 */
public class MySQLIndividualGenotypeAdaptor extends BasePersistentAdaptor
		implements IndividualGenotypeAdaptor {

	private final static String BASE_QUERY = "SELECT individual_genotype_id, variation_id, allele_1, allele_2, individual_id FROM   individual_genotype"; 
	
	public MySQLIndividualGenotypeAdaptor(MySQLVariationDriver vdriver) {
		super(vdriver,1000);
	}

	/**
	 * @see org.ensembl.driver.plugin.variation.BasePersistentAdaptor#createObject(java.sql.ResultSet)
	 */
	protected Object createObject(ResultSet rs) throws SQLException,
			AdaptorException {
		
		if (rs.isAfterLast()) return null;
		
		IndividualGenotype ig = new IndividualGenotypeImpl(vdriver, 
                                                       rs.getString("allele_1"), 
                                                       rs.getString("allele_2"),
                                                       rs.getLong("variation_id"),
                                                       rs.getLong("individual_id"));
		ig.setInternalID(rs.getLong("individual_genotype_id"));

		cache.put(ig, ig.getInternalID());

    rs.next();

		return ig;
	}

	/**
	 * @see org.ensembl.driver.variation.IndividualGenotypeAdaptor#fetch(long)
	 */
	public IndividualGenotype fetch(long internalID) throws AdaptorException {
		Object o = cache.get(internalID);
		if (o!=null)
			return (IndividualGenotype) o;
		
		String sql = BASE_QUERY + " WHERE  individual_genotype_id = " + internalID;
		return (IndividualGenotype) fetchByQuery(sql);
	}

	/**
	 * @see org.ensembl.driver.variation.IndividualGenotypeAdaptor#fetch(org.ensembl.datamodel.variation.Individual)
	 */
	public List fetch(Individual individual) throws AdaptorException {
		String sql = BASE_QUERY + " WHERE  individual_id = " + individual.getInternalID();
		return fetchListByQuery(sql);
	}

	/**
	 * @see org.ensembl.driver.Adaptor#getType()
	 */
	public String getType() throws AdaptorException {
		return TYPE;
	}

}
