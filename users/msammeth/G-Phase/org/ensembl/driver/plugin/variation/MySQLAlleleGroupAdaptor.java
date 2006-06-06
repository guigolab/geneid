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

import org.ensembl.datamodel.variation.AlleleGroup;
import org.ensembl.datamodel.variation.VariationGroup;
import org.ensembl.datamodel.variation.impl.AlleleGroupImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.variation.AlleleGroupAdaptor;
import org.ensembl.driver.variation.VariationDriver;
import org.ensembl.util.LongSet;

/**
 * This adaptor provides database connectivity for AlleleGroup objects.
 * AlleleGroups may be retrieved from the Ensembl variation database by several
 * means using this module.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 *  
 */
public class MySQLAlleleGroupAdaptor extends BasePersistentAdaptor implements
		AlleleGroupAdaptor {

	private final static String BASE_QUERY = "SELECT ag.allele_group_id, ag.variation_group_id, ag.population_id,"
			+ "          ag.name, s.name, ag.frequency, aga.allele, aga.variation_id"
			+ "   FROM   allele_group ag, source s"
			+ "   LEFT JOIN allele_group_allele aga"
			+ "   ON     aga.allele_group_id = ag.allele_group_id"
			+ "   WHERE  ag.source_id = s.source_id";

	/**
	 * Constructor.
	 * 
	 * @param vdriver
	 *            parent variation driver this adaptor belongs to.
	 */
	public MySQLAlleleGroupAdaptor(VariationDriver vdriver) {
		super(vdriver);
	}

	/**
	 * @see org.ensembl.driver.allele.AlleleGroupAdaptor#fetch(long)
	 */
	public AlleleGroup fetch(long internalID) throws AdaptorException {
		// left join allows allele groups without any alleles to be fetched
		String sql = BASE_QUERY + " AND ag.allele_group_id = " + internalID;
		return (AlleleGroup) fetchByQuery(sql);
	}

	/**
	 * @see org.ensembl.driver.variation.VariationGroupAdaptor#fetch(java.lang.String)
	 */
	public AlleleGroup fetch(String name) throws AdaptorException {
		String sql = BASE_QUERY + " AND    ag.name = \"" + name + "\"";
		return (AlleleGroup) fetchByQuery(sql);
	}

	/**
	 * @see org.ensembl.driver.variation.AlleleGroupAdaptor#fetch(org.ensembl.datamodel.variation.VariationGroup)
	 */
	public List fetch(VariationGroup variationGroup) throws AdaptorException {
		String sql = BASE_QUERY + " AND    ag.variation_group_id = "
				+ variationGroup.getInternalID()
				+ " ORDER BY ag.allele_group_id";
		return fetchListByQuery(sql);
	}

	/**
	 * @see org.ensembl.driver.variation.AlleleGroupAdaptor#TYPE
	 * @see org.ensembl.driver.Adaptor#getType()
	 */
	public String getType() throws AdaptorException {
		return TYPE;
	}

	protected Object createObject(ResultSet rs) throws SQLException,
			AdaptorException {
		if (rs.isAfterLast())
			return null;

		final long internalID = rs.getLong("allele_group_id");
		AlleleGroup ag = null;
		LongSet variationIDs = new LongSet();

		do {

			if (ag == null) {
				ag = new AlleleGroupImpl(vdriver, rs.getString("ag.name"), rs
						.getDouble("ag.frequency"), vdriver
						.getPopulationAdaptor().fetch(
								rs.getLong("ag.population_id")), rs
						.getString("s.name"), vdriver
						.getVariationGroupAdaptor().fetch(
								rs.getLong("ag.variation_group_id")));

				ag.setInternalID(internalID);
			}

			long variationID = rs.getLong("aga.variation_id");
			if (variationID > 1) {
				ag.addVariation(vdriver.getVariationAdaptor()
						.fetch(variationID), rs.getString("aga.allele"));
			}

		} while (rs.next() && rs.getLong("allele_group_id") == internalID);

		return ag;

	}
}
