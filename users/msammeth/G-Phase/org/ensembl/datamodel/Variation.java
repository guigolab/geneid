/*
    Copyright (C) 2002 EBI, GRL

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


package org.ensembl.datamodel;

import java.util.List;

/**
 * A variation such as a SNP.
 * 
 * <p>
 * Note that some variations do not map to the genome, in these
 * cases the location will be null.
 * </p>
 */
public interface Variation extends Locatable {
    String getFivePrimeSequence();

    void setFivePrimeSequence(String fivePrimeSequence);

    String getThreePrimeSequence();

    void setThreePrimeSequence(String threePrimeSequence);

    List getSynonyms();

    float getHeterozygosity();

    void setHeterozygosity(float heterozygosity);

    /**
     * @return error in range 0-1.
     */
    float getHeterozygosityStandardError();

    /**
     * @param heterozygosityStandardError error in range 0-1.
     */
    void setHeterozygosityStandardError(float heterozygosityStandardError);

    List getSubmittedVariations();

    void setSubmittedVariations(List submittedVariations);

		/**
		 * @return list of Strings.
		 * */
    List getAlleles();

    List getLocations();

    void setLocations(List locations);


    String getValidated();


    void setValidated(String validated);

    /**
     * The map weights.
     *
     * Map weights can have the values 1, 2, 3 or 10 where:
     * <pre>
     * mapweight =1 mapped to single position in genome
     * mapweight =2 mapped to 2 positions on a single chromosome
     * mapweight =3 mapped to 3-10 positions in genome (possible paralog hits)
     * mapweight =10 mapped to >10 positions in genome
     * </pre>
     * @return map weight.
     */
    int getMapWeights();

    /**
     * @see #getMapWeights() getMapWeights()
     */
    void setMapWeights(int mapWeights);

  void setLocation(Location location);

}
