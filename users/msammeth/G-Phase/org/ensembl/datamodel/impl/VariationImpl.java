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



package org.ensembl.datamodel.impl;

import java.util.ArrayList;
import java.util.List;

import org.ensembl.datamodel.SubmittedVariation;
import org.ensembl.datamodel.Variation;
import org.ensembl.datamodel.VariationFrequency;
import org.ensembl.driver.Driver;
import org.ensembl.util.StringUtil;

public class VariationImpl extends LocatableImpl implements Variation {

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


  public VariationImpl() {
    super();
  }

  public VariationImpl( Driver driver) {
    super( driver );
  }

  private String fivePrimeSequence;

  public String getFivePrimeSequence(){ return fivePrimeSequence; }

  public void setFivePrimeSequence(String fivePrimeSequence){ this.fivePrimeSequence = fivePrimeSequence; }

  private String threePrimeSequence;

  public String getThreePrimeSequence(){ return threePrimeSequence; }

  public void setThreePrimeSequence(String threePrimeSequence){ this.threePrimeSequence = threePrimeSequence; }


  /**
   * Synonym's are Strings.
   * @return list of zero or more synonyms collected from SubmittedVariations.
   */
  public List getSynonyms(){ 
    ArrayList synonyms = new ArrayList();
    if ( submittedVariations !=null) {
      final int nSubmittedVariations = submittedVariations.size();
      for(int i=0; i<nSubmittedVariations; ++i) {
        synonyms.add( ((SubmittedVariation)submittedVariations.get(i)).getSynonym() );
      }
    }
    return synonyms; 
  }

  private float heterozygosity;

  public float getHeterozygosity(){ return heterozygosity; }

  public void setHeterozygosity(float heterozygosity){ this.heterozygosity = heterozygosity; }

  private float heterozygosityStandardError;

  public float getHeterozygosityStandardError(){ return heterozygosityStandardError; }

  public void setHeterozygosityStandardError(float heterozygosityStandardError){ this.heterozygosityStandardError = heterozygosityStandardError; }

  private List submittedVariations;

  public List getSubmittedVariations(){ return submittedVariations; }

  public void setSubmittedVariations(List submittedVariations){ this.submittedVariations = submittedVariations; }

  /**
   * Allele's are Strings.
   * @return list of zero or more synonyms collected from SubmittedVariations.
   */
  public List getAlleles(){ 

    ArrayList alleles = new ArrayList();
    if ( submittedVariations !=null) {
      final int nSubmittedVariations = submittedVariations.size();
      for(int i=0; i<nSubmittedVariations; ++i) {

        final List freqs = ((SubmittedVariation)submittedVariations.get(i)).getFrequencies();
        final int nFreqs = freqs.size();
        for(int j=0; j<nFreqs; ++j) {
          String allele = ((VariationFrequency)freqs.get(j)).getAllele();
          if ( !alleles.contains( allele ) ) 
            alleles.add( allele );
        }
      }
    }
    return alleles; 
  }

  private List locations;

  public List getLocations(){ return locations; }

  public void setLocations(List locations){ this.locations = locations; }

  private int mapWeights;

  public int getMapWeights(){ return mapWeights; }

  public void setMapWeights(int mapWeights){ this.mapWeights = mapWeights; }

  private String validated;

  public String getValidated(){ return validated; }

  public void setValidated(String validated){ this.validated = validated; }

  public String toString() {    StringBuffer buf = new StringBuffer();
    buf.append("[");
    buf.append("internalID=").append(internalID).append(", ");
    buf.append("fivePrimeSequence=").append( StringUtil.briefOrUnset(fivePrimeSequence) ).append(", ");
    buf.append("threePrimeSequence=").append( StringUtil.briefOrUnset(threePrimeSequence) ).append(", ");
    buf.append("heterozygosity=").append(heterozygosity).append(", ");
    buf.append("heterozygosityStandardError=").append(heterozygosityStandardError).append(", ");

    buf.append("validated=").append(validated).append(", ");

    buf.append("locations=").append( StringUtil.sizeOrUnset(locations) ).append(", ");

    buf.append("submittedVariations=").append( StringUtil.sizeOrUnset(submittedVariations) ).append(", ");

    buf.append("location=").append( location ).append(", ");
    buf.append("]");

    return buf.toString();
  }
  

  

}
