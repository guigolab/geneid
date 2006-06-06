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



package org.ensembl.datamodel.impl;

import org.ensembl.datamodel.Population;
import org.ensembl.datamodel.SubmittedVariation;
import org.ensembl.datamodel.VariationFrequency;
import org.ensembl.util.StringUtil;

public class VariationFrequencyImpl implements VariationFrequency {

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


    public VariationFrequencyImpl() {
    }

  private String allele;
  private String otherAllele;


  public String getAllele(){ return allele; }

  public void setAllele(String allele){ this.allele = allele; }


  public String getOtherAllele(){ return otherAllele; }

  public void setOtherAllele(String otherAllele){ this.otherAllele = otherAllele; }

  private String type;

  public String getType(){ return type; }

  public void setType(String type){ this.type = type; }

  private float freq;

  public float getFreq(){ return freq; }

  public void setFreq(float freq){ this.freq = freq; }

  private float freqMax;

  public float getFreqMax(){ return freqMax; }

  public void setFreqMax(float freqMax){ this.freqMax = freqMax; }

  private float freqMin;

  public float getFreqMin(){ return freqMin; }

  public void setFreqMin(float freqMin){ this.freqMin = freqMin; }

  private int count;

  public int getCount(){ return count; }

  public void setCount(int count){ this.count = count; }

  private Population population;

  public Population getPopulation(){ return population; }

  public void setPopulation(Population population){ this.population = population; }

  private SubmittedVariation submittedVariation;

  public SubmittedVariation getSubmittedVariation(){ return submittedVariation; }

  public void setSubmittedVariation(SubmittedVariation submittedVariation){ this.submittedVariation = submittedVariation; }




  public String toString() {
    StringBuffer buf = new StringBuffer();
    buf.append("[");
    buf.append("allele=").append(allele).append(", ");
    buf.append("otherAllele=").append(otherAllele).append(", ");
    buf.append("type=").append(type).append(", ");
    buf.append("freq=").append(freq).append(", ");
    buf.append("freqMax=").append(freqMax).append(", ");
    buf.append("freqMin=").append(freqMin).append(", ");
    buf.append("count=").append(count).append(", ");
    buf.append("population=").append(population.getInternalID()).append(", ");
    buf.append("submittedVariation=").append(StringUtil.setOrUnset(submittedVariation) );
    buf.append("]");
    return buf.toString();
  }
}
