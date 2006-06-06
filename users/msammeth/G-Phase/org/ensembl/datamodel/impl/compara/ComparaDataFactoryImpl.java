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
package org.ensembl.datamodel.impl.compara;

import org.ensembl.datamodel.DnaDnaAlignFeature;
import org.ensembl.datamodel.FeaturePair;
import org.ensembl.datamodel.compara.ComparaDataFactory;
import org.ensembl.datamodel.compara.DnaFragment;
import org.ensembl.datamodel.compara.GenomeDB;
import org.ensembl.datamodel.compara.GenomicAlign;
import org.ensembl.datamodel.compara.MethodLink;
import org.ensembl.datamodel.impl.DnaDnaAlignFeatureImpl;
import org.ensembl.datamodel.impl.FeaturePairImpl;
import org.ensembl.driver.plugin.compara.ComparaMySQLDriver;

public class ComparaDataFactoryImpl implements ComparaDataFactory {
  private final ComparaMySQLDriver driver;

  public ComparaDataFactoryImpl(){
    driver = null;
  }

  public ComparaDataFactoryImpl(ComparaMySQLDriver driver) {
    this.driver = driver;
  }

  public GenomicAlign createGenomicAlign() {
    return new GenomicAlignImpl();
  }

  public GenomeDB createGenomeDB() {
    return new GenomeDBImpl();
  }
  
  public DnaFragment createDnaFragment() {
    return new DnaFragmentImpl();
  }
  
  public FeaturePair createFeaturePair() {
    return new FeaturePairImpl();
  }
  
  public DnaDnaAlignFeature createDnaDnaAlignFeature() {
    return new DnaDnaAlignFeatureImpl();
  }//end createDnaDnaAlignFeature
  
  public MethodLink createMethodLink() {
    return new MethodLinkImpl();
  }//end createMethodLink
} // EnsemblDataFactoryImpl 

