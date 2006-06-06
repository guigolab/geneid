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
package org.ensembl.datamodel.compara;

import org.ensembl.datamodel.Persistent;

/**
 * I am part of the compara-analysis. I keep information about
 * DnaDna alignments between species.
**/
public interface GenomicAlign extends Persistent{
  public int getQueryStrand(); 
  public void setQueryStrand(int strand);
  public DnaFragment getConsensusDnaFragment();
  public void setConsensusDnaFragment(DnaFragment dnaFrag);
  public DnaFragment getQueryDnaFragment();
  public void setQueryDnaFragment(DnaFragment dnaFrag);
  public int getConsensusDnaFragmentId();
  public void setConsensusDnaFragmentId(int dnaFrag);
  public int getQueryDnaFragmentId();
  public void setQueryDnaFragmentId(int dnaFrag);
  public String getCigarString();
  public void setCigarString(String cigarString);
  public double getScore();
  public void setScore(double score);
  public int getPercentageId();
  public void setPercentageId(int percId);
  public int getConsensusStart();
  public void setConsensusStart(int newValue);
  public int getQueryStart();
  public void setQueryStart(int newValue);
  public int getConsensusEnd();
  public void setConsensusEnd(int newValue);
  public int getQueryEnd();
  public void setQueryEnd(int newValue);
  public MethodLink getMethodLink();
  public int getMethodLinkInternalId();
  public void setMethodLink(MethodLink newValue);
  public void setMethodLinkInternalId(int newValue);
}
