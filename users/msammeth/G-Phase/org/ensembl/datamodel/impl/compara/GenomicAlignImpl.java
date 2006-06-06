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

import org.ensembl.datamodel.compara.DnaFragment;
import org.ensembl.datamodel.compara.GenomicAlign;
import org.ensembl.datamodel.compara.MethodLink;
import org.ensembl.datamodel.impl.PersistentImpl;

/**
 * I am part of the compara-analysis.
**/
public class GenomicAlignImpl extends PersistentImpl implements GenomicAlign
{

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


  int consensusDnaFragmentId;
  int queryDnaFragmentId;
  DnaFragment consensusDnaFragment;
  DnaFragment queryDnaFragment;
  int queryStrand;
  double score;
  int percentageId;
  String cigarString;
  int queryStart;
  int consensusStart;
  int queryEnd;
  int consensusEnd;
  MethodLink methodLink;
  int methodLinkInternalId;
  
  public int getQueryStrand() {
    return queryStrand;
  }//end getStrand
  
  public void setQueryStrand(int strand) {
    this.queryStrand = strand;
  }//end setStrand

  public DnaFragment getConsensusDnaFragment() {
    return consensusDnaFragment;
  }//end getConsensusDnaFragment
  
  public void setConsensusDnaFragment(DnaFragment dnaFrag) {
    this.consensusDnaFragment = dnaFrag;
  }//en dsetConsensusDnaFragment

  public DnaFragment getQueryDnaFragment() {
    return queryDnaFragment;
  }//end getQueryDnaFragment
  
  public void setQueryDnaFragment(DnaFragment dnaFrag) {
    this.queryDnaFragment = dnaFrag;
  }//end setQueryDnaFragment

  public int getConsensusDnaFragmentId() {
    return consensusDnaFragmentId;
  }//end getConsensusDnaFragmentId
  
  public void setConsensusDnaFragmentId(int dnaFrag) {
    this.consensusDnaFragmentId = dnaFrag;
  }//end setConsensusDnaFragmentId

  public int getQueryDnaFragmentId() {
    return queryDnaFragmentId;
  }//end getQueryDnaFragmentId
  
  public void setQueryDnaFragmentId(int dnaFrag) {
    this.queryDnaFragmentId = dnaFrag;
  }//end setQueryDnaFragmentId

  public String getCigarString() {
    return cigarString;
  }//end getCigarString
  
  public void setCigarString(String cigarString) {
    this.cigarString = cigarString;
  }//end setCigarString

  public double getScore() {
    return score;
  }//end getScore
  
  public void setScore(double score) {
    this.score = score;
  }//end getScore
  
  public int getPercentageId() {
    return percentageId;
  }//end getPercentageId
  
  public void setPercentageId(int percId) {
    percentageId = percId;
  }//end setPercentageId

  public String toString(){
    return 
      (new StringBuffer())
        .append("Genomic align(")
        .append(getInternalID())
        .append(")[")
        .append(getConsensusDnaFragment())
        .append(",")
        .append(getConsensusStart())
        .append("-")
        .append(getConsensusEnd())
        .append(",")
        .append(getQueryDnaFragment())
        .append(",")
        .append(getQueryStart())
        .append("-")
        .append(getQueryEnd())
        .append("]").toString();
  }//end toString
  
  public int getConsensusStart(){
    return consensusStart;
  }//end getConsensusStart
  
  public void setConsensusStart(int newValue){
    consensusStart = newValue;
  }//end setConsensusStart
  
  public int getQueryStart(){
    return queryStart;
  }//end getQueryStart
  
  public void setQueryStart(int newValue){
    queryStart = newValue;
  }//end setQueryStart
  
  public int getQueryEnd() {
    return queryEnd;
  }//end getQueryEnd

  public void setQueryEnd(int newValue) {
    queryEnd = newValue;
  }//end setQueryEnd

  public int getConsensusEnd() {
    return consensusEnd;
  }//end getConsensusEnd
  
  public void setConsensusEnd(int newValue) {
    consensusEnd = newValue;
  }//end getQueryEnd
  
  public MethodLink getMethodLink(){
    return methodLink;
  }//end getMethodLink
  
  public int getMethodLinkInternalId(){
    return methodLinkInternalId;
  }//end getMethodLinkInternalId
  
  public void setMethodLink(MethodLink newValue){
    methodLink = newValue;
  }//end setMethodLink
  
  public void setMethodLinkInternalId(int newValue){
    methodLinkInternalId = newValue;
  }//end setMethodLinkInternalId
}//end GenomicAlignImpl
