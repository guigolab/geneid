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
package org.ensembl.datamodel.impl;

import java.util.List;
import java.util.logging.Logger;

import org.ensembl.datamodel.ExternalDatabase;
import org.ensembl.datamodel.ExternalRef;
import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.Transcript;
import org.ensembl.datamodel.Translation;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.Driver;
import org.ensembl.driver.ExternalDatabaseAdaptor;


public class ExternalRefImpl extends PersistentImpl implements ExternalRef {

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


  
  private static final Logger logger = Logger.getLogger( ExternalRefImpl.class.getName() );
  
  public ExternalRefImpl(Driver driver) {
    super( driver );
  }
  
  public List getSynonyms(){
    return synonyms;
  }
  
  
  public String getDescription(){
    return description;
  }
  
  public void setDescription(String description){
    this.description = description;
  }
  
  public String getDisplayID(){
    return displayID;
  }
  
  public void setDisplayID(String displayID){
    this.displayID = displayID;
  }
  
  public String getPrimaryID(){
    return primaryID;
  }
  
  public void setPrimaryID(String primaryID){
    this.primaryID = primaryID;
  }
  
  public String getVersion(){
    return version;
  }
  
  public void setVersion(String version){
    this.version = version;
  }
  
  public void setSynonyms(List synonyms){
    this.synonyms = synonyms;
  }
  
  //     public int getType(){
  //         return type;
  //     }
  
  //     public void setType(int type){
  //         this.type = type;
  //     }
  
  private long queryInternalID;
  
  public long getQueryInternalID(){ return queryInternalID; }
  
  public void setQueryInternalID(long queryInternalID){ this.queryInternalID = queryInternalID; }
  
  private int targetIdentity;
  
  public int getTargetIdentity(){ return targetIdentity; }
  
  public void setTargetIdentity(int targetIdentity){ this.targetIdentity = targetIdentity; }
  
  private int queryIdentity;
  
  public int getQueryIdentity(){ return queryIdentity; }
  
  public void setQueryIdentity(int queryIdentity){ this.queryIdentity = queryIdentity; }
  
  
  public long getExternalDbId() {
    if (externalDatabase!=null )
      externalDbId = externalDatabase.getInternalID();
    return externalDbId;
  }
  
  public void setExternalDbId(long externalDbId) {
    this.externalDbId = externalDbId;
  }
  
  private ExternalDatabase externalDatabase;
  
  public ExternalDatabase getExternalDatabase(){
    if ( externalDatabase==null && externalDbId != -1 && driver!=null ) lazyLoadExternalDatabase();
    if ( externalDatabase==null ) {
      externalDbId = -1;
    }
    return externalDatabase;
  }
  
  public void setExternalDatabase(ExternalDatabase externalDatabase){ this.externalDatabase = externalDatabase; }
  
  public String toString() {
    StringBuffer buf = new StringBuffer();
    buf.append("[");
    buf.append("{").append(super.toString()).append("}, ");
    buf.append("externalDbId_id=").append(getExternalDbId()).append(", ");
    buf.append("dbprimary_id=").append(getPrimaryID()).append(", ");
    buf.append("display_id=").append(getDisplayID()).append(", ");
    buf.append("version=").append(getVersion()).append(", ");
    buf.append("description=").append(getDisplayID()).append(", ");
    buf.append("synonyms=").append(synonyms);
    buf.append("]");
    
    return buf.toString();
  }
  
  private void lazyLoadExternalDatabase() {
    try {
      // Load a new database with same internalID
      ExternalDatabaseAdaptor xDbAdaptor = (ExternalDatabaseAdaptor)driver.getAdaptor("external_database");
      externalDatabase = xDbAdaptor.fetch(this.externalDbId);
    } catch (AdaptorException e) {
      logger.warning(e.getMessage());
      externalDbId = -1;
    }
  }
  
  /**
   * Get the GO linkage type.
   * @return The linkage type if this is a go_xref. If not, an empty string is returned.
   */
  public String getGoLinkageType() {
    
    return goLinkageType;
    
  }
  
  /**
   * Set the GO linkage type.
   * @param linkageType The new linkage type.
   */
  public void setGoLinkageType(String linkageType) {
    
    goLinkageType = linkageType;
    
  }
  
  
  // One of the following three objects should be set.
  private Gene gene;
  private Translation translation;
  private Transcript transcript;
  
  private List synonyms;
  private int type;
  private String description;
  private String displayID;
  private String primaryID;
  private String version;
  private long externalDbId;
  private String goLinkageType;
  
}

