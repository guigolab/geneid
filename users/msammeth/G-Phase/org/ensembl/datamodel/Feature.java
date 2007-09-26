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

/**
 * Simple genomic feature. 
 */
public interface Feature extends Locatable {

    String getDisplayName();
    void setDisplayName(String displayName);

    String getDescription();

    void setDescription(String description);

     /**
     * Analysis which created this prediction.
     */

    Analysis getAnalysis();

  /**
   * @return analysis.internalID or 0 if analysis unset.
   */
  long getAnalysisID();

    /**
     * Analysis which created this prediction.
     */
    void setAnalysis(Analysis analysis);

    
    
}
