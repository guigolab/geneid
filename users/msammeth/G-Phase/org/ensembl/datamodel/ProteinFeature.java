/*
    Copyright (C) 2002  Frans Verhoef

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
 * An alignment between a part of the genome and a peptide.
 */
public interface ProteinFeature extends Locatable {
  Translation getTranslation();

  void setTranslation(Translation translation);

  /**
   * Returns translation internal id, this is the same
   * as translation.internalID if translation is available.
   *
   * @return internalID of the translation.
   */
  long getTranslationInternalID();

  /**
   * Sets translation internal id, also sets the
   * translation.internalID if translation is available.
   */
  void setTranslationInternalID(long translationInternalID);

  int getTranslationStart();

  void setTranslationStart(int translationStart);

  int getTranslationEnd();

  void setTranslationEnd(int translationEnd);

  int getPeptideEnd();

  void setPeptideEnd(int peptideEnd);

  int getPeptideStart();

  void setPeptideStart(int peptideStart);

    Analysis getAnalysis();

    void setAnalysis(Analysis analysis);

    String getDisplayName();

    void setDisplayName(String displayName);

    /**
     * @return score or Double.NaN if score not set 
     */
    double getScore();

    void setScore(double score);

    /**
     * @return evalue or Double.NaN if evalue not set 
     */
    double getEvalue();

    void setEvalue(double evalue);

    /**
     * @return percentage identity or Integer.MIN_VALUE of not set. 
     */
    double getPercentageIdentity();

    void setPercentageIdentity(double percentageIdentity);

    /** @link dependency */
    /*# Translation lnkTranslation; */
}
