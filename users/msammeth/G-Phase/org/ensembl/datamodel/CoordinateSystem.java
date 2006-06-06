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

import org.ensembl.datamodel.impl.PersistentImpl;
import org.ensembl.util.StringUtil;

/**
 * A genomic co-ordinate system.
 * 
 * <p>
 * In the database all Coordinate systems hava a name, version and are
 * either default or not default. When constructing a CoordinateSystem
 * the user will need to specify the coordinate system's
 * name and optionally it's version. If the name is unique or the desired version
 * is default in the database then the version is not needed. If the name
 * is ambiguous and a none default version is required then the user must
 * specify both the name and version.</p> 
 * 
 * 
 * <ul>This is a list of common coordinate system names.
 *
 * <li>"contig" - A contiguous sequence. ``N''s should be rare and of short
 * length. Can serve as your basic sequence holder.
 * </li>
 * 
 * <li>"clone" - Should have a real BAC or PAC or maybe YAC behind it.  Might
 * not be contiguous.
 * </li> 
 * 
 * <li>"supercontig" - Assembled from smaller contiguous
 * sequences.  May have small gaps (eg between read pairs).
 * </li>
 *
 * <li>"chromosome" - Use it only for real chromosomes.  Or for alternative
 * sequences of reference chromosomes.
 * </li>
 *
 * <li>"chunk" - Artificial coordinate system to hold sequence regions for
 * technical reasons. Create, when none of the other coordinate systems
 * can hold your sequence (eg. You only have full length chromosomes as
 * coordinate system but they are to long to store) or when you have 2 real
 * sequence containing coordinate systems.
 * </li>
 *
 * </ul>
 * 
 * Coordinate systems are used to partially specify <code>Location</code>s.
 *
 * @see Location
 */
public class CoordinateSystem
  extends PersistentImpl
  implements Cloneable, Comparable {

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



  private String name;
  private String version;
  private int rank;
  private boolean def;
  private boolean isSequenceLevel;
  private int hash = -1;

  /** 
   * Creates a new CoordinateSystem object.
   * Name, version and internalID are unset.
   * Default is false.
   */
  public CoordinateSystem() {
  }

  /** 
   * Creates a new CoordinateSystem object. 
   * Version is unset and default is false.
   * @param name name of the coordinate system.
   */
  public CoordinateSystem(String name) {
    this.name = name;
  }

  /** 
   * Creates a new CoordinateSystem object. 
   * 
   * internalID is unset.
   * @param name name of the coordinate system.
   * @param version coordinate system version
   * @param def whether the coordinate system is default.
   */
  public CoordinateSystem(String name, String version, boolean def) {

    this.name = name;
    this.version = version;
    this.def = def;

  }

  /** 
   * Creates a new CoordinateSystem object.
   * @param name name of the coordinate system.
   * @param version coordinate system version
   * @param def whether the coordinate system is default.
   * @param internalID internal ID of the coordinate system in the
   * database.
   */
  public CoordinateSystem(
    String name,
    String version,
    boolean def,
    long internalID) {

    this(name, version, def);
    this.internalID = internalID;

  }

  /** 
   * Get the name of this co-ordinate system.
   * @return Co-ordinate system name.
   */
  public String getName() {
    return name;
  }

  /** 
   * Set the name of this co-ordinate system.
   * @param name New name.
   *
   */
  public void setName(String name) {
    this.name = name;
    hash = -1;
  }

  /** 
   * Get the version of this co-ordinate system.
   * @return Co-ordinate system version.
   */
  public String getVersion() {
    return version;
  }

  /** 
   * Set the version of this co-ordinate system.
   * @param version New version.
   */
  public void setVersion(String version) {
    this.version = version;
    hash = -1;
  }

  /** 
   * Check if this co-ordinate system is the default for this species.
   * @return Co-ordinate system default.
   */
  public boolean isDefault() {
    return def;
  }

  /**
   * Set whether this co-ordinate system is the default.
   */
  public void setDefault(boolean b) {

    this.def = b;
    hash = -1;
  }

  /**
   * Test for equality with another Coordinate System.
   * Two co-ordinate systems are equal iff:
   *  -they have the same name
   *  -they have the same version (if set)
   * The comparison is case insensitive.
   * @param cs The CoordinateSystem to compare to.
   * @return true if the above conditions are met, false otherwise. 
   */
  public boolean equals(CoordinateSystem cs) {

    return (cs instanceof CoordinateSystem && cs.hashCode() == hashCode());
  }

  /**
   * Generate a hash code since we're overriding equals().
   * 
   * @return hash code.
   */
  public int hashCode() {

    if (hash == -1) {

      hash = 17;

      if (name != null) {
        hash = 37 * hash * name.hashCode();
      }

      if (version != null) {
        hash = 37 * hash * version.hashCode();
      }
    }
    return hash;
  }

  /**
   * Return a string representation of this object
   */
  public String toString() {

    StringBuffer buf = new StringBuffer();
    buf.append("[");
    buf.append("name=" + getName());
    buf.append(", version=" + getVersion());
    buf.append(", rank=" + getRank());
    buf.append(", sequenceLevel=" + isSequenceLevel());
    buf.append(",default=" + isDefault() + "]");

    return buf.toString();

  }

  /** 
   * Check if this co-ordinate system is the sequence-level co-ordinate system for this species.
   * @return Whether or not this is the sequence level co-ordinate system.
   */
  public boolean isSequenceLevel() {

    return isSequenceLevel;

  }

  /**
  * Set the flag indicating whether this co-ordinate system is the sequence level for this species.
  * @param f the new flag 
  */
  public void setSequenceLevel(boolean f) {

    isSequenceLevel = f;

  }

  /**
   * Orders by name then version.
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  public int compareTo(Object other) {

    if (other == null)
      return 1;

    CoordinateSystem otherCS = (CoordinateSystem) other;

    int tmp = name.compareTo(otherCS.name);
    if (tmp != 0)
      return tmp;

    return StringUtil.compare(version, otherCS.version);
  }

  /**
   * Return this coordinate system's rank.
   * @return the rank.
   */
  public int getRank() {
    return rank;
  }

  /**
   * Set this coordinate system's rank.
   * @param rank
   */
  public void setRank(int rank) {
    this.rank = rank;
  }

}
