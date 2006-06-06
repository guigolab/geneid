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
package org.ensembl.util;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.util.StringTokenizer;
import java.util.logging.Logger;

/**
 * Utility functions for working with JDBC.
 */
public class JDBCUtil {

  private static final Logger logger =
    Logger.getLogger(JDBCUtil.class.getName());

  /**
   * @return String containing each column in current row of result set.
   */
  public static String toString(ResultSet rs) {

    StringBuffer buf = new StringBuffer();
    buf.append("[");
    try {
      ResultSetMetaData metaData = rs.getMetaData();
      int nColumns = metaData.getColumnCount();
      for (int i = 1; i <= nColumns; ++i) {
        buf.append(metaData.getColumnName(i));
        buf.append(" = ");
        buf.append(rs.getString(i));
        if (i < nColumns)
          buf.append(" , ");
      }
    } catch (SQLException e) {
      buf.append(e.getMessage());
      e.printStackTrace();
    }
    buf.append("]");

    return buf.toString();
  }

  /**
   * Formats sql for reading.
   */
  public static String beautifySQL(String rawSQL) {

    StringBuffer buf = new StringBuffer();

    int indent = 0;
    StringTokenizer tokens = new StringTokenizer(rawSQL);
    while (tokens.hasMoreTokens()) {
      String token = tokens.nextToken();
      String lowerCaseToken = token.toLowerCase();
      if (lowerCaseToken.equals("select")
        || lowerCaseToken.equals("where")
        || lowerCaseToken.equals("from")
        || lowerCaseToken.equals("group")
        || lowerCaseToken.equals("order")
        || lowerCaseToken.equals("having"))
        buf.append("\n\n").append(token).append("\n").append("  ");
      else
        buf.append(' ').append(token);
    }

    return buf.toString();
  }

  public static void main(String[] args) throws Exception {
    String rawSql = null;
    if (args.length > 0) {
      StringBuffer buf = new StringBuffer();
      for (int i = 0; i < args.length; ++i)
        buf.append(args[i]).append(" ");
      rawSql = buf.toString();
      ;
    } else {
      System.out.println("Enter SQL:");
      rawSql = new BufferedReader(new InputStreamReader(System.in)).readLine();
    }
    System.out.println(beautifySQL(rawSql));
  }


  
} // MiscJDBC
