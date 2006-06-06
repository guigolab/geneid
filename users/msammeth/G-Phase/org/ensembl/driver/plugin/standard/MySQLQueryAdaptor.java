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
package org.ensembl.driver.plugin.standard;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.logging.Logger;

import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.ConfigurationException;
import org.ensembl.driver.QueryAdaptor;

/**
 * Provides the ability to execute SQL directly against the database
 * this adaptor is connected to.
 *
 * @author <a href="mailto: "Craig Melsopp</a>
 */
public class MySQLQueryAdaptor  extends BaseAdaptor implements QueryAdaptor{
  private static final Logger logger = Logger.getLogger(MySQLQueryAdaptor.class.getName());

    public MySQLQueryAdaptor (MySQLDriver driver){
	super(driver);
    }


    /**
     * Executes a piece of sql.
     *
     * Example usage: pass an SQL query to the adaptor for execution:
     * <br><code>ResultSet rs = (ResultSet)queryAdaptor.execute("SELECT * FROM
     * gene.") </code>
     * <br><code>Integer nRows = (Integer)queryAdaptor.execute("DELETE FROM
     * gene WHERE gene_internal_id=33.") </code>
     * * @param query object for which query.toString() returns an SQL statement.
     * @return a java.sql.ResultSet for select, and an Integer for update, insert or delete statements.
     * The Integer contains the number of rows affected.
     */
    public Object execute(Object query) throws AdaptorException {
      Connection conn = null;
      try {
        String sql = query.toString().toLowerCase();
        conn = getConnection();
				Object result = null;
				if ( sql.indexOf("update")!=-1
						|| sql.indexOf("insert")!=-1
						|| sql.indexOf("delete")!=-1 )
        	result = new Integer( conn.createStatement().executeUpdate(sql) );
        else 
					result = conn.createStatement().executeQuery(sql);
        return result;
      
      } catch (SQLException e) {
        throw new AdaptorException("Failed to execute query: "+query + ": "+ e.getMessage(),e);
      }finally {
        close( conn );
      }
    }


  public String getType(){
    return TYPE;
  }

  void configure() throws ConfigurationException {
  }
}// MySQLQueryAdaptor
