/**
 * Copyright (c) 2004-2006 Regents of the University of California.
 * See "license-prefuse.txt" for licensing terms.
 */
package prefuse.data.util;

import java.util.Iterator;

import prefuse.data.Node;

/**
 * A depth-first iterator over the subtree rooted at given node.
 * 
 * @author <a href="http://jheer.org">jeffrey heer</a>
 */
public class TreeNodeIterator implements Iterator {

    private Node m_node;
    private Node m_root;
    
    /**
     * Create a new TreeNodeIterator over the given subtree.
     * @param root the root of the subtree to traverse
     */
    public TreeNodeIterator(Node root) {
        m_root = root;
        m_node = root;
    }
    
    /**
     * @see java.util.Iterator#hasNext()
     */
    public boolean hasNext() {
        return m_node != null;
    }

    /**
     * @see java.util.Iterator#next()
     */
    public Object next() {
        Node n, c;
        if ( (c=m_node.getChild(0)) != null ) {
            // do nothing
        } else if ( m_node!=m_root && (c=m_node.getNextSibling()) != null ) {
            // do nothing
        } else {
            c = m_node.getParent();
            while ( c!=m_root && c != null ) {
                if ( (n=c.getNextSibling()) != null ) {
                    c = n;
                    break;
                }
                c = c.getParent();
            }
            if ( c == m_root )
                c = null;
        }
        n = m_node;
        m_node = c;
        return n;
    }

    /**
     * Throws an UnsupportedOperationException
     * @see java.util.Iterator#remove()
     */
    public void remove() {
        throw new UnsupportedOperationException("Remove not supported");
    }

} // end of class TreeNodeIterator
