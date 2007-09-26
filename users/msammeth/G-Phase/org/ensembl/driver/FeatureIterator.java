/*
 Copyright (C) 2001 EBI, GRL

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
package org.ensembl.driver;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.ensembl.datamodel.Location;

/**
 * Provides an iterator over features from a FeatureAdaptor.
 * 
 * <p>
 * The features are retrieved in batches and buffered. This allows for a
 * compromise between runtime speed and memory usage efficiency. Increasing the
 * buffer size potentially reduces the number of database calls required and
 * potentially increases the memory requirement.
 * </p>
 * 
 * <p>
 * Example usage: Retrieve all genes in the specified location using a buffer of
 * upto 100 genes.
 * </p>
 * 
 * <p>
 * <code>Iterator iter = new FeatureIterator(geneAdaptor,100,loc);</code>
 * </p>
 * 
 * <p>
 * The class is thread safe.
 * </p>
 * 
 * TODO concurrent database access. TODO support massive datasets such as SNPs
 * by retrieving ids in batches.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 * 
 */
public class FeatureIterator implements Iterator {

	private FeatureProducer featureProducer;

	/**
	 * Produces zero or more internalIDs of features from the adaptor that match
	 * the location filter. TODO add a buffer to enable batch id retrieval of
	 * large numbers of ids currently nieve fetch all id implementation (this
	 * will potentially require huge amounts of memory for things like snps).
	 * 
	 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
	 */
	private class IDProducer {

		private int index = 0;

		private long[] ids;

		private Location location;

		private FeatureAdaptor adaptor;

		private IDProducer(long[] internalIDs) {
			this.ids = internalIDs;
		}
		
		private IDProducer(FeatureAdaptor adaptor, Location location)
				throws AdaptorException {
			this.adaptor = adaptor;
			this.location = location;
			if (location == null)
				ids = adaptor.fetchInternalIDs();
			else
				ids = adaptor.fetchInternalIDs(location);
		}

		public long nextID() {
			if (hasNextID())
				return ids[index++];
			else
				throw new NoSuchElementException();
		}

		public boolean hasNextID() {
			return index < ids.length;
		}
	}

	/**
	 * Produces zero or more Features from the adaptor that match the location
	 * filter. Uses the IDProducer to create a stream of internalIDs to
	 * retrieve.
	 */
	private class FeatureProducer {

		private int featureBufferSize;

		private boolean loadChildren;

		private long[] bufID;

		private List buf;

		private int pos = 0;

		private int size = 0;

		private FeatureAdaptor adaptor;

		private IDProducer idProducer;

		private FeatureProducer(FeatureAdaptor adaptor, int featureBufferSize,
				boolean loadChildren, Location location)
				throws AdaptorException {

			this.adaptor = adaptor;
			this.featureBufferSize = featureBufferSize;
			this.loadChildren = loadChildren;
			bufID = new long[featureBufferSize];
			this.idProducer = new IDProducer(adaptor, location);
		}

		/**
		 * @param adaptor2
		 * @param featureBufferSize2
		 * @param loadChildren2
		 * @param internalIDs
		 */
private FeatureProducer(FeatureAdaptor adaptor, 
				int featureBufferSize, 
				boolean loadChildren, 
				long[] internalIDs) {
    this.adaptor = adaptor;
    this.featureBufferSize = featureBufferSize;
    this.loadChildren = loadChildren;		
    bufID = new long[featureBufferSize];
		this.idProducer = new IDProducer(internalIDs);
	}
		private boolean hasNextFeature() {
			if (pos >= size)
				updateBuffer();
			return pos < size;
		}

		private Object nextFeature() {
			if (hasNextFeature())
				return buf.get(pos++);
			else
				throw new NoSuchElementException();
		}

		private void updateBuffer() {

			Arrays.fill(bufID, 0);
			int c = 0;
			while (idProducer.hasNextID() && c < featureBufferSize)
				bufID[c++] = idProducer.nextID();

			if (c > 0) {

				long[] tmp = bufID;
				if (c < featureBufferSize) {
					tmp = new long[c];
					System.arraycopy(bufID, 0, tmp, 0, c);
				}
				try {
					buf = adaptor.fetch(tmp, loadChildren);
					size = buf.size();
				} catch (AdaptorException e) {
					throw new RuntimeException(e);
				}
			} else {
				size = 0;
			}
			pos = 0;

		}

	}

	/**
	 * Create an iterator over all features in the adaptor. The order is
	 * unspecified.
	 * 
	 * @param adaptor
	 *            adaptor to retrieve features from.
	 * @throws AdaptorException
	 */
	public FeatureIterator(FeatureAdaptor adaptor, int featureBufferSize,
			boolean loadChildren) throws AdaptorException {
		this(adaptor, featureBufferSize, loadChildren, (Location) null);
	}

	/**
	 * Create an iterator over features in the adaptor that overlap the
	 * specified location. The order is unspecified.
	 * 
	 * @param adaptor
	 *            adaptor to retrieve features from.
	 * @param featureBufferSize
	 * @param loadChildren
	 *            whether to preload child data.
	 * @param location
	 *            location filter.
	 * @throws AdaptorException
	 */
	public FeatureIterator(FeatureAdaptor adaptor, int featureBufferSize,
			boolean loadChildren, Location location) throws AdaptorException {

		this.featureProducer = new FeatureProducer(adaptor, featureBufferSize,
				loadChildren, location);
	}

	/**
	 * Create an iterator over features with the specified internalIDs. The
	 * order is unspecified.
	 * 
	 * @param adaptor
	 *            adaptor to retrieve features from.
	 * @param featureBufferSize
	 * @param loadChildren
	 *            whether to preload child data.
	 * @param internalIDs
	 *            internal IDs of the features to retrieve
	 * @throws AdaptorException
	 */
	public FeatureIterator(FeatureAdaptor adaptor, int featureBufferSize,
			boolean loadChildren, long[] internalIDs) throws AdaptorException {

		this.featureProducer = new FeatureProducer(adaptor, featureBufferSize,
				loadChildren, internalIDs);
	}

	/**
	 * Not implemented.
	 * 
	 * @throws UnsupportedOperationException
	 *             this method is not supported
	 */
	public void remove() throws UnsupportedOperationException {
		throw new UnsupportedOperationException();

	}

	/**
	 * Returns whether there is another item to be returned via next().
	 * 
	 * @return true if the iterator has another feature.
	 */
	public synchronized boolean hasNext() {
		return featureProducer.hasNextFeature();
	}

	/**
	 * Returns the next element in the iteration.
	 * 
	 * @return the next Feature in the iteration.
	 * @throws NoSuchElementException
	 *             iteration has no more elements.
	 * @throws RuntimeException
	 *             problem occured reading the data from the database.
	 */
	public synchronized Object next() {
		return featureProducer.nextFeature();
	}

}
