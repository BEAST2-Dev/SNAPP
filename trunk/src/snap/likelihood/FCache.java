
/*
 * File FCache.java
 *
 * Copyright (C) 2010 Remco Bouckaert, David Bryant remco@cs.auckland.ac.nz
 *
 * This file is part of SnAP.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * SnAP is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  SnAP is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SnAP; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */
package snap.likelihood;

import java.util.Arrays;
import java.util.Vector;

import snap.FMatrix;
import snap.NodeData;

/** cache for storing F matrices used for Site Probability calculation **/
public class FCache {

	/** class used to store F matrix + cache ID **/
	public class CacheObject {
		FMatrix m_F;
		public FMatrix getF() {
			return m_F;
		}
		int m_nCacheID;
		CacheObject(FMatrix F, int nCacheID) {
			m_F = F;
			m_nCacheID = nCacheID;
		} // c'tor
	} // CacheObject
	
	public class CacheObject2 extends CacheObject {
		int m_nCacheID2;
		CacheObject2(FMatrix F, int nCacheID, int nCacheID2) {
			super(F, nCacheID);
			m_nCacheID2 = nCacheID2;
		} // c'tor
	}
		
	/** matrix for storing leaf F-matrices **/
	CacheObject [][][] m_leafCache;
	/** vector for storing top-of-branch F-matrices **/
	Vector<CacheObject> m_TopOfBranche;
	Vector<CacheObject>[] m_TopOfBrancheID;
	/** two-dimensional vector array for storing bottom-of-branch F-matrices **/
	Vector<Vector<CacheObject2>> m_BottomOfBranche;
	Vector<CacheObject2>[] m_BottomOfBrancheID;
	
	@SuppressWarnings("unchecked")
	public FCache(int nNodeNrMax, int nRedsMax) {
		m_leafCache = new CacheObject[nNodeNrMax][nRedsMax][nRedsMax];
		m_TopOfBranche = new Vector<CacheObject>(1024,128);
		m_TopOfBrancheID = new Vector[nNodeNrMax];
		for (int i = 0; i < nNodeNrMax; i++) {
			m_TopOfBrancheID[i] = new Vector<CacheObject>(1024,128);
		}
		m_BottomOfBranche = new Vector<Vector<CacheObject2>>(1024,128);
		m_BottomOfBrancheID = new Vector[nNodeNrMax];
		for (int i = 0; i < nNodeNrMax; i++) {
			m_BottomOfBrancheID[i] = new Vector<CacheObject2>(1024,128);
		}
		g_nID = 0;
	} //c'tor

	/** remove all elements from cache **/
	public void clear() {
		//System.err.println("FCache::clear() dont know how to clear m_leafCache");
		for (int i = 0; i < m_TopOfBranche.size(); i++) {
			CacheObject o = m_TopOfBranche.elementAt(i);
			if (o != null) {
				o.m_F = null;
			}
		}
		//m_TopOfBranche.removeAllElements();
		for (int iNode = 0; iNode < m_BottomOfBrancheID.length; iNode++) {
			clearNode(iNode);
		}
		g_nID = 0;
	} // clear

	public void clearNode(int iNode) {
		Vector<CacheObject> topCacheObjects = m_TopOfBrancheID[iNode];
		for (int i = 0; i < topCacheObjects.size(); i++) {
			topCacheObjects.elementAt(i).m_F = null;
		}

		Vector<CacheObject2> bottomeCacheObjects = m_BottomOfBrancheID[iNode];
		for (int i = 0; i < bottomeCacheObjects.size(); i++) {
			bottomeCacheObjects.elementAt(i).m_F = null;
		}
	} // clearNode;
	
	/** counter to generate unique identifiers for pedigrees in cache **/
	static int g_nID = 0;
	/** generate the next pedigree identifier **/
	static synchronized int nextID() {
		return g_nID++;
	} // nextID

	CacheObject getLeafF(NodeData node, int numReds, int totalCount, boolean bHasDominantMarkers, SiteProbabilityCalculator spc) {
		
		if (m_leafCache[node.getNr()][numReds][totalCount] == null) {
			// it's not in the cache yet, so create the object
			//System.err.println("CACHE leaf = " + node.getNr() + " nReds = " + numReds + " total = " + totalCount);
			
			
			spc.doLeafLikelihood(node, numReds, totalCount, bHasDominantMarkers, false);
			//FMatrix Fb = node.cloneFb();
			FMatrix Fb = node.getFb();
			CacheObject o = new CacheObject(Fb, nextID());
			m_leafCache[node.getNr()][numReds][totalCount] = o; 
	    } else {
		    //System.err.println("RETRIEVE leaf = " + node.getNr() + " nReds = " + numReds + " total = " + totalCount);
	    }
							   
		return m_leafCache[node.getNr()][numReds][totalCount];
	} // getLeafF
	
	CacheObject getTopOfBrancheF(int nCacheID, NodeData node, double u, double v, double rate, Double [] coalescenceRate, SiteProbabilityCalculator spc) throws Exception {
		while (nCacheID >= m_TopOfBranche.size()) {
			m_TopOfBranche.add(null);
		}
		//System.err.println("Fetch nCacheID = " + nCacheID);
		CacheObject o = m_TopOfBranche.elementAt(nCacheID);//m_TopOfBrancheMap.get(nCacheID);
		if (o == null) {
			// it's not in the cache yet, so create the object
			spc.doTopOfBranchLikelihood(node, u, v, rate, coalescenceRate, false);
			//FMatrix Ft = node.cloneFt(); 
			FMatrix Ft = node.getFt(); 
			o = new CacheObject(Ft, nextID());
			m_TopOfBrancheID[node.getNr()].add(o);
			//System.err.println("Set TOB nCacheID = " + nCacheID);
			m_TopOfBranche.set(nCacheID, o);
		} else if (o.m_F == null) {
			// it's removed from the cache, so recalculate the F matrix
			spc.doTopOfBranchLikelihood(node, u, v, rate, coalescenceRate, false);
			o.m_F = node.getFt();
		}
		return o;
	} // getTopOfBrancheF


	CacheObject getBottomOfBrancheF(int nCacheID1, int nCacheID2, NodeData u1, NodeData u2, NodeData parent, SiteProbabilityCalculator spc) {
//		if (false) {
//			SiteProbabilityCalculator.doInternalLikelihood(u1, u2, parent, false);
//			return new CacheObject(parent.getFb(), -1);
//		}
		// try to fetch result from cache
		while (m_BottomOfBranche.size() <= nCacheID1) {
			m_BottomOfBranche.add(null);
		}
		Vector<CacheObject2> nodeCache2 = m_BottomOfBranche.elementAt(nCacheID1);
		if (nodeCache2 == null) {
			nodeCache2 = new Vector<CacheObject2>(1024,128);
			m_BottomOfBranche.set(nCacheID1, nodeCache2);
			//System.err.println("try to fetch from cache with IDs swapped " + nCacheID1 + "," + nCacheID2 );
			// not in cache, try to fetch from cache with IDs swapped
			// TODO: this does not get any hits as long as branch lengths (=t*gamma) are not part of the Key.
			// TODO: Fix this!
				spc.doInternalLikelihood(u1, u2, parent, false);
				//FMatrix Fb = parent.cloneFb();
				FMatrix Fb = parent.getFb();
				CacheObject2 o = new CacheObject2(Fb, nextID(), nCacheID2);
				m_BottomOfBrancheID[parent.getNr()].add(o);
				nodeCache2.add(o);
				return o;
		} 
		for (int i = 0; i < nodeCache2.size(); i++) {
			CacheObject2 o = nodeCache2.elementAt(i);
			if (o.m_nCacheID2 == nCacheID2) {
				if (o.m_F == null) {
					// it was removed from the cache, so recalculate it
					//System.err.println("it was removed from the cache, so recalculate it " + nCacheID1 + "," + nCacheID2 );
					spc.doInternalLikelihood(u1, u2, parent, false);
					o.m_F = parent.getFb();
				}
				return o;
			}
		}
		// it's not in the cache yet, so create the object
		//System.err.println("it's not in the cache yet, so create the object " + nCacheID1 + "," + nCacheID2 );
		spc.doInternalLikelihood(u1, u2, parent, false);
		//FMatrix Fb = parent.cloneFb();
		FMatrix Fb = parent.getFb();
		CacheObject2 o = new CacheObject2(Fb, nextID(), nCacheID2);
		m_BottomOfBrancheID[parent.getNr()].add(o);
		nodeCache2.add(o);
		return o;
	} // getBottomOfBrancheF
	
	/** print F matrix, for debugging purposes **/
	void printF(double[][]F) {
		for (int i = 1; i < F.length; i++) {
			System.err.println(Arrays.toString(F[i]));
		}
	}
} // class FCache
