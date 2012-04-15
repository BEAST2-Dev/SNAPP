
/*

 * File NodeData.java
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
package snap;




import java.io.Serializable;

import beast.core.Description;
import beast.evolution.tree.Node;



@Description("This node has population parameter (coalescenceRate=2/theta) " +
		"and some specific members for performing SnAP (SNP and AFLP) " +
		"analysis.")
public class NodeData extends Node implements Serializable  {
	private static final long serialVersionUID = 1L;

	//private	double m_fCoalescenceRate = 1;
	
	//Number of individuals at or below this node.
//	public	int m_n;
	public int getSize() {
		return m_Fb.getSize();
	}
	//Number of species at or below this node
	//public int nspecies;
	//public int mintaxa; //Minimum id taxa below this node.
	//Vector of probabilities for the number of ancestral lineages at this node.
//	public double [] m_Nt;
//	public double [] m_Nb;
	//This is the likelihood multiplied by the lineage probabilities Pr(Ry | n,r ) x Pr(n)  in the paper.
	FMatrix m_Fb;
	public FMatrix getFb() {return m_Fb;}
	public void initFb(int n, int nReds) {
		m_Fb = new FMatrix(n, nReds);
	}
	public void initFb(int n, double [] FAsArray) {
		m_Fb = new FMatrix(n, FAsArray);
	}
	public void initFb(FMatrix F) {
		m_Fb = F;
	}
	FMatrix m_Ft;
	public FMatrix getFt() {return m_Ft;}
	public void initFt(FMatrix F) {
		m_Ft = F;
	}
	//public double [][] F; 


	//height.... only used for debugging
	//public double height;

	/* NB: length comes from basic_newick tree **/
	//public double length;

	/** Tree related methods **/
	//SSSTreeLikelihood m_tree;
	public NodeData getChild(int i) {
		if (i == 0) {
			return (NodeData) getLeft();
		} else if (i== 1) {
			return (NodeData) getRight();
		}
		return null;
	}
	public int getNrOfChildren() {
		if (getLeft() == null) {
			return 0;
		}
		if (getRight() == null) {
			return 1;
		}
		return 2;
	}

	public NodeData() {
		//m_n = -1;
		//cerr<<"YYY\t\tallocating\t"<<this<<endl;
		//m_children = new Vector<NodeData>();
		m_Fb = new FMatrix();
		m_Ft = new FMatrix();
//		m_nTaxonID = 0;
	}

	public NodeData(int nmax) {
		//m_n = nmax;
		//cerr<<"YYY\t\tallocating\t"<<this<<endl;
		//m_children = new Vector<NodeData>();
		m_Fb = new FMatrix();
		m_Ft = new FMatrix();
		resize(nmax);
//		m_nTaxonID = 0;
	}

	public void resize(int nmax) {
			//m_n = nmax;
			if (nmax>=1) {
				//m_Nt = new double[m_n+1];
				//m_Nb = new double[m_n+1];
				m_Ft.resize(nmax);
				m_Fb.resize(nmax);
			}
		}

	public double t() {return getLength();}
//	public void set_t(double t) {m_fLength = t;}

//	public double coalescenceRate() {
//		return m_fCoalescenceRate;
//	}
//	public void set_coalescenceRate(double fCoalescenceRate) {m_fCoalescenceRate = fCoalescenceRate;}

	public NodeData copyx() throws CloneNotSupportedException {
		NodeData node = new NodeData();
		node.m_fHeight = m_fHeight;
		node.m_iLabel = m_iLabel;
		node.m_sMetaData = m_sMetaData;
		node.setParent(null);
		if (getLeft() != null) {
			node.setLeft(((NodeData)getLeft()).copyx());
			node.setRight(((NodeData)getRight()).copyx());
			node.getLeft().setParent(node);
			node.getRight().setParent(node);
		}
		//node.m_n = m_n;
		return node;
	}

	public NodeData copy() {
		NodeData node = new NodeData();
		node.m_fHeight = m_fHeight;
		node.m_iLabel = m_iLabel;
		node.m_sMetaData = m_sMetaData;

//		node.set_coalescenceRate(m_fCoalescenceRate);
		//node.m_n = m_n;
		//if (m_Nt != null) {
			//node.m_Nt = new double[m_Nt.length];
			//System.arraycopy(m_Nt, 0, node.m_Nt, 0, m_Nt.length);
			//node.m_Nb = new double[m_Nb.length];
			//System.arraycopy(m_Nb, 0, node.m_Nb, 0, m_Nb.length);
		//}
//		node.Ft = new FMatrix(Ft);
//		node.Fb = new FMatrix(Fb);

		node.setParent(null);
		if (getLeft() != null) {
			NodeData left = ((NodeData)getLeft()).copy(); 
			node.setLeft(left);
			left.setParent(node);
			if (getRight() != null) {
				NodeData right = ((NodeData)getRight()).copy();
				node.setRight(right);
				right.setParent(node);
			}
		}
		return node;
	}

	public String getNewickMetaData() {
		return "";//"[theta=" + 2.0/coalescenceRate() + ']';
	}

	public void resizeF(int n) {
		m_Fb.resize(n);
		m_Ft.resize(n);
	} // resizeF
	public FMatrix cloneFbx() {
		return new FMatrix(m_Fb);
	} // cloneF
	public FMatrix cloneFtx() {
		return new FMatrix(m_Ft);
	} // cloneF
	public void assignFb(FMatrix _F) {
		m_Fb.assign(_F);
	} // cloneF
	public void assignFt(FMatrix _F) {
		m_Ft.assign(_F);
	} // cloneF

	int m_nCacheIDB;
	public int getCacheIDB() {
		return m_nCacheIDB;
	}
	public void setCacheIDB(int cacheIDB) {
		m_nCacheIDB = cacheIDB;
	}

	int m_nCacheIDT;
	public int getCacheIDT() {
		return m_nCacheIDT;
	}
	public void setCacheIDT(int cacheIDT) {
		m_nCacheIDT = cacheIDT;
	}
	
	/** used for lazy updating **/
	//boolean m_bIsDirty = true;
	@Override
	public void setMetaData(String sPattern, Object fValue) {
//		if (sPattern.equals("coalescenceRate")) {
//			m_fCoalescenceRate = (Double) fValue;
//		}
//		if (sPattern.equals("theta")) {
//			m_fCoalescenceRate = 2.0/(Double) fValue;
//		}
		//super.setMetaData(sPattern, fValue);
	}
	
	@Override
	public Object getMetaData(String sPattern) {
		if (m_sMetaData != null && m_sMetaData.indexOf(sPattern+"=")>=0) {
			int i = m_sMetaData.indexOf(sPattern+"=") + sPattern.length() + 1;
			String sStr = m_sMetaData.substring(i, m_sMetaData.length());
			if (sStr.indexOf(',')>=0) {
				sStr = sStr.substring(0, sStr.indexOf(','));
			}
			double f = Double.parseDouble(sStr);
			return f;
		}
		
		
//		if (sPattern.equals("coalescenceRate")) {
//			return m_fCoalescenceRate;
//		}
//		if (sPattern.equals("theta")) {
//			 return 2.0/m_fCoalescenceRate;
//		}
		return super.getMetaData(sPattern);
	}

//	@Override
//	public void assignFromFragile(Node other) {
//		m_fGamma = ((NodeData)other).m_fGamma;
//	}

	public String toString() {
		return toNewick(null);
	}
}

