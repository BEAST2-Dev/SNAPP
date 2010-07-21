
/*
 * File GammaParameter.java
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


import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

@Description("Represents population size associated with each node in a tree.")
public class GammaParameter extends RealParameter {
	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree associated with this (array) parameter");
	public Input<Boolean> m_bInitFromTree = new Input<Boolean>("initFromTree", "whether to initialize from starting tree values (if true), or vice versa (if false)");
	public Input<String> m_pPattern = new Input<String>("pattern", "pattern of metadata element associated with this parameter in the tree");
	
	public GammaParameter() {
	}


	@Override
	public void initAndValidate() throws Exception {
		Tree tree = m_pTree.get();
		
		if (m_bInitFromTree.get() == true) {
			values = new Double[tree.getNodeCount()];
			tree.getMetaData(tree.getRoot(), values, m_pPattern.get());
			m_nDimension.setValue(new Integer(values.length), this);
		} else {
			values = new Double[tree.getNodeCount()];
			m_nDimension.setValue(new Integer(values.length), this);
			for (int i = 0; i < values.length; i++) {
				values[i] = (Double) m_pValues.get();
			}
			tree.setMetaData(tree.getRoot(), values, m_pPattern.get());
		}
    	m_bIsDirty = new boolean[m_nDimension.get()];
	}
	
	@Override
    public GammaParameter copy() {
    	GammaParameter copy = new GammaParameter();
    	copy.m_sID = m_sID;
    	copy.values = new Double[values.length];
    	System.arraycopy(values, 0, copy.values, 0, values.length);
        copy.m_bIsDirty = new boolean[values.length];
    	copy.m_fLower = m_fLower;
    	copy.m_fUpper = m_fUpper;
    	//TODO: verify index is set properly
    	//copy.index = getIndex();
    	copy.m_pTree = m_pTree;
    	copy.m_pPattern = m_pPattern;
    	return copy;
    }
    
    public void prepare() {
		syncTree(m_pTree.get().getRoot(), values, m_pPattern.get());
	}
    
	void syncTree(Node node, Double [] fValues, String sPattern) {
		node.setMetaData(sPattern, fValues[Math.abs(node.getNr())]);
		if (!node.isLeaf()) {
			syncTree(node.m_left, fValues, sPattern);
			syncTree(node.m_right, fValues, sPattern);
		}
	}

    public String toString() {
    	StringBuffer buf = new StringBuffer();
    	buf.append(m_sID);
    	buf.append(": ");
    	for (int i = 0; i < values.length; i++) {
    		buf.append(values[i] + " ");
    	}
		if (m_pTree != null) {
			buf.append("Associated with " + m_pPattern.get() + " for tree '" + m_pTree.get().getID() + "'");
		}
    	return buf.toString();
    }
}
