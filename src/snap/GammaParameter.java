
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


import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("Represents population size associated with each node in a tree.")
public class GammaParameter extends RealParameter {
	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree associated with this (array) parameter", Validate.REQUIRED);
	public Input<Boolean> m_bInitFromTree = new Input<Boolean>("initFromTree", "whether to initialize from starting tree values (if true), or vice versa (if false)");
	public Input<String> m_pPattern = new Input<String>("pattern", "pattern of metadata element associated with this parameter in the tree");
	
	public GammaParameter() {
	}
	
	// do not try this at home: this is here as long as variable Traits are not 
	// properly handled in Beast II 
	static Tree m_tree;
	
	@Override
	public void initAndValidate() throws Exception {
    	if (lowerValueInput.get() != null) {
    		m_fLower = lowerValueInput.get();
    	} else {
    		m_fLower = Double.NEGATIVE_INFINITY;
    	}
    	if (upperValueInput.get() != null) {
    		m_fUpper = upperValueInput.get();
    	} else {
    		m_fUpper = Double.POSITIVE_INFINITY;
    	}

    	m_tree = m_pTree.get();
		
		if (m_bInitFromTree.get() == true) {
			values = new Double[m_tree.getNodeCount()];
			m_tree.getMetaData(m_tree.getRoot(), values, m_pPattern.get());
			dimensionInput.setValue(new Integer(values.length), this);
		} else {
			values = new Double[m_tree.getNodeCount()];
			dimensionInput.setValue(new Integer(values.length), this);

	    	List<Double> sValues = valuesInput.get();
	        for (int i = 0; i < values.length; i++) {
	            values[i] = new Double(sValues.get(i % sValues.size()));
	        }
			m_tree.setMetaData(m_tree.getRoot(), values, m_pPattern.get());
		}
    	m_bIsDirty = new boolean[dimensionInput.get()];
    	m_pTree.setValue(null, this);
	}
	
	@Override
    public GammaParameter copy() {
    	GammaParameter copy = new GammaParameter();
    	copy.ID = ID;
    	copy.values = new Double[values.length];
    	System.arraycopy(values, 0, copy.values, 0, values.length);
        copy.m_bIsDirty = new boolean[values.length];
    	copy.m_fLower = m_fLower;
    	copy.m_fUpper = m_fUpper;
    	//TODO: verify index is set properly
    	copy.index = index;
    	copy.m_pTree = m_pTree;
    	copy.m_pPattern = m_pPattern;
    	return copy;
    }
    
    public void prepare() {
		//syncTree(m_pTree.get().getRoot(), values, m_pPattern.get());
    	Tree tree = (Tree)m_tree.getCurrent();
    	Node node = tree.getRoot();
		syncTree(node, values, m_pPattern.get());
	}
    
	void syncTree(Node node, Double [] fValues, String sPattern) {
		node.setMetaData(sPattern, fValues[Math.abs(node.getNr())]);
		if (!node.isLeaf()) {
			syncTree(node.getLeft(), fValues, sPattern);
			if (node.getRight() != null) {
				syncTree(node.getRight(), fValues, sPattern);
			}
		}
	}

//    public String toString() {
//    	StringBuffer buf = new StringBuffer();
//    	buf.append(m_sID);
//    	buf.append(": ");
//    	for (int i = 0; i < values.length; i++) {
//    		buf.append(values[i] + " ");
//    	}
//		if (m_pTree != null) {
//			buf.append("Associated with " + m_pPattern.get() + " for tree '" + m_tree.getID() + "'");
//		}
//    	return buf.toString();
//    }
}
