
/*
 * File NodeSwapper.java
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
package snap.operators;



import java.util.Vector;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;



@Description("Randomly selects two nodes and swap them in the tree.")
public class NodeSwapper extends Operator {
	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations");
	
	int m_nNodeCount = -1;
	
	@Override
	public void initAndValidate() {
		m_nNodeCount = m_pTree.get().getNodeCount();
	}

	@Override
	public double proposal() { // throws Exception {
		double hastingsRatio = 1.0;
		Tree tree = m_pTree.get(this);
		Node [] nodes = tree.getNodesAsArray();

		//First select a triple (x,y,m) where x and y are a random pair of leaves and m is the mrca of x and y.


		//choose two random leaves.
		int ntax = m_nNodeCount/2+1;
		int xid = Randomizer.nextInt(ntax);
		int yid = Randomizer.nextInt(ntax);

		Node x = nodes[xid];
		Node y = nodes[yid];
		
		//Check that x and y have different parents
		if (x.getParent().getNr() == y.getParent().getNr())
			return Double.NEGATIVE_INFINITY; //No move.

		//Find the MRCA
		Vector<Node> ancerstorsx = new Vector<Node>();
		Node xa = x;
		while (!xa.isRoot()) {
			ancerstorsx.add(xa);
			xa = xa.getParent();
		}
		ancerstorsx.add(xa);
		Vector<Node> ancerstorsy = new Vector<Node>();
		Node ya = y;
		while (!ya.isRoot()) {
			ancerstorsy.add(ya);
			ya = ya.getParent();
		}
		ancerstorsy.add(ya);
		Node mrca = null;
		int nSizeX = ancerstorsx.size();
		int nSizeY = ancerstorsy.size();
		int i = 0;
		mrca = ancerstorsx.elementAt(nSizeX - 1);
		while (ancerstorsx.elementAt(nSizeX - i - 1).getNr() == ancerstorsy.elementAt(nSizeY - i - 1).getNr() && i < nSizeX && i < nSizeY) {
			mrca = ancerstorsx.elementAt(nSizeX - i - 1);
			i++;
		}
	
//		throw new Exception("nspecies field is increasing up the tree" is no longer true!!! FIX THIS);
//		Node xa=x, ya=y;
//		double mheight = 0.0; //length of path up to mrca
//		while(xa.getNr() != ya.getNr()) {
//			if (xa.getNr() <= ya.getNr()) {
//				mheight += xa.getLength();
//				xa=xa.getParent();
//			} else {
//				ya = ya.getParent();
//			}
//		}
//		Node mrca = xa;

		
		/*
		 Let mheight be the height of the MRCA, xheight & yheight be the heights of x and y.
		 Let minb be the length of the shortest branch adjancent to, and below, the mrca.
		 
		 We want to cut the tree at a height U[max(xheight,yheight), mheight-minb]
		 
		 */
		
		
		
		
		double mheight = mrca.getHeight();
		double xheight = x.getHeight();
		double yheight = y.getHeight();
		double minb = Math.min(mrca.getLeft().getLength(), mrca.getRight().getLength());
		double cutheight = Randomizer.nextDouble()*(mheight - minb - Math.max(xheight,yheight)) + Math.max(xheight,yheight);
		
		
		//System.err.println("mheight = "+mheight+" "+xheight+" "+yheight+" "+minb+" "+cutheight);
		
		/* Identify the nodes on the paths from x to y that are directly below the cut point */
		xa = x;
		while(xheight+(xa.getLength())<cutheight) {
			xheight+=(xa.getLength());
			xa = xa.getParent();
		}
		ya=y;
		while(yheight+(ya.getLength())<cutheight) {
			yheight+=ya.getLength();
			ya = ya.getParent();
		}
		
		/* Compute the lengths of the branches after the swap. This is (segment below cut) + (segment above cut) */
		double xlength = (cutheight - xheight) + (yheight + ya.getLength() - cutheight);
		double ylength = (cutheight - yheight) + (xheight + xa.getLength() - cutheight);

		/* Perform the swap and fix heights */
		swap(nodes, xa, ya);
		
		xa.setHeight(xa.getParent().getHeight() - xlength);
		ya.setHeight(ya.getParent().getHeight() - ylength);
		
		
		return Math.log(hastingsRatio);
		
	}
	
	void swap(Node [] nodes, Node x, Node y) {
		int ix = x.getNr();
		int iy = y.getNr();
		int ixp = x.getParent().getNr();
		int iyp = y.getParent().getNr();
		if (nodes[ixp].getLeft().getNr() == ix) {
			nodes[ixp].setLeft(y);
		} else {
			nodes[ixp].setRight(y);
		}
		if (nodes[iyp].getLeft().getNr() == iy) {
			nodes[iyp].setLeft(x);
		} else {
			nodes[iyp].setRight(x);
		}
		y.setParent(nodes[ixp]);
		x.setParent(nodes[iyp]);
	}

} // class NodeSwapper
