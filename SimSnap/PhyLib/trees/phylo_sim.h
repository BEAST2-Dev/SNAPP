/**
 * phylo_sim.h
 * 
 * Utilities for generating random phylogenies.
 * 
 */

#ifndef PHYLO_SIM_H_
#define PHYLO_SIM_H_

#include"../global/stdIncludes.h"
#include"../utilities/phylibException.h"
#include "phylo.h"
#include"../utilities/random.h"

namespace Phylib {
	template<typename T> phylo<T> sim_uniform_phylo(int ntax) {
		if (ntax<=0)
			throw PhylibException("Simulating tree with 0 or negative taxa!");
		phylo<T> tree;
		
		//Initial tree has one leaf.
		tree.insert_child(tree.header());
		tree.root()->id = 0;
		int nnodes = 1;
		for (int i = 1; i<ntax; i++) {
			typename phylo<T>::iterator p = tree.root(); //position to insert next leaf
			int pos = std::rand()%nnodes; //Number of steps from root for this position
			for (int j=0; j<pos; j++) {
				p = p.next_pre();
			}
			//Insert node above p.
			tree.insert_sibling(p); //add sibling
			tree.graft_child(p.right(), tree, p); //move node to below sibling
			//Insert new leaf as sibling of p.
			tree.insert_sibling(p);
			p.right()->id = i;
			nnodes+=2;
		}
		return tree;
	}
	
	template<typename T> phylo<T> sim_caterpillar_phylo(int ntax) {
		
		typedef typename phylo<T>::iterator ITERATOR;
		
		if (ntax<=0)
			throw PhylibException("Simulating tree with 0 or negative taxa!");
		phylo<T> tree;
		
		vector<int> leaves(ntax);
		for(int i=0;i<ntax;i++)
			leaves[i] = i;
		//std::random_shuffle(leaves.begin(),leaves.end());
		//Initial tree has one leaf.
		tree.insert_child(tree.header());
		tree.root()->id = leaves[0];
		ITERATOR bottom = tree.root(); //Bottom root in the tree.
		
		for(int i=1;i<ntax;i++) {
			//Insert node below bottom
			ITERATOR newnode1 = tree.insert_child(bottom);
			ITERATOR newnode2 = tree.insert_sibling(newnode1);
			newnode1->id = bottom->id;
			bottom->id = -1;
			newnode2->id = leaves[i];
			bottom = newnode2;
		}
		return tree;
	}
	
	template<typename T> phylo<T> _sim_balanced_phylo_recurse(int ntax) {
		if (ntax == 1) {
			phylo<T> tree;
			tree.insert_child(tree.header());
			tree.root()->id = 1;
			return tree;
		}
		else {
			int left = (int) std::ceil(ntax/2);
			int right = ntax - left;
			phylo< T> treeL = _sim_balanced_phylo_recurse<T>(left);
			phylo< T> treeR = _sim_balanced_phylo_recurse<T>(right);
			phylo< T> tree;
			tree.insert_child(tree.header()); //Root
			tree.graft_child(tree.root(),treeR);
			tree.graft_child(tree.root(),treeL);
			return tree;
		}
	}
	
	template<typename T> phylo<T> sim_rooted_balanced(int ntax) {
		
		phylo<T> tree = _sim_balanced_phylo_recurse<T>(ntax);
		vector<int> taxaIds(ntax);
		for(int i=0;i<ntax;i++)
			taxaIds[i] = i;
		std::random_shuffle(taxaIds.begin(),taxaIds.end());
		
		int i=0;
		for(typename phylo<T>::iterator p = tree.root(); !p.null(); p = p.next_pre()) {
			if (p.leaf())
				p->id = taxaIds[i++];
		}
		return tree;
		
	}
	
	template<typename T> phylo<T> sim_unrooted_balanced(int ntax) {
		if (ntax<3)
			return sim_rooted_balanced<T>(ntax);
		
		vector<int> taxaIds(ntax);
		for(int i=0;i<ntax;i++)
			taxaIds[i] = i;
		std::random_shuffle(taxaIds.begin(),taxaIds.end());
		
		phylo<T> tree;
		typename phylo<T>::iterator r = tree.insert_child(tree.header());
		phylo<T> tree1, tree2, tree3;
		int n3 = ntax/3;
		int n2 = (ntax - n3)/2;
		int n1 = ntax - n2 - n3;
		
		tree1 = sim_rooted_balanced<T>(n1);
		tree2 = sim_rooted_balanced<T>(n2);
		tree3 = sim_rooted_balanced<T>(n3);
		tree.graft_child(r,tree3);
		tree.graft_child(r,tree2);
		tree.graft_child(r,tree1);
		
		int i=0;
		for(typename phylo<T>::iterator p = tree.root(); !p.null(); p = p.next_pre()) {
			if (p.leaf())
				p->id = taxaIds[i++];
		}
		
		return tree;
	}
	
		
	

}

#endif /*PHYLO_SIM_H_*/




