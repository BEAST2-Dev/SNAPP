/*
 *  simulation.cpp
 *  SingleSiteSorter
 *
 *  Created by David Bryant on 29/08/08.
 *  Copyright 2008 David Bryant
 *
 */

#include "sortingSimulation.h"

using namespace Phylib;
/**
 Copies the tree structure and initialises the fields for gene tree simulation.
 
 rate is the mutation rate, needed to convert theta values into gamma values
 **/
phylo<simNodeData> initialiseSimTree(phylo<basic_newick>& tree, double rate) {
	typedef phylo<simNodeData>::iterator S_ITER;

	
	phylo<simNodeData> simTree;
	Phylib::copy(tree,simTree);
	
	//Convert the theta values into gamma values
	// The expected divergence between two individuals is 
	// theta = 2*rate/gamma.
	for(S_ITER s = simTree.leftmost_leaf();!s.null(); s = s.next_post()) 
		s->gamma = 2.0*rate/(s->theta);
	
	
	return simTree;
}


/**
 Simulates a gene tree. 
 We use GIllespie's algorithm up each branch of the species tree to simulate the coalescent.
 **/
phylo<geneTreeNode> simulateGeneTree(phylo<simNodeData>& speciesTree, const vector<uint> sampleSizes) {
	typedef phylo<simNodeData>::iterator S_ITER;
	typedef phylo<geneTreeNode>::iterator G_ITER;
	typedef list<phylo<geneTreeNode> >::iterator L_ITER;
	
	uint nleaves = 0;
	
	for(S_ITER s = speciesTree.leftmost_leaf(); !s.null(); s = s.next_post()) {
		if (s.leaf()) {
			//Create new lineages at the base of the branch, one for each individual in the sample.
			int count=sampleSizes[s->id];
			s->lineages.clear();
			for(int i=0;i<count;i++) {
				phylo<geneTreeNode> leafNode;
				G_ITER newLeaf = leafNode.insert_child(leafNode.header());
				newLeaf->id = nleaves;
				newLeaf->length = 0.0;
				newLeaf->species = s->id;
				nleaves++;
				s->lineages.push_back(leafNode);
			}
		} else {
			//Splice the lineages from the children to the bottom of this branch.
			s->lineages.clear();
			for(S_ITER c = s.left(); !c.null(); c=c.right()) 
				s->lineages.splice(s->lineages.end(),c->lineages);
		}
		
		
		//Evolve the lineages up the branch.
		double height_in_branch = 0.0;
		int k = s->lineages.size();
		double branch_length = s->length;
		
		for( ; ; ) {
			if (k==1) { //All lineages coalesced
				if(s.root()) {//If we're at the root of the species tree, and coalesced to one lineage, we're done.
					phylo<geneTreeNode> geneTree;
					geneTree.swap(s->lineages.front());
					s->lineages.clear();
					return geneTree;
				}
				else {
					//Skip straight to the top of the branch
					for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) 
						(*node).root()->length+=branch_length - height_in_branch;
					break;
				}	
			}
			
			//Waiting time until next coalescent event.
			double wait = random_exp(2.0 / ( (double)k*(k-1.0) * s->gamma));
			
			if (!s.root() && height_in_branch+wait>=branch_length) {
				//Reached the top of the branch
				for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) 
					(*node).root()->length+=branch_length - height_in_branch;
				break;
			}
			
			for(L_ITER node = s->lineages.begin();node!=s->lineages.end();node++) 
				(*node).root()->length+=wait;
			height_in_branch+=wait;
			
			//Choose a random pair (a,b),  a<b, to coalesce.
			int a = random_num(k);
			int b = random_num(k-1);
			if (b>=a) 
				b++;
			else { //Swap a and b so that a<b.
				int c=a; a=b; b = c;		
			}	
			
			//Update the gene tree and list of lineages.
			phylo<geneTreeNode> newNode;
			newNode.insert_child(newNode.header());
			newNode.root()->length = 0.0;
			
			L_ITER ptr = s->lineages.begin();
			std::advance(ptr,a);
			phylo<geneTreeNode>& child_a = *ptr;
			std::advance(ptr,b-a);
			phylo<geneTreeNode>& child_b = *ptr;
			
			newNode.graft_child(newNode.root(),child_a);
			newNode.graft_child(newNode.root(),child_b);
			ptr = s->lineages.begin();
			advance(ptr,b);
			s->lineages.erase(ptr);
			ptr = s->lineages.begin();
			advance(ptr,a);
			s->lineages.erase(ptr);
			s->lineages.push_back(newNode);
			k--;
			
			//TODO: Clean up this list pointer rubbish.
		}
	}
	throw PhylibException("Error in gene tree simulation");
	phylo<geneTreeNode> nullTree;
	return nullTree;
}

typedef phylo<geneTreeNode>::iterator G_ITER;

/**
 Pre-compute the mutation probabilities using the 2-state model
 **/
static void computeTransitionProbs(phylo<geneTreeNode>& G, double u, double v) {
	double pi_0 = v/(u+v), pi_1 = 1.0-pi_0;
	
	for(G_ITER p = G.leftmost_leaf();p!=G.root();p=p.next_post()) {
		double t = p->length;
		double x = 1.0 - std::exp(-(u+v)*t);
		p->P[1][0] = pi_0*x;
		p->P[1][1] = 1.0 - pi_0*x;
		p->P[0][1] = pi_1*x;
		p->P[0][0] = 1.0 - pi_1*x;
	}
}

/**
 Simulate a two state character on the tree
 **/
static void simulateCharacter(phylo<geneTreeNode>& G, double pi_0) {
	G.root()->state = (randu()<pi_0)?0:1;
	
	for(G_ITER p = G.root().next_pre();!p.null();p=p.next_pre()) {
		uint i = p.par()->state;
		p->state = (randu()<p->P[i][0])?0:1;
	}
}

/**
 Check to see whether a character (given by allele counts) is constant or not.
 **/
static bool checkIfConstant(const vector<uint>& redCount, const vector<uint>& sampleSizes) {
	bool allGreen, allRed;
	allGreen = allRed = true;
	for(uint i=0;i<redCount.size();i++) {
		uint r = redCount[i];
		if (r>0)
			allGreen = false;
		if (r<sampleSizes[i])
			allRed = false;
		if (!allGreen && !allRed)
			break;
	}
	//cerr<<((allRed||allGreen)?"constant":"variable")<<endl;
	return (allRed || allGreen);
}

/**
 Wrapper for simulateSingleSite without the proportionConstant variable
 **/
void simulateSingleSite(phylo<simNodeData>& speciesTree, double u, double v, const vector<uint>& sampleSizes, vector<uint>& redCounts, bool rejectConstant, bool outputTree) {
	uint numberAttempts;
	simulateSingleSite(speciesTree,u,v,sampleSizes,redCounts, rejectConstant, numberAttempts, outputTree);
}
	
/**
 Simulate a gene tree and then simulate a single binary character on it. Returns allele counts for each species.
 Uses a rejection algorithm to simulate non-constant characters. The numberAttempts is the number of characters we had to simulate to a polymorphic one.
 **/
 void simulateSingleSite(phylo<simNodeData>& speciesTree, double u, double v, const vector<uint>& sampleSizes, vector<uint>& redCounts, bool rejectConstant, uint& numberAttempts, bool outputTree) {
	
			
	//Now evolve the markers.
	double pi_0 = v/(u+v);
	uint nSpecies = sampleSizes.size();
	numberAttempts = 0;
	 
	bool allConstant = true;

	do {
		numberAttempts++;
		//First generate the species tree
		phylo<geneTreeNode> G = simulateGeneTree(speciesTree, sampleSizes);
		
		
		computeTransitionProbs(G,u,v);
		
		simulateCharacter(G, pi_0);
		redCounts.resize(nSpecies);
		std::fill(redCounts.begin(),redCounts.end(),0);
		for(G_ITER p = G.root();!p.null();p=p.next_pre()) {
			if (p.leaf() && p->state==0) 
				redCounts[p->species]++;
		}
		
		if (rejectConstant) 
			allConstant = checkIfConstant(redCounts, sampleSizes);	
		
		if (outputTree && (!rejectConstant || !allConstant) ) {
			print_newick(cout, G, true);
			cout<<endl;
		}
		
		G.clear();

	} while(allConstant && rejectConstant);	
}

/**
 Simulate multiple sites.
 
 For each site, simulates a gene tree according to the multi-species coalescent, then simulates a single binary character on that gene tree.
 If rejectConstant is set to true, the process will repeat at each site until a non-constant character is obtained. This could cause problems
 if the mutation rate is really low.
 **/

/**wrapper**/
void simulateMultipleSites(phylo<basic_newick>& tree, double u, double v, const vector<uint>& sampleSizes, int nSites, bool rejectConstant, vector<vector<uint> >& redCounts, bool outputTrees) {
	double proportionConstant;
	simulateMultipleSites(tree,u,v,sampleSizes,nSites,rejectConstant,redCounts,proportionConstant, outputTrees);
}


void simulateMultipleSites(phylo<basic_newick>& tree, double u, double v, const vector<uint>& sampleSizes, int nSites, bool rejectConstant, vector<vector<uint> >& redCounts, double& proportionConstant, bool outputTrees) {
	redCounts.resize(nSites);
	double rate = 2.0*u*v/(u+v);
	
	phylo<simNodeData> simTree = initialiseSimTree(tree,rate);
	
	cout<<"SimTree = ";
	print_newick(cout,simTree,true,true);
	cout<<endl;
	
	uint numberAttempts, numberAttemptsTotal;
	numberAttemptsTotal = 0;
	
	for(int i=0;i<nSites;i++) {
		simulateSingleSite(simTree, u, v, sampleSizes, redCounts[i], rejectConstant,numberAttempts, outputTrees);
		numberAttemptsTotal+=numberAttempts;
	}
	
	simTree.clear();
	if (rejectConstant)
		proportionConstant = 1.0-(double)nSites/(double)numberAttemptsTotal;
	else
		proportionConstant = -1.0;
}

















