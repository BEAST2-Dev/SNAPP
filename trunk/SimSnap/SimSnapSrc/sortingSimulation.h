/*
 *  simulation.h
 *  SingleSiteSorter
 *
 *  Created by David Bryant on 29/08/08.
 *  Copyright 2008 David Bryant
 *
 */

#ifndef SORTING_SIMULATION_H
#define SORTING_SIMULATION_H

#include "phylib.h"

using namespace Phylib;


/**
 Node type for gene trees 
 **/
class geneTreeNode : public basic_newick {
public: 
	int state;
	int species;
	double P[2][2]; //Transition probabilities.
};

/**
 Node type for species trees
 **/
class simNodeData : public basic_newick {
public:
	double theta;
	double gamma;	
	list< phylo<geneTreeNode> > lineages;
	
	
	simNodeData() : basic_newick(), gamma(1.0) {}
	simNodeData(int id, double len, double gamma_val): basic_newick(id, len), gamma(gamma_val) {}
	simNodeData(const basic_newick& data) {
		basic_newick::copy(data);
		theta = atof(data.meta_data.c_str());
		gamma = 0.0;
	}
		
	void copy(const simNodeData& data) {
		basic_newick::copy(data);
		gamma = data.gamma;
		theta = data.theta;
		lineages.clear();
		lineages.insert(lineages.begin(),data.lineages.begin(),data.lineages.end());
	}
	
	simNodeData& operator=(const simNodeData& data) {
		if (this!=&data) 
			copy(data);
		return *this;
	}
	
};

/**
 Copies the tree structure and initialises the fields for gene tree simulation.
 **/
phylo<simNodeData> initialiseSimTree(phylo<basic_newick>& tree);

/**
 Simulates a gene tree. 
 **/
phylo<geneTreeNode> simulateGeneTree(phylo<simNodeData>& speciesTree, const vector<uint> sampleSizes);

/**
 Simulate a single SNP. if rejectConstant = true, a non-constant site (polymorphism) will be generated.
 **/
void simulateSingleSite(phylo<simNodeData>& speciesTree, double u, double v, const vector<uint>& sampleSizes, vector<uint>& redCounts, bool rejectConstant, bool outputTree = false);
void simulateSingleSite(phylo<simNodeData>& speciesTree, double u, double v, const vector<uint>& sampleSizes, vector<uint>& redCounts, bool rejectConstant, uint& numberAttempts, bool outputTree = false);
/**
 Simulate multiple unlinked sites.
 
 If rejectConstant is true, only non Constant sites are returned. We use a simple rejection algorithm (any better ideas?)
 **/

void simulateMultipleSites(phylo<basic_newick>& tree, double u, double v, const vector<uint>& sampleSizes, int nSites, bool rejectConstant, vector<vector<uint> >& redCounts, bool outputTrees = false);
void simulateMultipleSites(phylo<basic_newick>& tree, double u, double v, const vector<uint>& sampleSizes, int nSites, bool rejectConstant, vector<vector<uint> >& redCounts, double& proportionConstant, bool outputTrees = false);

#endif


