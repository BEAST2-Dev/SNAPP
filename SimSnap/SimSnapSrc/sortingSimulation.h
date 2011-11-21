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

extern string g_simtree;

/**
 Node type for gene trees 
 **/
class geneTreeNode : public basic_newick {
public: 
	int state;
	int species;
	double height;
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
	int numberCoalescences;
	int numberCoalescencesTotal;

	
	simNodeData() : basic_newick(), gamma(1.0) , numberCoalescences(0) {}
	simNodeData(int id, double len, double gamma_val): basic_newick(id, len), gamma(gamma_val), numberCoalescences(0) {}
	simNodeData(const basic_newick& data) {
		basic_newick::copy(data);
		
		
		size_t equalPos = data.meta_data.find_first_of("=");
		if (equalPos!=string::npos) {
			if (data.meta_data.compare(0,5,"theta") == 0)
				theta = atof((data.meta_data.substr(equalPos+1)).c_str());
			else if (data.meta_data.compare(0,15,"coalescenceRate")==0)
				theta = 2.0 / atof((data.meta_data.substr(equalPos+1)).c_str());
			else	{		
				theta = atof((data.meta_data.substr(equalPos+1)).c_str());
				cerr<<"Error reading in tree"<<endl;
			}
		}
		else
			theta = atof(data.meta_data.c_str());
		
		//TODO: Make this more robust.
		
		gamma = 0.0;
		numberCoalescences = 0;
	}
		
	void copy(const simNodeData& data) {
		basic_newick::copy(data);
		gamma = data.gamma;
		theta = data.theta;
		numberCoalescences = data.numberCoalescences;
		
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
phylo<geneTreeNode> simulateGeneTree(phylo<simNodeData>& speciesTree, const vector<uint>& sampleSizes);

/**
 Simulate a single SNP. if rejectConstant = true, a non-constant site (polymorphism) will be generated.
 **/
void simulateSingleSite(phylo<simNodeData>& speciesTree, double u, double v, const vector<uint>& sampleSizes, vector<uint>& redCounts, bool rejectConstant, bool onlyRootMutation, bool outputTree = false, int site=0);
void simulateSingleSite(phylo<simNodeData>& speciesTree, double u, double v, const vector<uint>& sampleSizes, vector<uint>& redCounts, bool rejectConstant, bool onlyRootMutation, uint& numberAttempts, bool outputTree = false, int site=0);
/**
 Simulate multiple unlinked sites.
 
 If rejectConstant is true, only non Constant sites are returned. We use a simple rejection algorithm (any better ideas?)
 **/

void simulateMultipleSites(phylo<basic_newick>& tree, double u, double v, const vector<uint>& sampleSizes, int nSites, bool rejectConstant, bool onlyRootMutation, vector<vector<uint> >& redCounts, bool outputTrees = false);
void simulateMultipleSites(phylo<basic_newick>& tree, double u, double v, const vector<uint>& sampleSizes, int nSites, bool rejectConstant, bool onlyRootMutation, vector<vector<uint> >& redCounts, double& proportionConstant, bool outputTrees = false);

#endif


