/*
 *  simSnap.cpp
 *  SimSnap
 *
 *  Created by David Bryant on 8/03/10.
 *  Copyright 2010 University of Auckland. All rights reserved.
 *
 */


//TODO: extend tree parser so that [theta=xxx] is allowed.

#include "sortingSimulation.h"
#include "characterData.h"
//#include "posteriorCheck.h"


void printUsage(ostream& os) {
	os<<"SimSnap\n\nSimulates SNPs on a species tree.\n";
	os<<"Usage:\n\n\tSimSnap [-nc] [nsites] filename\n\n";
	os<<"\tFlags are:\n";
	os<<"\t-n \tOutput nexus file (default is to output snap .xml format)\n";
	os<<"\t-c \tInclude constant sites (default is to simulate only polymorphic sites)\n";
	os<<"INPUT FILE FORMAT:\n\n";
	os<<"<number of species>\n";
	os<<"<species1-name>\t<sample size species 1>\n";
	os<<"<species2-name>\t<sample size species 2>\n";
	os<<"...\n";
	os<<"<species_n-name>\t<sample size species n>\n";
	os<<"<mutation rate u> <mutation rate v>\n";
	os<<"<number of trees>\n";
	os<<"Trees in newick notation, with branch lengths. Theta values indicated in square brackets following node.\n\n";
	os<<"EXAMPLE\n\n";
	os<<"3\n";
	os<<"A 3\nB 4\nC 2\n";
	os<<"0.2 0.6\n2\n";
	os<<"((A[0.1]:0.3,B[0.3]:0.3)[0.2]:0.1,C[0.2]:0.4)[0.3]"<<endl;
	os<<"((C[0.1]:1.3,B[0.3]:1.3)[0.2]:0.1,A[0.2]:1.4)[0.3]\n"<<endl;
	exit(1);
}

class ArgumentParser {
public: 
	bool outputXML;
	bool excludeConst;
	bool outputTrees;
	string inputfile;
	int nsites;
	ArgumentParser(int argc, char* argv[]) {

		outputXML=true;
		excludeConst = true;
		outputTrees = false;
		
		
		nsites = 0;
		inputfile = "";
		
		if (argc<3)
			printUsage(cerr);
		
		int arg = 1;
		
		//First read in the flags.
		string flags = string(argv[arg]);
		if (flags[0]=='-') {
			if (flags.find_first_not_of("-nct")!=string::npos)
				printUsage(cerr);
			
			if (flags.find_first_of('n')!=string::npos)
				outputXML = false;
			if (flags.find_first_of('c')!=string::npos)
				excludeConst = false;
			if (flags.find_first_of('t')!=string::npos)
				outputTrees = true;
			
			arg++;
			if (argc < 4)
				printUsage(cerr);
		}
		
		nsites = atoi(argv[arg]); 
		arg++;
		inputfile = string(argv[arg]);
	}
};
		
void output_xml(ostream& os, const vector<string>& taxa, const vector<uint>& sampleSizes, double u, double v, const vector<vector<uint> >&alleleCounts, const string fileroot="test") {
	os<<"<snap version='2.0' namespace='snap:beast.util'>\n";
	os<<"\n";
	os<<"<map name='snapprior'>snap.likelihood.SnAPPrior</map>\n";
	os<<"<map name='snaptreelikelihood'>snap.likelihood.SnAPTreeLikelihood</map>\n";

	os<<"\n\n\n\t<!-- n = "<<alleleCounts.size()<<" -->\n";
	//Compute largest sample size.
	int statecount = *max_element(sampleSizes.begin(),sampleSizes.end());
	os<<"\t<data id='alignment' dataType='integerdata' statecount='"<<statecount + 1<<"'>\n";
	for(uint i=0;i<taxa.size();i++) {
		os<<"\t\t<sequence taxon='"<<taxa[i]<<"' totalcount='"<<sampleSizes[i]<<"'>\n";
		for(uint j=0;j<alleleCounts.size()-1;j++)
			os<<alleleCounts[j][i]<<",";
		os<<alleleCounts[alleleCounts.size()-1][i]<<",\n";
		os<<"\t\t</sequence>\n";
	}
	os<<"\t</data>\n\n\n";

double alpha=2;
double beta=200;
double lambda=10;


os <<"\n";
os <<"<run id='mcmc' spec='snap.MCMC' chainLength='100000' preBurnin='0' stateBurnin='10000'>\n";
os <<"        <state>\n";
os <<"          <stateNode spec='GammaParameter' id='gamma' initFromTree='false' pattern='gamma' value='10'>\n";
os <<"            <tree idref='tree'/>\n";
os <<"          </stateNode>\n";
os <<"\n";
os <<"          <parameter name='stateNode' id='v' value='"<<v<<"' lower='0.0'/>\n";
os <<"          <parameter name='stateNode' id='u' value='"<<u<<"' lower='0.0'/>\n";
os <<"          <parameter name='stateNode' id='alpha'  value='"<<alpha<<"' lower='0.0'/>\n";
os <<"          <parameter name='stateNode' id='beta'   value='"<<beta<<"' lower='0.0'/>\n";
os <<"          <parameter name='stateNode' id='lambda' value='"<<lambda<<"' lower='0.0'/>\n";
os <<"\n";
os <<"          <tree name='stateNode' spec='ClusterTree' id='tree' nodetype='snap.NodeData' clusterType='upgma'>\n";
os <<"               <input name='taxa' idref='alignment'/>\n";
os <<"          </tree>\n";
os <<"        </state>\n";
os <<"\n";
os <<"        <distribution id='posterior' spec='beast.core.util.CompoundDistribution'>\n";
os <<" 	          <snapprior name='distribution' id='prior'>\n";
os <<"	              <input name='alpha' idref='alpha'/>\n";
os <<"    		      <input name='beta' idref='beta'/>\n";
os <<"	    	      <input name='lambda' idref='lambda'/>\n";
os <<"	              <input name='gamma' idref='gamma'/>\n";
os <<"		          <input name='tree' idref='tree'/>\n";
os <<"            </snapprior>\n";
os <<"            <snaptreelikelihood name='distribution' id='treeLikelihood'>\n";
os <<"                <input name='data' idref='alignment'/>\n";
os <<"                <input name='tree' idref='tree'/>\n";
os <<"                <input name='mutationRateU' idref='u'/>\n";
os <<"                <input name='mutationRateV' idref='v'/>\n";
os <<"	              <input name='gamma' idref='gamma'/>\n";
os <<"            </snaptreelikelihood>\n";
os <<"        </distribution>\n";
os <<"\n";
os <<"        <stateDistribution idref='prior'/>\n";
os <<"\n";
os <<"    	<operator spec='operators.NodeSwapper' weight='0.5'>\n";
os <<"	        <tree name='tree' idref='tree'/>\n";
os <<"    	</operator>\n";
os <<"        <operator spec='operators.NodeBudger' weight='0.5' size='0.5'>\n";
os <<"            <tree name='tree' idref='tree'/>\n";
os <<"        </operator>\n";
os <<"	    <operator spec='beast.evolution.operators.ScaleOperator' scaleFactor='0.25' weight='0.5'>\n";
os <<"	        <tree name='tree' idref='tree'/>\n";
os <<"    	</operator>\n";
os <<"        <operator spec='operators.GammaMover' pGammaMove='0.5' weight='8'>\n";
os <<"	        <parameter name='gamma' idref='gamma'/>\n";
os <<"        </operator>\n";
os <<"        <operator spec='operators.RateMixer' scaleFactors='0.25' weight='1'>\n";
os <<"	        <tree name='tree' idref='tree'/>\n";
os <<"	        <parameter name='gamma' idref='gamma'/>\n";
os <<"        </operator>\n";
os <<"\n";
	//Settings for output of MCMC chain
os <<"        <logger logEvery='100'>\n";
os <<"			  <model idref='posterior'/>\n";
os <<"            <parameter name='log' idref='u'/>\n";
os <<"            <parameter name='log' idref='v'/>\n";
os <<"            <distribution name='log' idref='prior'/>\n";
os <<"            <distribution name='log' idref='treeLikelihood'/>\n";
os <<"            <distribution name='log' idref='posterior'/>\n";
os <<"	          <parameter name='log' idref='gamma'/>\n";
os <<"	          <log spec='snap.ThetaLogger'>\n";
os <<"		          <gamma idref='gamma'/>\n";
os <<"	          </log>\n";
os <<"	          <log spec='beast.evolution.tree.TreeHeightLogger'>\n";
os <<"		         <tree idref='tree'/>\n";
os <<"	          </log>\n";
os <<"        </logger>\n";
os <<"        <logger logEvery='100' fileName='"<<fileroot<<".$(seed).log'>\n";
os <<"	          <model idref='posterior'/>\n";
os <<"            <parameter name='log' idref='u'/>\n";
os <<"            <parameter name='log' idref='v'/>\n";
os <<"            <distribution name='log' idref='prior'/>\n";
os <<"            <distribution name='log' idref='treeLikelihood'/>\n";
os <<"            <distribution name='log' idref='posterior'/>\n";
os <<"			<parameter name='log' idref='gamma'/>\n";
os <<"			<log spec='snap.ThetaLogger'>\n";
os <<"				<gamma idref='gamma'/>\n";
os <<"			</log>\n";
os <<"			<log spec='beast.evolution.tree.TreeHeightLogger'>\n";
os <<"				<tree idref='tree'/>\n";
os <<"			</log>\n";
os <<"        </logger>\n";
os <<"        <logger logEvery='100' fileName='"<<fileroot<<".$(seed).trees'>\n";
os <<"            <tree name='log' idref='tree'/>\n";
os <<"        </logger>\n";
os <<"</run>\n";
os <<"\n";
os <<"\n";
os <<"</snap>\n";
}

void output_nexus(ostream& os, const vector<string>& species, const vector<uint>& sampleSizes, double u, double v, const vector<vector<uint> >&alleleCounts) {
	
	/**
	 Outputs a nexusfile. Taxa (individual names) are generated from the species.
	 **/
	
	int nspecies = species.size();
	int ntax = 0;
	for(uint i=0;i<sampleSizes.size();i++)
		ntax+=sampleSizes[i];
	
	os<<"#NEXUS\n\n";
	
	os<<"[\nOutput from SimSnap\n\nMutation rates:\n\tu = "<<u<<"\tv = "<<v<<"\n]\n\n";
	
	
	os<<"begin taxa;\n";
	os<<"\tdimensions ntax = "<<ntax<<";\n";
	os<<"\ttaxlabels\n";
	for(uint i=0;i<nspecies;i++) 
		for(uint j=0;j<sampleSizes[i];j++) 
			os<<"\t\t'"<<species[i]<<"_"<<(j+1)<<"'\n";
	os<<"\t;\nend;\n\n";
	
	uint nchar = alleleCounts.size();
	os<<"begin characters;\n";
	os<<"\tdimensions nchar = "<<nchar<<";\n";
	os<<"\tformat\n\tdatatype=STANDARD\n\tmissing=?\n\tgap=-\n\tsymbols=\"01\"\n";
	os<<"\tlabels=left\n\ttranspose=no\n\tinterleave=no;\n";
	os<<"\tMATRIX\n";
	
	//Build up the character matrix first. For each marker, randomly assign
	//the alleles conditional on the frequency.
	vector< vector<bool> > matrix;
	matrix.resize(ntax);
	for(int i=0;i<ntax;i++)
		matrix[i].resize(nchar);
	for(int k=0;k<nchar;k++) {
		int id = 0;
		for (int i=0;i<nspecies;i++) {
			vector<bool> ch(sampleSizes[i]);
			fill(ch.begin(),ch.end(),false);
			for(int j=0;j<alleleCounts[k][i];j++)
				ch[j]=true;
			random_shuffle(ch.begin(),ch.end());
			for(int j=0;j<ch.size();j++) {
				matrix[id][k] = ch[j];
				id++;
			}
		}
	}
	
	//Now output the matrix block
		
	int id = 0;
	for(uint i=0;i<nspecies;i++) 
		for(uint j=0;j<sampleSizes[i];j++) {
			os<<"\t\t'"<<species[i]<<"_"<<(j+1)<<"'\t";
			for(uint k=0;k<nchar;k++) {
				if (matrix[id][k])
					os<<"1";
				else
					os<<"0";
			}
			id++;
			os<<"\n";
		}
	os<<"\t;\n";
	os<<"end;\n\n";
	
	os<<"begin traits;\n";
	os<<"\tdimensions nTraits = 1;\n";
	os<<"\tformat labels = yes separator = TAB missing=?;\n";
	os<<"\ttraitlabels species;\n";
	os<<"\tmatrix\n";
	for(uint i=0;i<nspecies;i++) 
		for(uint j=0;j<sampleSizes[i];j++) 
			os<<"\t\t'"<<species[i]<<"_"<<(j+1)<<"'\t'"<<species[i]<<"'\n";
	os<<"\t;\nend;\n\n"<<endl;
}


int main(int argc, char* argv[]) {
	
	ArgumentParser ap(argc,argv);
	
	//Read in the input file.
	//First line is number of taxa
	int nspecies;
	// nr of trees to generate sequences for
	int ntrees;
	ifstream is(ap.inputfile.c_str());
	if (!is) {
		cerr<<"Error reading in file "<<ap.inputfile<<"\n\n";
		printUsage(cerr);
	}
	is>>nspecies;
	if (nspecies<=0)
		printUsage(cerr);

	//Get the prefix from the filename.
	size_t dotpos = ap.inputfile.find_last_of('.');
	string fileroot;
	if (dotpos != string::npos) 
		fileroot = ap.inputfile.substr(0,dotpos);
	else
		fileroot = ap.inputfile;
	
	
	
	//string filemain = 
		
	//Now read in the taxon names and sample sizes
	vector<uint> sampleSizes(nspecies);
	vector<string> species(nspecies);
	for(int i=0;i<nspecies;i++) {
		is>>species[i];
		is>>sampleSizes[i];
	}
	double u,v;
	is>>u;
	is>>v;

	is>>ntrees;
	if (ntrees<=0)
		printUsage(cerr);
	
	cout<<"Simulating data from "<<ntrees<<" species trees\n";
	cout<<"Mutation rates u = "<<u<<" v = "<<v<<"\n";
	cout<<"Number of sites = "<<ap.nsites<<"\n";
	cout<<"Number of trees = "<<ntrees<<"\n";
	cout<<"Output to "<<((ap.outputXML)?"XML":"nexus")<<"\n";
	cout<<((ap.excludeConst)?"Exclude":"Include")<<" constant characters\n\n";
	cout<<((ap.outputTrees)?"Output":"Don't output")<<" trees\n";
	cout<<"Species and sample sizes:\n";
	for(int i=0;i<nspecies;i++)
		cout<<(i+1)<<" "<<species[i]<<"\t"<<sampleSizes[i]<<"\n";
	cout<<endl;
	
	cout<<"Fileroot = "<<fileroot<<endl;
	
	
	for (int iTree = 0; iTree < ntrees; iTree++) {
		//Now read in tree string
		string treeString;
		is>>treeString;
		//TODO: read in semicolons, or at least check for them.

		phylo<basic_newick> tree;
		read_newick(treeString, tree, species, 0.0);
		
		print_newick(cout,tree,true,true);
		cout<<endl;
		
		ostringstream s1;
		s1 << fileroot<<"_tree_"<<(iTree+1);
		if (ap.outputXML)
		  s1<< ".xml";
		else
		  s1 << ".nex";		
		string sFile = s1.str();
		
		cout << "Writing " << sFile << endl;		
		ofstream * os =  new ofstream(sFile.c_str());
		vector<vector<uint> > alleleCounts;
		simulateMultipleSites(tree, u, v, sampleSizes, ap.nsites, ap.excludeConst, alleleCounts, ap.outputTrees);
		
		if (ap.outputXML) {
        	(*os)<<"<!-- Generated with SimSnap -->\n";
            (*os)<<"<!-- -->\n";
            (*os) << "<!-- input tree: ";
            print_newick(*os,tree,true,true);
            (*os) << "-->\n";
			output_xml(*os,species,sampleSizes,u,v,alleleCounts,fileroot);
		} else {
			output_nexus(*os,species,sampleSizes,u,v,alleleCounts);
        }  
        (*os).close();
	}
		
	
}

/*			//Perform posterior predictive check
 vector<vector<uint> > alleleCounts;
 if (iTree==0) {
 (*os)<<"Tree \t SimulatedKappa \t";
 SummaryStatistics::print_header_row(*os,species);
 }
 double trueKappa;	
 simulateMultipleSites(tree, u, v, sampleSizes, ap.nsites, ap.exclude, alleleCounts, trueKappa);
 (*os)<<iTree<<"\t"<<trueKappa<<"\t";
 SummaryStatistics::print_summary_row(*os, sampleSizes, alleleCounts);
 (*os)<<endl;*/


