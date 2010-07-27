/*
 *  simSnap.cpp
 *  SimSnap
 *
 *  Created by David Bryant on 8/03/10.
 *  Copyright 2010 University of Auckland. All rights reserved.
 *
 */

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
		
void output_xml(ostream& os, const vector<string>& taxa, const vector<uint>& sampleSizes, double u, double v, const vector<vector<uint> >&alleleCounts) {
	os<<"<snap>\n";
	os<<"\t<taxa id='taxa'>\n";
	for(uint i=0;i<taxa.size();i++) {
		os<<"\t\t<taxon id='"<<taxa[i]<<"'>\n";
		os<<"\t\t<attr name='totalcount'>"<<sampleSizes[i]<<"</attr>\n";
		os<<"\t\t</taxon>\n";
	}
	os<<"\t</taxa>\n";
	
	os<<"\n\n\n\t<!-- n = "<<alleleCounts.size()<<" -->\n";
	//Compute largest sample size.
	int statecount = *max_element(sampleSizes.begin(),sampleSizes.end());
	os<<"\t<alignment id='alignment' dataType='integerdata' statecount='"<<statecount + 1<<"'>\n";
	for(uint i=0;i<taxa.size();i++) {
		os<<"\t\t<sequence>\n\t\t<taxon idref='"<<taxa[i]<<"'/>\n";
		for(uint j=0;j<alleleCounts.size()-1;j++)
			os<<alleleCounts[j][i]<<",";
		os<<alleleCounts[alleleCounts.size()-1][i]<<",\n";
		os<<"\t\t</sequence>\n";
	}
	os<<"\t</alignment>\n\n\n";

double alpha=2;
double beta=0.2;
double lambda=10;
	
os<<""<<endl;
os<<"  <parameter id='v' value='"<<v<<"' lower='0.0'/>"<<endl;
os<<"  <parameter id='u' value='"<<u<<"' lower='0.0'/>"<<endl;
os<<"  <parameter id='alpha'  value='"<<alpha<<"' lower='0.0'/>"<<endl;
os<<"  <parameter id='beta'   value='"<<beta<<"' lower='0.0'/>"<<endl;
os<<"  <parameter id='lambda' value='"<<lambda<<"' lower='0.0'/>"<<endl;
os<<"  <parameter id='gamma' dimension='"<<taxa.size()*2 - 1<<"' value='"<<gamma<<"'/>"<<endl;
os<<""<<endl;
os<<""<<endl;
os<<"  <patterns id='patterns' from='1'>"<<endl;
os<<"         <alignment idref='alignment'/>"<<endl;
os<<"  </patterns>"<<endl;
os<<""<<endl;
os<<""<<endl;
os<<"  <upgmaTree id='startingTree1'>"<<endl;
os<<"         <distanceMatrix correction='none'>"<<endl;
os<<"                 <patterns idref='patterns'/>"<<endl;
os<<"         </distanceMatrix>"<<endl;
os<<"  </upgmaTree>"<<endl;
os<<""<<endl;
os<<"  <constantSize id='initialDemo' units='substitutions'>"<<endl;
os<<"         <populationSize>"<<endl;
os<<"                <parameter id='initialDemo.popSize' value='0.0001'/>"<<endl;
os<<"         </populationSize>"<<endl;
os<<"  </constantSize>"<<endl;
os<<"  <coalescentTree id='startingTree'>"<<endl;
os<<"        <taxa idref='taxa'/>"<<endl;
os<<"        <constantSize idref='initialDemo'/>"<<endl;
os<<"  </coalescentTree>"<<endl;
os<<""<<endl;
os<<"  <report>"<<endl;
os<<"    Starting Tree: <tree idref='startingTree'/>"<<endl;
os<<"  </report>"<<endl;
os<<""<<endl;
os<<""<<endl;
os<<"  <sssModel id='sss'>"<<endl;
os<<"         <alignment idref='alignment'/>"<<endl;
os<<"         <mutationRateV>"<<endl;
os<<"                 <parameter idref='v'/>"<<endl;
os<<"         </mutationRateV>"<<endl;
os<<"         <mutationRateU>"<<endl;
os<<"                 <parameter idref='u'/>"<<endl;
os<<"         </mutationRateU>"<<endl;
os<<"         <gamma><!-- dimension = #external nodes + #internal nodes = #taxa + #taxa-1 -->"<<endl;
os<<"                  <parameter idref='gamma'/>"<<endl;
os<<"         </gamma>"<<endl;
os<<"  </sssModel>"<<endl;
os<<""<<endl;
os<<"  <treeModel id='treeModel'>"<<endl;
os<<"         <tree idref='startingTree'/>"<<endl;
os<<"         <rootHeight>"<<endl;
os<<"                 <parameter id='treeModel.rootHeight'/>"<<endl;
os<<"         </rootHeight>"<<endl;
os<<"         <nodeHeights internalNodes='true'>"<<endl;
os<<"                 <parameter id='treeModel.internalNodeHeights'/>"<<endl;
os<<"         </nodeHeights>"<<endl;
os<<"         <nodeHeights internalNodes='true' rootNode='true'>"<<endl;
os<<"                 <parameter id='treeModel.allInternalNodeHeights'/>"<<endl;
os<<"         </nodeHeights>"<<endl;
os<<"  </treeModel>"<<endl;
os<<""<<endl;
os<<""<<endl;
os<<"  <sssTreeLikelihood id='treeLikelihood'>"<<endl;
os<<"        <patterns idref='patterns'/>"<<endl;
os<<"        <treeModel idref='treeModel'/>"<<endl;
os<<"        <sssModel idref='sss'/>"<<endl;
os<<""<<endl;
os<<"	<treeLength><parameter id='treeLength'/></treeLength>"<<endl;
os<<"        <treeHeight><parameter id='treeHeight'/></treeHeight>"<<endl;
os<<"        <expNumRootLineages><parameter id='expNumRootLineages'/></expNumRootLineages>"<<endl;
os<<"        <meanVarTheta><parameter id='meanVarTheta'/></meanVarTheta>"<<endl;
os<<"        <ThetaRoot><parameter id='ThetaRoot'/></ThetaRoot>"<<endl;
os<<"        <ESS><parameter id='ESS'/></ESS>"<<endl;
os<<"  </sssTreeLikelihood>"<<endl;
os<<""<<endl;
os<<""<<endl;
os<<"    <operators id='operators'>"<<endl;
os<<"	<nodeSwapper weight='0.5'>"<<endl;
os<<"		<treeModel idref='treeModel'/>"<<endl;
os<<"	</nodeSwapper>"<<endl;
os<<"    <subtreeSlide weight='0.5' size='0.5'>"<<endl;
os<<"        <treeModel idref='treeModel'/>"<<endl;
os<<"    </subtreeSlide>"<<endl;
os<<"	<scaleOperator scaleFactor='0.25' weight='0.5'>"<<endl;
os<<"		<parameter idref='treeModel.rootHeight'/>"<<endl;
os<<"	</scaleOperator>"<<endl;
os<<"        <gammaMover pGammaMove='0.5' weight='8'>"<<endl;
os<<"		<parameter idref='gamma'/>"<<endl;
os<<"        </gammaMover>"<<endl;
os<<"        <ratemixer scaleFactors='0.25' weight='1'>"<<endl;
os<<"                <mutationRateU><parameter idref='u'/></mutationRateU>"<<endl;
os<<"                <mutationRateV><parameter idref='v'/></mutationRateV>"<<endl;
os<<"        </ratemixer>"<<endl;
os<<"    </operators>"<<endl;
os<<""<<endl;
os<<"    <mcmc id='mcmc' chainLength='1000' preBurnin='0' autoOptimize='true'>"<<endl;
os<<"        <posterior id='posterior'>"<<endl;
os<<"	    <prior id='prior'>"<<endl;
os<<"		<sssPrior lowerGamma='0.0'>"<<endl;
os<<"			<alpha> <parameter idref='alpha'/> </alpha>"<<endl;
os<<"			<beta>  <parameter idref='beta'/> </beta>"<<endl;
os<<"			<lambda><parameter idref='lambda'/> </lambda>"<<endl;
os<<"	                <mutationRateU><parameter idref='u'/></mutationRateU>"<<endl;
os<<"	                <mutationRateV><parameter idref='v'/></mutationRateV>"<<endl;
os<<"	                <gamma><parameter idref='gamma'/></gamma>"<<endl;
os<<"		        <treeModel idref='treeModel'/>"<<endl;
os<<"		</sssPrior>"<<endl;
os<<"	    </prior>"<<endl;
os<<"            <likelihood id='likelihood'>"<<endl;
os<<"                <treeLikelihood idref='treeLikelihood'/>"<<endl;
os<<"            </likelihood>"<<endl;
os<<"        </posterior>"<<endl;
os<<"        <operators idref='operators'/>"<<endl;
os<<"        <log logEvery='10'>"<<endl;
os<<"            <parameter idref='u'/>"<<endl;
os<<"            <parameter idref='v'/>"<<endl;
os<<"            <column label='Prior' sf='6' width='12'>"<<endl;
os<<"         	   <prior idref='prior'/>"<<endl;
os<<"	    </column>"<<endl;
os<<"            <column label='Likelihood' dp='4' width='12'>"<<endl;
os<<"	            <likelihood idref='likelihood'/>"<<endl;
os<<"	    </column>"<<endl;
os<<"            <column label='Posterior' dp='4' width='12'>"<<endl;
os<<"                <posterior idref='posterior'/>"<<endl;
os<<"            </column>"<<endl;
os<<""<<endl;
os<<"            <column label='meanVarTheta' sf='6' width='12'>"<<endl;
os<<"	            <parameter idref='meanVarTheta'/>"<<endl;
os<<"            </column>"<<endl;
os<<"            <column label='ThetaRoot' sf='6' width='12'>"<<endl;
os<<"	            <parameter idref='ThetaRoot'/>"<<endl;
os<<"            </column>"<<endl;
os<<"            <parameter idref='treeLengthG'/>"<<endl;
os<<"            <parameter idref='treeHeightG'/>"<<endl;
os<<"        </log>"<<endl;
os<<"        <log logEvery='100' fileName='test.$(seed).log'>"<<endl;
os<<"            <prior idref='prior'/>"<<endl;
os<<"            <likelihood idref='likelihood'/>"<<endl;
os<<"            <posterior idref='posterior'/>"<<endl;
os<<"            <parameter idref='u'/>"<<endl;
os<<"            <parameter idref='v'/>"<<endl;
os<<"            <parameter idref='treeLength'/>"<<endl;
os<<"            <parameter idref='treeHeight'/>"<<endl;
os<<"            <parameter idref='treeLengthG'/>"<<endl;
os<<"            <parameter idref='treeHeightG'/>"<<endl;
os<<"            <parameter idref='expNumRootLineages'/>"<<endl;
os<<"            <parameter idref='meanVarTheta'/>"<<endl;
os<<"            <parameter idref='ThetaRoot'/>"<<endl;
os<<"        </log>"<<endl;
os<<"        <logTree logEvery='100' nexusFormat='true' fileName='test.$(seed).trees'>"<<endl;
os<<"            <treeModel idref='treeModel'/>"<<endl;
os<<"        </logTree>"<<endl;
os<<"    </mcmc>"<<endl;
os<<""<<endl;
os<<""<<endl;
os<<"</snap>"<<endl;




/*	
    os<<"\t<sssModel id = 'sss'>\n";
	os<<"\t\t<alignment idref='alignment'/>\n";
	os<<"\t\t<mutationRateV>\n";
	os<<"\t\t<parameter id = 'v' value='"<<v<<"'/>\n";
	os<<"\t\t</mutationRateV>\n";
	os<<"\t\t<mutationRateU>\n";
	os<<"\t\t<parameter id = 'u' value='"<<u<<"'/>\n";
	os<<"\t\t</mutationRateU>\n";
	os<<"\t\t<gamma> <parameter id = 'gamma' dimension = \'"<<taxa.size()*2 - 1<<"\' value = '10.0'/> </gamma>\n";
	os<<"\t</sssModel>\n";

	os<<""<<endl;
	os<<"  <patterns id='patterns' from='1'>"<<endl;
	os<<"         <alignment idref='alignment'/>"<<endl;
	os<<"  </patterns>"<<endl;
	os<<""<<endl;
	os<<"  <upgmaTree id='startingTree'>"<<endl;
	os<<"         <distanceMatrix correction='none'>"<<endl;
	os<<"                 <patterns idref='patterns'/>"<<endl;
	os<<"         </distanceMatrix>"<<endl;
	os<<"  </upgmaTree>"<<endl;
	os<<""<<endl;
	os<<"<treeModel id='treeModel'>"<<endl;
	os<<"         <tree idref='startingTree'/>"<<endl;
	os<<"         <rootHeight>"<<endl;
	os<<"                 <parameter id='treeModel.rootHeight'/>"<<endl;
	os<<"         </rootHeight>"<<endl;
	os<<"         <nodeHeights internalNodes='true'>"<<endl;
	os<<"                 <parameter id='treeModel.internalNodeHeights'/>"<<endl;
	os<<"         </nodeHeights>"<<endl;
	os<<"         <nodeHeights internalNodes='true' rootNode='true'>"<<endl;
	os<<"                 <parameter id='treeModel.allInternalNodeHeights'/>"<<endl;
	os<<"         </nodeHeights>"<<endl;
	os<<"  </treeModel>"<<endl;
	os<<""<<endl;
	os<<"  <sssTreeLikelihood id='treeLikelihood'>"<<endl;
	os<<"        <patterns idref='patterns'/>"<<endl;
	os<<"        <treeModel idref='treeModel'/>"<<endl;
	os<<"        <sssModel idref='sss'/>"<<endl;
	os<<""<<endl;
	os<<"	<treeLength><parameter id='treeLength'/></treeLength>"<<endl;
	os<<"        <treeHeight><parameter id='treeHeight'/></treeHeight>"<<endl;
	os<<"        <expNumRootLineages><parameter id='expNumRootLineages'/></expNumRootLineages>"<<endl;
	os<<"        <meanVarTheta><parameter id='meanVarTheta'/></meanVarTheta>"<<endl;
	os<<"        <ThetaRoot><parameter id='ThetaRoot'/></ThetaRoot>"<<endl;
	os<<"        <ESS><parameter id='ESS'/></ESS>"<<endl;
	os<<"  </sssTreeLikelihood>"<<endl;
	os<<""<<endl;
	os<<""<<endl;
	os<<"    <operators id='operators'>"<<endl;
	os<<"	<nodeSwapper weight='0.5'>"<<endl;
	os<<"		<treeModel idref='treeModel'/>"<<endl;
	os<<"	</nodeSwapper>"<<endl;
	os<<"        <subtreeSlide weight='0.5' size='0.5'>"<<endl;
	os<<"            <treeModel idref='treeModel'/>"<<endl;
	os<<"        </subtreeSlide>"<<endl;
	os<<"	<scaleOperator scaleFactor='0.25' weight='0.5'>"<<endl;
	os<<"		<parameter idref='treeModel.rootHeight'/>"<<endl;
	os<<"	</scaleOperator>"<<endl;
	os<<"        <gammaMover pGammaMove='0.5' weight='8'>"<<endl;
	os<<"		<parameter idref='gamma'/>"<<endl;
	os<<"        </gammaMover>"<<endl;
	os<<"        <ratemixer scaleFactors='0.25' weight='1'>"<<endl;
	os<<"                <mutationRateU><parameter idref='u'/></mutationRateU>"<<endl;
	os<<"                <mutationRateV><parameter idref='v'/></mutationRateV>"<<endl;
	os<<"        </ratemixer>"<<endl;
	os<<"    </operators>"<<endl;
	os<<""<<endl;
	os<<"    <mcmc id='mcmc' chainLength='100000' preBurnin='0' autoOptimize='true'>"<<endl;
	os<<"        <posterior id='posterior'>"<<endl;
	os<<"	    <prior id='prior'>"<<endl;
	os<<"		<sssPrior lowerGamma='0.0'>"<<endl;
	os<<"			<alpha> <parameter id='alpha'  value='2' lower='0.0'/> </alpha>"<<endl;
	os<<"			<beta>  <parameter id='beta'   value='16' lower='0.0'/> </beta>"<<endl;
	os<<"			<lambda><parameter id='lambda' value='2' lower='0.0'/> </lambda>"<<endl;
	os<<"	                <mutationRateU><parameter idref='u'/></mutationRateU>"<<endl;
	os<<"	                <mutationRateV><parameter idref='v'/></mutationRateV>"<<endl;
	os<<"	                <gamma><parameter idref='gamma'/></gamma>"<<endl;
	os<<"		        <treeModel idref='treeModel'/>"<<endl;
	os<<"		</sssPrior>"<<endl;
	os<<"	    </prior>"<<endl;
	os<<"            <likelihood id='likelihood'>"<<endl;
	os<<"                <treeLikelihood idref='treeLikelihood'/>"<<endl;
	os<<"            </likelihood>"<<endl;
	os<<"        </posterior>"<<endl;
	os<<"        <operators idref='operators'/>"<<endl;
	os<<"        <log logEvery='1000'>"<<endl;
	os<<"            <parameter idref='u'/>"<<endl;
	os<<"            <parameter idref='v'/>"<<endl;
	os<<"            <column label='Prior' sf='6' width='12'>"<<endl;
	os<<"         	   <prior idref='prior'/>"<<endl;
	os<<"	    </column>"<<endl;
	os<<"            <column label='Likelihood' dp='4' width='12'>"<<endl;
	os<<"	            <likelihood idref='likelihood'/>"<<endl;
	os<<"	    </column>"<<endl;
	os<<"            <column label='Posterior' dp='4' width='12'>"<<endl;
	os<<"                <posterior idref='posterior'/>"<<endl;
	os<<"            </column>"<<endl;
	os<<""<<endl;
	os<<"            <column label='meanVarTheta' sf='6' width='12'>"<<endl;
	os<<"	            <parameter idref='meanVarTheta'/>"<<endl;
	os<<"            </column>"<<endl;
	os<<"            <column label='ThetaRoot' sf='6' width='12'>"<<endl;
	os<<"	            <parameter idref='ThetaRoot'/>"<<endl;
	os<<"            </column>"<<endl;
	os<<"            <parameter idref='treeLengthG'/>"<<endl;
	os<<"            <parameter idref='treeHeightG'/>"<<endl;
	os<<"        </log>"<<endl;
	os<<"        <log logEvery='1000' fileName='test.$(seed).log'>"<<endl;
	os<<"            <prior idref='prior'/>"<<endl;
	os<<"            <likelihood idref='likelihood'/>"<<endl;
	os<<"            <posterior idref='posterior'/>"<<endl;
	os<<"            <parameter idref='u'/>"<<endl;
	os<<"            <parameter idref='v'/>"<<endl;
	os<<"            <parameter idref='treeLength'/>"<<endl;
	os<<"            <parameter idref='treeHeight'/>"<<endl;
	os<<"            <parameter idref='treeLengthG'/>"<<endl;
	os<<"            <parameter idref='treeHeightG'/>"<<endl;
	os<<"            <parameter idref='expNumRootLineages'/>"<<endl;
	os<<"            <parameter idref='meanVarTheta'/>"<<endl;
	os<<"            <parameter idref='ThetaRoot'/>"<<endl;
	os<<"        </log>"<<endl;
	os<<"        <logTree logEvery='1000' nexusFormat='true' fileName='test.$(seed).trees'>"<<endl;
	os<<"            <treeModel idref='treeModel'/>"<<endl;
	os<<"        </logTree>"<<endl;
	os<<"    </mcmc>"<<endl;

	os<<"\n</beast>"<<endl;
*/	
	
	
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
		
		if (ap.outputXML) 
			output_xml(*os,species,sampleSizes,u,v,alleleCounts);
		else
			output_nexus(*os,species,sampleSizes,u,v,alleleCounts);
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


