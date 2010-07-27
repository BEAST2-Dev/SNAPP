/*
 *  posteriorCheck.h
 *  
 *
 *  Created by David Bryant on 22/06/10.
 *  Copyright 2010 University of Auckland. All rights reserved.
 *
 */


/**
 Compute and output posterior predictive check statistics from a chain.
 
 **/
class SummaryStatistics {
public:
	static void print_header_row(ostream& os, vector<string>& species) { 
		os<<"KappaEst\t";
		int nspecies = species.size();
		for (int i=0;i<nspecies;i++) 
			os<<"MeanAllele "<<species[i]<<"\t";
		for (int i=0;i<nspecies;i++) 
			os<<"VarAllele "<<species[i]<<"\t";
		os<<"Covariances\t";
		for(int i=0;i<nspecies;i++)
			os<<"S "<<species[i]<<"\t";
	}
	static void print_summary_row(ostream& os, vector<uint>& sampleSizes, vector<vector<uint> >& alleleCounts) { 
		
		//First estimate the number of constant markers that were removed.
        int nspecies = sampleSizes.size();
        int nchar = alleleCounts[0].size();
        int npairs =0;
        double total = 0.0;
        for(int i=0;i<nspecies;i++) {
            int ni = sampleSizes[i];
            if (ni<2)
                continue;
            for(int j=i+1;j<nspecies;j++) {
                int fi,fj,fij;
                fi = fj = fij=0;
                int nj = sampleSizes[j];
                for(int k=0;k<nchar;k++) {
                    int iCount = alleleCounts[i][k];
                    int jCount = alleleCounts[j][k];
                    if (0<iCount && iCount<ni)
                        fi++;
                    if (0<jCount && jCount<nj)
                        fj++;
                    if ((0<iCount && iCount<ni) && (0<jCount && jCount<nj))
                        fij++;
                }
                if (fi>0 && fj>0){
                    double kappa = 1 - (double)fij/(double)(fi*fj);
                    npairs++;
                    total += kappa;
                }
            }
        }
        double kappa;
        if (npairs>0)
			kappa = (total/npairs);
        else
            kappa = 0.9;
		
		//Now allele count means and variances
		vector<double> alleleMeans(nspecies);
		vector<double> alleleVars(nspecies);
		for(int i=0;i<nspecies;i++) {
			double sum,sum2;
			sum = sum2 = 0.0;
			for(int k=0;k<nchar;k++) {
				double x = (double)alleleCounts[i][k];
				sum+=x;
				sum2+=x*x;
			}
			alleleMeans[i] = sum/(double)nchar;
			alleleVars[i] = sum2/(double)nchar - (sum/(double)nchar)*(sum/(double)nchar);
		}
		
		//Now the average covariance of allele counts across all species.
		double sum = 0.0;
		for(int k=0;k<nchar;k++) {
			double pairsum = 0;
			for(int i=0;i<nspecies;i++)
				for(int j=i+1;j<nspecies;j++)
					pairsum+=alleleCounts[i][k] * alleleCounts[j][k];
			sum+=pairsum;
		}
		double averageCovar = sum/((double)nchar * nspecies*(nspecies-1)/2);
		
		//Now the proportion of segregating sites for each species
		vector<double> pSegregating(nspecies);
		for(int i=0;i<nspecies;i++) {
			sum=0.0;
			int ni = sampleSizes[i];
			for(int k=0;k<nchar;k++) {
				int x = alleleCounts[i][k];
				if (x>0 && x<ni)
					sum+=1.0;
			}
			pSegregating[i]=sum/(double)nchar;
		}
		
		
		
		os<<kappa<<"\t";
		for(int i=0;i<nspecies;i++)
			os<<alleleMeans[i]<<"\t";
		for(int i=0;i<nspecies;i++)
			os<<alleleVars[i]<<"\t";
		os<<averageCovar<<"\t";
		for(int i=0;i<nspecies;i++)
			os<<pSegregating[i]<<"\t";
		
		
		
	}
	
	
	
};


