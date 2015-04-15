package snap;

import java.util.List;

import snap.likelihood.SnAPTreeLikelihood;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.util.CompoundDistribution;
import beast.inference.HeatedMCMC;

@Description("Heated chain for Metropolis coupled MCMC. "
		+ "Chains are heated by subsampling numbers of lineages")
public class SNAPPHeatedMCMC extends HeatedMCMC {
	int chainNr = 0;
	
	@Override
	public void setChainNr(int i, int resampleEvery) throws Exception {
		if (i > 0) {
			posterior = posteriorInput.get();
			List<Distribution> distrs = ((CompoundDistribution)posterior).pDistributions.get();
			int j = 0;
			while (!distrs.get(j).getID().equals("likelihood")) {
				j++;
			}
			if (j >= distrs.size()) {
				throw new RuntimeException("Expected item with id='likelihood' in posterior");
			}
			distrs = ((CompoundDistribution)distrs.get(j)).pDistributions.get();
			SnAPTreeLikelihood likelihood = (SnAPTreeLikelihood) distrs.get(0);
			Data data = (Data) likelihood.dataInput.get();
			
			
			
			SubSampledData subSample = new SubSampledData();
			try {
				subSample.initByName("data", data, 
						"proportion", 4.0 / (i + 4.0),
						"taxonset", data.m_taxonsets.get());
			} catch (Exception e) {
				e.printStackTrace();
				throw new RuntimeException("Failed to initialise subsample in SNAPPHeatedMCMC");
			}
			
			
			
			likelihood.dataInput.setValue(subSample, likelihood);
			likelihood.initAndValidate();
		}
		this.resampleEvery = resampleEvery;
		
		this.chainNr = i;
	} // setChainNr

	@Override
	public void optimiseRunTime(long startTime, long endTime, long endTimeMainChain) {
		double factor = ((double) endTimeMainChain - startTime) / ((double)endTime - startTime);
		this.resampleEvery = 1 + (int)(this.resampleEvery * factor);
		//System.err.println(this.resampleEvery);
	}

} // SNAPPHeatedMCMC
