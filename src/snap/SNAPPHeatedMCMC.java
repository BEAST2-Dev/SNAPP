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

	@Override
	public void setChainNr(int i, int resampleEvery) throws Exception {
		if (i > 0) {
			List<Distribution> distrs = ((CompoundDistribution)posterior).pDistributions.get();
			int j = 0;
			while (distrs.get(j).getID().equals("likelihood")) {
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
						"proportion", 2.0/ Math.pow(2,i),
						"taxonset", data.m_taxonsets.get());
			} catch (Exception e) {
				e.printStackTrace();
				throw new RuntimeException("Failed to initialise subsample in SNAPPHeatedMCMC");
			}
			
			
			
			likelihood.dataInput.setValue(subSample, likelihood);
			likelihood.initAndValidate();
		}
		this.resampleEvery = resampleEvery;
	} // setChainNr

} // SNAPPHeatedMCMC
