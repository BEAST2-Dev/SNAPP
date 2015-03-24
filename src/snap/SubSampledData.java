package snap;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.TaxonSet;
import beast.util.Randomizer;

@Description("Alignment where within species lineages are sampled in order to reduce the dataset")
public class SubSampledData extends Data {
	public Input<Data> dataInput = new Input<Data>("data", "alignment to sample from", Validate.REQUIRED);
	public Input<Integer> minLineagesInput = new Input<Integer>("minLineages", "minimum number of lineages to sample from each species", 1);
	public Input<Double> proportionInput = new Input<Double>("proportion", "proportion of lineages to sample from each species, provided the minimum number is reached", 0.333);
	
	@Override
	public void initAndValidate() throws Exception {
		List<TaxonSet> taxa = m_taxonsets.get();
		List<TaxonSet> subtaxa = new ArrayList<>();
		
		double sampleProportion = proportionInput.get();
		int minLineages = minLineagesInput.get();
		
		for (int i = 0; i < taxa.size(); i++) {
			TaxonSet orgSet = taxa.get(i);
			if (minLineages > orgSet.getTaxonCount()) {
				throw new RuntimeException("minLineages must not exceed the minimum number of lineages per species. " +
						orgSet.getID() + " has only " + orgSet.getTaxonCount() + " lineages, but minLineages is " + minLineages);
			}
			
			TaxonSet subset = new TaxonSet();
			subset.setID(orgSet.getID()+ "subset");
			
			
			boolean [] done = new boolean[orgSet.getTaxonCount()];
			for (int j = 0; j < Math.max(minLineages, sampleProportion * orgSet.getTaxonCount()); j++) {
				int k = Randomizer.nextInt(orgSet.getTaxonCount());
				while (done[k]) {
					k = Randomizer.nextInt(orgSet.getTaxonCount());
				}
				subset.taxonsetInput.setValue(orgSet.getTaxonId(k), subset);
				done[k] = true;
			}
			subtaxa.add(subset);
		}
		
		m_taxonsets.get().clear();
		m_taxonsets.get().addAll(subtaxa);
		
		super.initAndValidate();
	}

}

