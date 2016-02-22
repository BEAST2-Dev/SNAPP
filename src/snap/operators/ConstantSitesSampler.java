package snap.operators;


import snap.distribution.GammaDist;
import snap.likelihood.SnAPTreeLikelihood;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.util.Randomizer;

@Description("Gibbs sampler for the number of constant sites, sampling from a negative multinomial distribution")
public class ConstantSitesSampler extends Operator {

	public Input<IntegerParameter> ascSiteCountInput = new Input<IntegerParameter>("ascSiteCount", "Counts for number of ascertained sites", Validate.REQUIRED);
    public Input<SnAPTreeLikelihood> treeLikelihoodInput = new Input<SnAPTreeLikelihood>("treelikelihood", "SNAPP tree likelihood for the tree", Validate.REQUIRED);
	
	IntegerParameter ascSiteCount;
	SnAPTreeLikelihood treeLikelihood;
	
	@Override
	public void initAndValidate() {
		ascSiteCount = ascSiteCountInput.get();
		if (ascSiteCount.getDimension() != 2) {
			throw new RuntimeException("Expected ascSiteCount to have dimension 2, not " + ascSiteCount.getDimension());
		}
		treeLikelihood = treeLikelihoodInput.get();
	}

	@Override
	public double proposal() {
		int n0 = ascSiteCount.getValue(0);
		int n1 = ascSiteCount.getValue(1);
		int n = n0 + n1;
		double alpha = n;
		//double beta = 1.0/alpha;
		int d = 6; // digits of precission
		double m = GammaDist.inverseF(alpha, 1.0, d, Randomizer.nextDouble());
		
		// alternative implementation -- find out which one is fastest
//		GammaDistribution gamma = new GammaDistributionImpl(alpha, 1/beta);
//		double m = -1;
//		try {
//			m = gamma.inverseCumulativeProbability(Randomizer.nextDouble());
//		} catch (MathException e) {
//			return Double.NEGATIVE_INFINITY;
//		}
		
		double q0  = treeLikelihood.getSitesProbs(-2);
		double q1  = treeLikelihood.getSitesProbs(-1);

		double lambda0 = q0 / (q0 + q1);
		double lambda1 = q1 / (q0 + q1);
		
		
		n0 = (int) Randomizer.nextPoisson(lambda0 * m);
		n1 = (int) Randomizer.nextPoisson(lambda1 * m);
		if (n0 <= ascSiteCount.getLower() || n1 <= ascSiteCount.getLower() || n0 > ascSiteCount.getUpper() || n1 > ascSiteCount.getUpper()) {
			return Double.NEGATIVE_INFINITY;
		}
		ascSiteCount.setValue(0, n0);
		ascSiteCount.setValue(1, n1);
		//System.out.println("NEW: " + n0 + " " + n1 + " " + ascSiteCount);
		return Double.POSITIVE_INFINITY;
	}

}
