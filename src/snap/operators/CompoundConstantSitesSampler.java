package snap.operators;


import snap.distribution.GammaDist;
import snap.likelihood.SnAPTreeLikelihood;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.util.Randomizer;

@Description("Compound sampler for the number of constant sites, sampling from a negative multinomial distribution")
public class CompoundConstantSitesSampler extends Operator {

	public Input<IntegerParameter> ascSiteCountInput = new Input<IntegerParameter>("ascSiteCount", "Counts for number of ascertained sites", Validate.REQUIRED);
    public Input<SnAPTreeLikelihood> treeLikelihoodInput = new Input<SnAPTreeLikelihood>("treelikelihood", "SNAPP tree likelihood for the tree", Validate.REQUIRED);
	public Input<Operator> operatorInput = new Input<>("operator", "base operator used in combination with this operator", Validate.REQUIRED);
	
	IntegerParameter ascSiteCount;
	SnAPTreeLikelihood treeLikelihood;
	Operator operator;
	
	@Override
	public void initAndValidate() {
		ascSiteCount = ascSiteCountInput.get();
		if (ascSiteCount.getDimension() != 2) {
			throw new RuntimeException("Expected ascSiteCount to have dimension 2, not " + ascSiteCount.getDimension());
		}
		treeLikelihood = treeLikelihoodInput.get();
		operator = operatorInput.get();
	}

	@Override
	public double proposal() {
		int siteCount = treeLikelihood.dataInput.get().getSiteCount();

		double q0  = treeLikelihood.getSitesProbs(-2);
		double q1  = treeLikelihood.getSitesProbs(-1);
		double logHR = operator.proposal();
		double [] constProbs = treeLikelihood.calcNewConstProbs();
		double q0new = constProbs[0];
		double q1new = constProbs[1];
		
		int n0 = ascSiteCount.getValue(0);
		int n1 = ascSiteCount.getValue(1);
		//int n = n0 + n1;
		double alpha = siteCount;
		//double beta = 1.0/alpha;
		int d = 6; // digits of precision
		double m = GammaDist.inverseF(alpha, 1.0, d, Randomizer.nextDouble());
		

		double lambda0 = q0new / (1 - q0new - q1new);
		double lambda1 = q1new / (1 - q0new - q1new);
		
		int n0new = (int) Randomizer.nextPoisson(lambda0 * m);
		int n1new = (int) Randomizer.nextPoisson(lambda1 * m);
		if (n0new <= ascSiteCount.getLower() || n1new <= ascSiteCount.getLower() || n0new > ascSiteCount.getUpper() || n1new > ascSiteCount.getUpper()) {
			System.err.println("PROBLEM: Generated strange Ascertained site count");
			return Double.NEGATIVE_INFINITY;
		}
		ascSiteCount.setValue(0, n0new);
		ascSiteCount.setValue(1, n1new);
		//System.out.println("NEW: " + n0new + " " + n1new + " " + ascSiteCount);
		logHR += siteCount * Math.log((1.0 - q0 - q1)/(1.0 - q0new - q1new));
		logHR += n0 * Math.log(q0) + n1 * Math.log(q1) - n0new * Math.log(q0new) - n1new * Math.log(q1new);
		
		return logHR;
	}

}
