package snap.likelihood;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;

public class SnapSubstitutionModel extends CalculationNode implements SubstitutionModel {
	public Input<RealParameter> m_pU = new Input<RealParameter>("mutationRateU", "mutation rate from red to green?");
	public Input<RealParameter> m_pV = new Input<RealParameter>("mutationRateV", "mutation rate from green to red?");
	public Input<RealParameter> m_pCoalescenceRate = new Input<RealParameter>("coalescenceRate", "population size parameter with one value for each node in the tree");
	
    @Override
    public void initAndValidate() {}
    
	@Override
	public void getTransitionProbabilities(double substitutions, double[] matrix) {}

	@Override
	public double[] getFrequencies() {return null;}

	@Override
	public EigenDecomposition getEigenDecomposition() {return null;}

}
