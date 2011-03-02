package snap.likelihood;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;

public class SnapSubstitutionModel extends SubstitutionModel.Base {
	public Input<RealParameter> m_pU = new Input<RealParameter>("mutationRateU", "mutation rate from red to green?");
	public Input<RealParameter> m_pV = new Input<RealParameter>("mutationRateV", "mutation rate from green to red?");
	public Input<RealParameter> m_pCoalescenceRate = new Input<RealParameter>("coalescenceRate", "population size parameter with one value for each node in the tree");
	
	public SnapSubstitutionModel() {
		frequencies.setRule(Validate.OPTIONAL);
	}
	
    @Override
    public void initAndValidate() {}
    
	@Override
    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {}

	@Override
	public double[] getFrequencies() {return null;}

	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {return null;}

	@Override
	public boolean canReturnComplexDiagonalization() {return false;}

}
