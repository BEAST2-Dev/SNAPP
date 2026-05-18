package snap.spec.likelihood;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.evolution.substitutionmodel.Base;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealScalar;
import beast.base.core.Log;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.tree.Node;

public class SnapSubstitutionModel extends Base {
	public Input<RealScalar<Real>> m_pU = new Input<>("mutationRateU", "Instantaneous rate of mutating from the 0 allele to the 1 alelle");
	public Input<RealScalar<Real>> m_pV = new Input<>("mutationRateV", "Instantaneous rate of mutating from the 1 allele to the 0 alelle");
	public Input<RealVectorParam<PositiveReal>> thetaInput = new Input<>("theta", "population size parameter with one value for each node in the tree");
	public Input<RealVectorParam<PositiveReal>> m_pCoalescenceRate = new Input<>("coalescenceRate", "population size parameter with one value for each node in the tree", Validate.XOR, thetaInput);
	
	public SnapSubstitutionModel() {
		frequenciesInput.setRule(Validate.OPTIONAL);
	}
	
    @Override
    public void initAndValidate() {
    	double u = m_pU.get().get();
    	double v  = m_pV.get().get();
    	if (Math.abs(2*u*v-u-v) > 1e-6) {
    		Log.warning.println("WARNING: Mutation rates are not normalised. "
    				+ "This means that the tree height may not be in units of substitution any more. "
    				+ "To ensure mutation rates are normalised, check that 2 * u * v = u + v, where v and u are the mutation rates.");
    	}
    }
    
	@Override
    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {}

	@Override
	public double[] getFrequencies() {return null;}

	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {return null;}

	@Override
	public boolean canReturnComplexDiagonalization() {return false;}

	@Override
	public boolean canHandleDataType(DataType dataType) {return true;}

	@Override
	public double[] getRateMatrix(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

}
