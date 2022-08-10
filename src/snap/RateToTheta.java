package snap;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.BEASTObject;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;

@Description("Converts coalescence rates to thetas, using theta=2/coalescent rate")
public class RateToTheta extends BEASTObject implements Function {
	public Input<RealParameter> m_coalescenceRate = new Input<RealParameter>("coalescenceRate","to represent coalescence rate paramater", Validate.REQUIRED);

	@Override
	public void initAndValidate() {
	}
	
	@Override
	public int getDimension() {
		return m_coalescenceRate.get().getDimension();
	}

	@Override
	public double getArrayValue() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getArrayValue(int iDim) {
		return 2.0/m_coalescenceRate.get().getArrayValue(iDim);
	}

}
