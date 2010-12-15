package snap;

import beast.core.Description;
import beast.core.Input;
import beast.core.Plugin;
import beast.core.Valuable;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

@Description("Converts coalescence rates to thetas, using theta=2/coalescent rate")
public class RateToTheta extends Plugin implements Valuable {
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
