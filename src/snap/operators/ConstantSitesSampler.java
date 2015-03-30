package snap.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;

@Description("Samples the numberof constant sites from a negative multinomial distribution")
public class ConstantSitesSampler extends Operator {

	public Input<IntegerParameter> ascSiteCountInput = new Input<IntegerParameter>("ascSiteCount", "Counts for number of ascertained sites", Validate.REQUIRED);
	
	IntegerParameter ascSiteCount;
	
	@Override
	public void initAndValidate() throws Exception {
		ascSiteCount = ascSiteCountInput.get();
		if (ascSiteCount.getDimension() != 2) {
			throw new RuntimeException("Expected ascSiteCount to have dimension 2, not " + ascSiteCount.getDimension());
		}
	}

	@Override
	public double proposal() {
		int n0 = ascSiteCount.getValue(0);
		int n1 = ascSiteCount.getValue(0);
		double logHR = 0.0;
		
		// TODO sample n0, n1 
		
		ascSiteCount.setValue(0, n0);
		ascSiteCount.setValue(1, n1);
		return logHR;
	}

}
