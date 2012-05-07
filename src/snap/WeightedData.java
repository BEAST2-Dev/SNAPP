package snap;

import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;

@Description("Alignment with weights attached to patterns")
public class WeightedData extends Data {
	public Input<String> weightInput = new Input<String>("weights", "comma separated list of weights, one for each site in the alignment", Validate.REQUIRED);

	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();
		
		String [] sStr = weightInput.get().split(",");
		
		int siteCount = getSiteCount() + m_nWeight[m_nWeight.length - 2] + m_nWeight[m_nWeight.length - 1];
		if (sStr.length != siteCount) {
			throw new Exception("Number of weights (" + sStr.length + ") does not match number of sites (" + siteCount + ") in alignment");
		}
		Arrays.fill(m_nWeight, 0);
		for (int i = 0; i < sStr.length; i++) {
			int userWeight = Integer.parseInt(sStr[i]);;
			m_nWeight[getPatternIndex(i)] += userWeight;
		}
	} // initAndValidate
	
}
