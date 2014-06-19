package beast.evolution.datatype;

import beast.core.Description;
import beast.evolution.datatype.DataType.Base;

@Description("Datatype for integer sequences")
public class IntegerData2 extends Base {
	
	public IntegerData2() {
		stateCount = -1;
		mapCodeToStateSet = null;
		codeLength = -1;
		codeMap = null;
	}
	
	@Override
	public String getTypeDescription() {
		return "integerdata";
	}
}
