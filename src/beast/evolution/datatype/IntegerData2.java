package beast.evolution.datatype;

import beast.core.Description;

@Description("Datatype for integer sequences")
public class IntegerData2 extends IntegerData {
	
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

    @Override
    public boolean isAmbiguousState(int state) {
        return (state < 0);
    }
}
