package beast.evolution.datatype;

import beast.core.Description;
import beast.evolution.datatype.DataType.Base;

@Description("Datatype for integer sequences")
public class IntegerData2 extends Base {
	
	public IntegerData2() {
		m_nStateCount = -1;
		m_mapCodeToStateSet = null;
		m_nCodeLength = -1;
		m_sCodeMap = null;
	}
	
	@Override
	public String getDescription() {
		return "integerdata";
	}
}
