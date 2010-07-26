package snap;

import java.io.PrintStream;

import beast.core.Input;
import beast.core.Loggable;
import beast.core.Plugin;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class ThetaLogger extends Plugin implements Loggable {
	public Input<RealParameter> m_pGamma = new Input<RealParameter>("gamma","reports 2 over the value of the parameter.", Validate.REQUIRED);

	
	@Override 
	public void initAndValidate() {
	}
	

	@Override
	public void init(State state, PrintStream out) throws Exception {
		RealParameter param = (RealParameter) m_pGamma.get();
        int nValues = param.getDimension();
        if (nValues == 1) {
            out.print(param.getID() + "\t");
        } else {
            for (int iValue = 0; iValue < nValues; iValue++) {
                out.print("2/"+param.getID() + iValue + "\t");
            }
        }
	}

	@Override
	public void log(int nSample, State state, PrintStream out) {
        RealParameter var = (RealParameter) m_pGamma.get();
        int nValues = var.getDimension();
        for (int iValue = 0; iValue < nValues; iValue++) {
            out.print((2.0/var.getValue(iValue)) + "\t");
        }
	}

	@Override
	public void close(PrintStream out) {
		// nothing to do
	}
}
