package snap;

import java.io.PrintStream;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.BEASTObject;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;



@Description("Logger that reports coalescent rates as theta (using theta=2/rate)")
public class ThetaLogger extends BEASTObject implements Loggable {
	public Input<RealParameter> m_coalescenceRate = new Input<RealParameter>("coalescenceRate","reports 2 over the value of the parameter.", Validate.REQUIRED);

	
	@Override 
	public void initAndValidate() {
	}
	

	@Override
	public void init(PrintStream out) throws Exception {
		RealParameter param = (RealParameter) m_coalescenceRate.get();
        int nValues = param.getDimension();
        if (nValues == 1) {
            out.print(param.getID() + "\t");
        } else {
            for (int iValue = 0; iValue < nValues; iValue++) {
                //out.print("2/"+param.getID() + iValue + "\t");
				out.print("theta" + iValue + "\t"); //The 2/ coalescence rate confuses R (and, at times, me)
            }
        }
	}

	@Override
	public void log(int nSample, PrintStream out) {
        RealParameter var = (RealParameter) m_coalescenceRate.get();
        int nValues = var.getDimension();
        for (int iValue = 0; iValue < nValues; iValue++) {
            out.print((2.0/var.getValue(iValue)) + "\t"); //WARNING: this will be a bug when we allow u and v to change. Value should be 2uv/((u+v)*rate) 
        }
	}

	@Override
	public void close(PrintStream out) {
		// nothing to do
	}
}
