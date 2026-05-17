package snap;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.BEASTObject;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;



@Description("Logger that reports coalescent rates as theta (using theta=2/rate)")
public class ThetaLogger extends BEASTObject implements Loggable , Function {
	public Input<RealParameter> m_coalescenceRate = new Input<RealParameter>("coalescenceRate","reports 2 over the value of the parameter.", Validate.REQUIRED);

	
	@Override 
	public void initAndValidate() {
	}
	

	@Override
	public void init(PrintStream out) {
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
	public void log(long nSample, PrintStream out) {
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


	@Override
	public int getDimension() {
		return m_coalescenceRate.get().getDimension();
	}


	@Override
	public double getArrayValue(int dim) {
        RealParameter var = (RealParameter) m_coalescenceRate.get();
        //WARNING: this will be a bug when we allow u and v to change. Value should be 2uv/((u+v)*rate) 
        return (2.0/var.getValue(dim));
	}
}
