package snap.spec;


import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.BEASTObject;
import beast.base.core.Input.Validate;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealVectorParam;



@Description("Logger that reports coalescent rates as theta (using theta=2/rate)")
public class ThetaLogger extends BEASTObject implements Loggable , Function {
	public Input<RealVectorParam<PositiveReal>> m_coalescenceRate = new Input<>("coalescenceRate","reports 2 over the value of the parameter.", Validate.REQUIRED);

	public ThetaLogger() {
		
	}

	@Override 
	public void initAndValidate() {
	}
	

	@Override
	public void init(PrintStream out) {
		RealVectorParam<PositiveReal> param = m_coalescenceRate.get();
        int nValues = param.size();
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
		RealVectorParam<PositiveReal> var = m_coalescenceRate.get();
        int nValues = var.size();
        for (int iValue = 0; iValue < nValues; iValue++) {
            out.print((2.0/var.get(iValue)) + "\t"); //WARNING: this will be a bug when we allow u and v to change. Value should be 2uv/((u+v)*rate) 
        }
	}

	@Override
	public void close(PrintStream out) {
		// nothing to do
	}


	@Override
	public int getDimension() {
		return m_coalescenceRate.get().size();
	}


	@Override
	public double getArrayValue(int dim) {
		RealVectorParam<PositiveReal> var = m_coalescenceRate.get();
        //WARNING: this will be a bug when we allow u and v to change. Value should be 2uv/((u+v)*rate) 
        return (2.0/var.get(dim));
	}
}
