package snap.spec.operators;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.util.Randomizer;

/*
 Change the base frequencies without changing the average rate.

 Actually, we move the pi_0 parameter pi_0 = v/(u+v) and update from there.

 pi_0 ranges from 0 to 1.



 v' = v + X*(u+v)
 u' = u - X*(u+v)

 where X is uniform on [-pi_wind,pi_wind]

 Hastings ratio is 1.0

 J = [X 1+X u+v ; 1-X -X u+v; 0 0 -1]
 which has determinant

 -X^2 - (1+X)(1-X) = -X^2 -1 + X^2 = 1

 */
@Description("Operation that moves both U and V frequency parameters without changing the average rate")
public class MutationMover extends Operator {
	public Input<RealScalarParam<PositiveReal>> m_u = new Input<>("u","frequency of the number of reds", Validate.REQUIRED); 
	public Input<RealScalarParam<PositiveReal>> m_v = new Input<>("v","frequency of the number of greens", Validate.REQUIRED); 
	public Input<Double> m_window = new Input<Double>("window", "window size to move in (defaults to 1)", 1.0); 
    double pi_window;

	@Override
	public void initAndValidate() {
		pi_window = m_window.get();
	}

	@Override
	public double proposal() {
		RealScalarParam<PositiveReal> u = m_u.get();
		RealScalarParam<PositiveReal> v = m_v.get();
            double pi_0 = v.get()/(u.get()+v.get());
            double x = pi_window*(2.0*Randomizer.nextDouble()-1.0);
            //cout<<"old pi_0 = "<<pi_0<<" new pi_0 = "<<pi_0 + x<<endl;

            if ((pi_0+x >= 0.0) && (pi_0 + x <= 1.0)) {
                    pi_0 += x;
                    if (outsideBounds(1.0/(2*pi_0), u) || outsideBounds(1.0/(2*(1-pi_0)), v)) {
                    	return Double.NEGATIVE_INFINITY;
                    }
                    u.set(1.0/(2*pi_0));
                    v.set(1.0/(2*(1-pi_0)));
                    //cout<<"new u,v = "<<toState.u<<"\t"<<toState.v<<endl;
                    return 0;
            }
            return Double.NEGATIVE_INFINITY;
	}

    private boolean outsideBounds(double d, RealScalarParam<PositiveReal> realParameter) {
		if (d < realParameter.getLower() || d > realParameter.getUpper()) {
			return true;
		}
		return false;
	}
}
