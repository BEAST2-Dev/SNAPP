package snap.likelihood;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;

@Description("Log distribution normalised by taking out contribution by constant sites")
public class NormalisedDistributionLogger extends BEASTObject implements Loggable {
    public Input<SnAPTreeLikelihood> treeLikelihoodInput = new Input<SnAPTreeLikelihood>("treelikelihood", "SNAPP tree likelihood for the tree", Validate.REQUIRED);
    public Input<Distribution> distributionInput = new Input<Distribution>("distribution", "distirbution to be logged", Validate.REQUIRED);

    SnAPTreeLikelihood treeLikelihood;
    Distribution distribution;
    
    @Override
	public void initAndValidate() {
    	treeLikelihood = treeLikelihoodInput.get();
    	distribution = distributionInput.get();
	}
    
	@Override
	public void init(PrintStream out) {
        out.print(getID() + "\t");
    }

    @Override
    public void log(final long nSample, final PrintStream out) {
        out.print(distribution.getCurrentLogP() - treeLikelihood.getAscSitesLogP() + "\t");
    }

    @Override
    public void close(final PrintStream out) {
        // nothing to do
    }


}
