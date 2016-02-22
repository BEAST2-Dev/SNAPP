package snap.likelihood;

import java.io.PrintStream;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;

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
    public void log(final int nSample, final PrintStream out) {
        out.print(distribution.getCurrentLogP() - treeLikelihood.getAscSitesLogP() + "\t");
    }

    @Override
    public void close(final PrintStream out) {
        // nothing to do
    }


}
