package snap;

import java.io.PrintStream;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.StateNode;
import beast.core.BEASTObject;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;



@Description("Logs a tree with branch lengths 2*(branch length in SNAPP) / (theta value for branch)")
public class CoalescentUnitTreeLogger extends BEASTObject implements Loggable {
	public Input<Tree> treeInput = new Input<Tree>("tree", "tree to report on.", Validate.REQUIRED);
	public Input<RealParameter> coalescenceRateInput = new Input<RealParameter>("coalescenceRate", "population size parameter with one value for each node in the tree");

	Tree tree;
	RealParameter coalescenceRate;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		coalescenceRate = coalescenceRateInput.get();
	}

	@Override
	public void init(PrintStream out) {
		tree.init(out);
	}

	@Override
	public void log(long nSample, PrintStream out) {
		// make sure we get the current version of the inputs
       	coalescenceRate = (RealParameter) ((StateNode) coalescenceRate).getCurrent();
        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
		tree.getRoot().sort();
		out.print(toNewick(tree.getRoot()));
        out.print(";");
	}

	
	String toNewick(Node node) {
		StringBuffer buf = new StringBuffer();
		if (node.getLeft() != null) {
			buf.append("(");
			buf.append(toNewick(node.getLeft()));
			if (node.getRight() != null) {
				buf.append(',');
				buf.append(toNewick(node.getRight()));
			}
			buf.append(")");
		} else {
			buf.append(node.getNr() + 1);
		}
	    buf.append(":").append(node.getLength() * coalescenceRate.getValue(node.getNr()));
		return buf.toString();
	}
    
    @Override
	public void close(PrintStream out) {
		tree.close(out);
	}
}
