package snap;

import java.io.PrintStream;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.BEASTObject;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;



@Description("Logger to report sum of node heights of a tree")
public class TreeLengthLogger extends BEASTObject implements Loggable {
	public Input<Tree> m_tree = new Input<Tree>("tree", "tree to report on.", Validate.REQUIRED);

	@Override
	public void initAndValidate() {
		// nothing to do
	}

	@Override
	public void init(PrintStream out) throws Exception {
		final Tree tree = m_tree.get();
		out.print(tree.getID() + ".length\t");
	}

	@Override
	public void log(int nSample, PrintStream out) {
        Tree tree = (Tree) m_tree.get();
        double heightsum = tree.getRoot().getHeight();
        heightsum += heightSum(tree.getRoot());
		out.print(heightsum + "\t");
	}

    double heightSum(Node node) {
        if (node.isLeaf()) {
            return 0;
        } else {
            double h = node.getHeight();
            h += heightSum(node.getLeft());
            if (node.getRight() != null) {
                h += heightSum(node.getRight());
            }
            return h;
        }
    } // heightSum
    
    @Override
	public void close(PrintStream out) {
		// nothing to do
	}
}
