package snap;


import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.BEASTObject;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;



@Description("Logger to report height of all internal nodes in a tree")
public class TreeNodeLogger extends BEASTObject implements Loggable, Function {
	public Input<Tree> m_tree = new Input<Tree>("tree", "tree to report height for.", Validate.REQUIRED);

	@Override
	public void initAndValidate() {
		// nothing to do
	}

	@Override
	public void init(PrintStream out) {
		final Tree tree = m_tree.get();
		initLog(out, tree.getRoot());
	}

	void initLog(PrintStream out, Node node) {
		if (!node.isLeaf()) {
			out.print("node." +node.getNr()+ ".height\t");
			initLog(out, node.getRight());
			initLog(out, node.getLeft());
		}
	}

	@Override
	public void log(long nSample, PrintStream out) {
		final Tree tree = m_tree.get();
		logTree(out, tree.getRoot());
	}
	
	void logTree(PrintStream out, Node node) {
		if (!node.isLeaf()) {
			out.print(node.getHeight()+ "\t");
			logTree(out, node.getRight());
			logTree(out, node.getLeft());
		}
	}
	

	@Override
	public void close(PrintStream out) {
		// nothing to do
	}

	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public double getArrayValue() {
		return m_tree.get().getRoot().getHeight();
	}

	@Override
	public double getArrayValue(int iDim) {
		return m_tree.get().getRoot().getHeight();
	}
}
