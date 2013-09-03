package snap;


import java.io.PrintStream;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.BEASTObject;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;



@Description("Logger to report height of all internal nodes in a tree")
public class TreeNodeLogger extends BEASTObject implements Loggable, Function {
	public Input<Tree> m_tree = new Input<Tree>("tree", "tree to report height for.", Validate.REQUIRED);

	@Override
	public void initAndValidate() {
		// nothing to do
	}

	@Override
	public void init(PrintStream out) throws Exception {
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
	public void log(int nSample, PrintStream out) {
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
