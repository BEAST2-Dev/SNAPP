package snap;

import java.io.PrintStream;

import snap.util.AncestralTree;
import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;


@Description("Logs heights of all nodes in an AncestralTree")
public class AncestralTreeHeightLogger extends BEASTObject implements Loggable {
	public Input<AncestralTree> treeInput = new Input<>("tree", "tree that need to be logged", Validate.REQUIRED);
	
	Tree tree;
	@Override
	public void initAndValidate() {
		tree = treeInput.get();

	}

	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < tree.getNodeCount(); i++) {
			out.append("node" + i + "\t");
		}

	}

	@Override
	public void log(long nSample, PrintStream out) {
		for (Node node : tree.getNodesAsArray()) {
			out.append(node.getHeight() + "\t");
		}

	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub

	}

}
