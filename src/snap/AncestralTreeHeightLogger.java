package snap;

import java.io.PrintStream;

import snap.util.AncestralTree;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;


@Description("Logs heights of all nodes in an AncestralTree")
public class AncestralTreeHeightLogger extends BEASTObject implements Loggable {
	public Input<AncestralTree> treeInput = new Input<>("tree", "tree that need to be logged", Validate.REQUIRED);
	
	Tree tree;
	@Override
	public void initAndValidate() {
		tree = treeInput.get();

	}

	@Override
	public void init(PrintStream out) throws Exception {
		for (int i = 0; i < tree.getNodeCount(); i++) {
			out.append("node" + i + "\t");
		}

	}

	@Override
	public void log(int nSample, PrintStream out) {
		for (Node node : tree.getNodesAsArray()) {
			out.append(node.getHeight() + "\t");
		}

	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub

	}

}
