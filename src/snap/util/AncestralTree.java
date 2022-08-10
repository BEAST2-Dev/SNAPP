package snap.util;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;


@Description("Tree for single population with possibly ancestral nodes")
public class AncestralTree extends Tree {
	public Input<Integer> groupSizeInput = new Input<>("groupSize", "number of groups = number of branches in tree", 2);
	public Input<Double> initialHeightInput = new Input<>("initialHeight", "initial height of the tree", 0.1);
	
	
	@Override
	public void initAndValidate() {
        // make a chain of single node child nodes
        final List<String> sTaxa = m_taxonset.get().asStringList();
        if (sTaxa.size() != 1) {
        	throw new RuntimeException("Expected only a single taxon");
        }
        
        Node left = newNode();
        left.setNr(0);
        left.setHeight(0);
        left.setID(sTaxa.get(0));
        for (int i = 1; i <= groupSizeInput.get(); i++) {
            final Node parent = newNode();
            parent.setNr(i);
            parent.setHeight(i * initialHeightInput.get()/groupSizeInput.get());
            parent.setLeft(left);
            left.setParent(parent);
            left = parent;
        }
        root = left;
        leafNodeCount = 1;
        nodeCount = groupSizeInput.get() + 1;
        internalNodeCount = groupSizeInput.get();

        // RRB: unlikely we need this?
        processTraits(m_traitList.get());

        initArrays();
    }

}
