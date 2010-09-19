/**
 * @version $Id: NexusBlock.java,v 1.9 2007-09-11 12:30:59 kloepper Exp $
 *
 * @author Daniel Huson
 *
 */

package nexus;

import javax.swing.tree.DefaultMutableTreeNode;
import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;


/**
 * NexusBlock base class
 */
abstract public class NexusBlock {
    private boolean valid = true;
      /**
     * Identification string
     */
      final static String NAME = "NEXUSBLOCK";

    /**
     * set valid state
     *
     * @param valid
     */
    public void setValid(boolean valid) {
        this.valid = valid;
    }

    /**
     * get valid state
     *
     * @return valid?
     */
    public boolean isValid() {
        return valid;
    }

    /**
     * write a block, blocks should override this
     *
     * @param w
     * @param taxa
     * @throws IOException
     */
    public abstract void write(Writer w, Taxa taxa) throws IOException;

    /**
     * Return the name of this block. Equivalent to 'BlockClass.NAME'
     *
     * @return String name.
     */
    public static String getBlockName() {
        return NAME;
    }


    /**
     * Returns a tree structure containing the data in the block. Subclasses should override
     * this if they want more sophisticated trees. Default is to present entire block as
     * one object.
     *
     * @param taxa
     * @return DefaultMutableTreeNode root of subtree
     */
    public DefaultMutableTreeNode getDataTree(Taxa taxa) {
        Writer w = new StringWriter();
        try {
            write(w, taxa);
        } catch (IOException e) {
            return new DefaultMutableTreeNode(NAME + "<output error>");
        }
        return new DefaultMutableTreeNode(w.toString());
    }
}

// EOF
