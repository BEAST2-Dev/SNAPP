package snap.app.inputeditor;

import java.util.List;

import javax.swing.Box;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.ListInputEditor;
import beastfx.app.util.FXUtils;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import snap.likelihood.SnAPPrior;

public class SNAPPPriorListEditor extends ListInputEditor {
    public SNAPPPriorListEditor(BeautiDoc doc) {
		super(doc);
	}

	private static final long serialVersionUID = 1L;

    public Class<?> baseType() {
        return SnAPPrior.class;
    }
    
    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpand, boolean bAddButtons) {
		m_bAddButtons = bAddButtons;
    	m_bExpandOption = bExpand;
        m_input = input;
        m_beastObject = plugin;
		this.itemNr = itemNr;

        m_listBox = FXUtils.newVBox();
        // list of inputs 
        for (Object o : (List<?>) input.get()) {
            if (o instanceof SnAPPrior) {
            	SnAPPrior plugin2 = (SnAPPrior) o;
            	doc.getInputEditorFactory().addInputs(m_listBox, plugin2, this, null, doc);
            }
        }
		getChildren().add(m_listBox);
        updateState();
    }

}
