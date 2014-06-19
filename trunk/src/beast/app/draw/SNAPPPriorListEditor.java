package beast.app.draw;

import java.util.List;

import javax.swing.Box;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.ListInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
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
        m_plugin = plugin;
		this.itemNr = itemNr;

        m_listBox = Box.createVerticalBox();
        // list of inputs 
        for (Object o : (List<?>) input.get()) {
            if (o instanceof SnAPPrior) {
            	SnAPPrior plugin2 = (SnAPPrior) o;
            	doc.getInpuEditorFactory().addInputs(m_listBox, plugin2, this, null, doc);
            }
        }
		add(m_listBox);
        updateState();
    }

}
