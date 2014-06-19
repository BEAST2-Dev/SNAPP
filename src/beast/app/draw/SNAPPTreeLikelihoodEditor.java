package beast.app.draw;

import java.util.List;

import javax.swing.Box;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.app.draw.ListInputEditor;
import beast.app.draw.ParameterInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.evolution.sitemodel.SiteModel;
import snap.likelihood.SnAPTreeLikelihood;
import snap.likelihood.SnapSubstitutionModel;

public class SNAPPTreeLikelihoodEditor extends ListInputEditor {
    public SNAPPTreeLikelihoodEditor(BeautiDoc doc) {
		super(doc);
	}

	private static final long serialVersionUID = 1L;

    public Class<?> baseType() {
        return SnAPTreeLikelihood.class;
    }
    
    SnapSubstitutionModel substModel;
    
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
            if (o instanceof SnAPTreeLikelihood) {
            	SnAPTreeLikelihood plugin2 = (SnAPTreeLikelihood) o;
            	substModel = (SnapSubstitutionModel) ((SiteModel.Base) plugin2.siteModelInput.get()).substModelInput.get();
            	doc.getInpuEditorFactory().addInputs(m_listBox, substModel, this, null, doc);
            	doc.getInpuEditorFactory().addInputs(m_listBox, plugin2, this, null, doc);
            }
        }
		add(m_listBox);
        updateState();
    }
    

    public InputEditor createMutationRateVEditor() throws Exception {
    	ParameterInputEditor editor = (ParameterInputEditor) doc.getInpuEditorFactory().createInputEditor(substModel.m_pV, substModel, doc);
    	editor.m_isEstimatedBox.setVisible(false);
    	return editor;
    }

}
