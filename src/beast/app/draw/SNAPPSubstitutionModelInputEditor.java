package beast.app.draw;


import java.lang.reflect.InvocationTargetException;

import javax.swing.Box;
import javax.swing.JButton;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.app.draw.ParameterInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.sitemodel.SiteModel;
import snap.Data;
import snap.likelihood.SnAPTreeLikelihood;
import snap.likelihood.SnapSubstitutionModel;

public class SNAPPSubstitutionModelInputEditor extends InputEditor.Base {

	public SNAPPSubstitutionModelInputEditor(BeautiDoc doc) {
		super(doc);
	}

	private static final long serialVersionUID = 1L;

    public Class<?> type() {
        return SnapSubstitutionModel.class;
    }
    
    SnAPTreeLikelihood treeLikelihood;
    SnapSubstitutionModel substModel;
    Data data;
    JButton muButton;
    
    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpand, boolean bAddButtons) {
		m_bAddButtons = bAddButtons;
        m_input = input;
        m_beastObject = plugin;
        substModel = (SnapSubstitutionModel) ((SiteModel) plugin).getSubstitutionModel();
        
		this.itemNr = itemNr;

		Box m_listBox = Box.createVerticalBox();
        // list of inputs
		for (Object o : m_beastObject.getOutputs()) {
			if (o instanceof SnAPTreeLikelihood) {
				m_beastObject = (SnAPTreeLikelihood) o;
			}
		}
		treeLikelihood = (SnAPTreeLikelihood) m_beastObject;
		
		
    	if (treeLikelihood.dataInput.get() instanceof Data) {
    		data = (Data) treeLikelihood.dataInput.get();
    	}

    	muButton = new JButton("Calc mutation rates");
    	muButton.setToolTipText("Calcaulate mutation rates based on data in the alignment");
    	muButton.addActionListener(e -> setUpMutationRates());
    	m_listBox.add(muButton);
    	
		try {
	    	ParameterInputEditor editor;
			editor = (ParameterInputEditor) doc.getInputEditorFactory().createInputEditor(substModel.m_pU, substModel, doc);
			m_listBox.add(editor);
			editor = (ParameterInputEditor) doc.getInputEditorFactory().createInputEditor(substModel.m_pV, substModel, doc);
	    	editor.m_isEstimatedBox.setVisible(false);
	    	m_listBox.add(editor);
		} catch (NoSuchMethodException | SecurityException | ClassNotFoundException | InstantiationException
				| IllegalAccessException | IllegalArgumentException | InvocationTargetException e1) {
			e1.printStackTrace();
		}

    	
    	// doc.getInputEditorFactory().addInputs(m_listBox, substModel, this, null, doc);

    	BooleanInputEditor polyEditor = new BooleanInputEditor(doc);
    	polyEditor.init(treeLikelihood.m_usenNonPolymorphic, treeLikelihood, -1, ExpandOption.FALSE, false);
    	m_listBox.add(polyEditor);

    	BooleanInputEditor muAtRootEditor = new BooleanInputEditor(doc);
    	muAtRootEditor.init(treeLikelihood.mutationOnlyAtRoot, treeLikelihood, -1, ExpandOption.FALSE, false);
    	m_listBox.add(muAtRootEditor);

    	BooleanInputEditor showPatterLikelihoodEditor = new BooleanInputEditor(doc);
    	showPatterLikelihoodEditor.init(treeLikelihood.showPatternLikelihoodsAndQuit, treeLikelihood, -1, ExpandOption.FALSE, false);
    	m_listBox.add(showPatterLikelihoodEditor);

    	BooleanInputEditor useLikelihoodCorrectionEditor = new BooleanInputEditor(doc);
    	useLikelihoodCorrectionEditor.init(treeLikelihood.useLogLikelihoodCorrection, treeLikelihood, -1, ExpandOption.FALSE, false);
    	m_listBox.add(useLikelihoodCorrectionEditor);

		add(m_listBox);
        //updateState();
    }
    

    private Object setUpMutationRates() {
    	double proportionZeros = data.getProportionZeros();
    	double muU = 1 / (2.0 * (1.0 - proportionZeros));
    	double muV = 1 / (2.0 * proportionZeros);
    	RealParameter pU = substModel.m_pU.get();
    	pU.valuesInput.setValue(muU + "", pU);
    	RealParameter pV = substModel.m_pV.get();
    	pV.valuesInput.setValue(muV + "", pV);
    	refreshPanel();
		return null;
	}


	public InputEditor createMutationRateVEditor() throws Exception {
    	ParameterInputEditor editor = (ParameterInputEditor) doc.getInputEditorFactory().createInputEditor(substModel.m_pV, substModel, doc);
    	editor.m_isEstimatedBox.setVisible(false);
    	return editor;
    }

}
