package snap.app.inputeditor;


import java.lang.reflect.InvocationTargetException;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BooleanInputEditor;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.inputeditor.ParameterInputEditor;
import beastfx.app.util.FXUtils;
import javafx.scene.control.Button;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.VBox;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.sitemodel.SiteModel;
import snap.Data;
import snap.likelihood.SnAPTreeLikelihood;
import snap.likelihood.SnapSubstitutionModel;

public class SNAPPSubstitutionModelInputEditor extends InputEditor.Base {

	public SNAPPSubstitutionModelInputEditor(BeautiDoc doc) {
		super(doc);
	}

    public Class<?> type() {
        return SnapSubstitutionModel.class;
    }
    
    SnAPTreeLikelihood treeLikelihood;
    SnapSubstitutionModel substModel;
    Data data;
    Button muButton;
    
    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpand, boolean bAddButtons) {
		m_bAddButtons = bAddButtons;
        m_input = input;
        m_beastObject = plugin;
        substModel = (SnapSubstitutionModel) ((SiteModel) plugin).getSubstitutionModel();
        
		this.itemNr = itemNr;

		VBox m_listBox = FXUtils.newVBox();
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

    	muButton = new Button("Calc mutation rates");
    	muButton.setTooltip(new Tooltip("Calcaulate mutation rates based on data in the alignment"));
    	muButton.setOnAction(e -> setUpMutationRates());
    	m_listBox.getChildren().add(muButton);
    	
		try {
	    	ParameterInputEditor editor;
			editor = (ParameterInputEditor) doc.getInputEditorFactory().createInputEditor(substModel.m_pU, substModel, doc);
			m_listBox.getChildren().add(editor);
			editor = (ParameterInputEditor) doc.getInputEditorFactory().createInputEditor(substModel.m_pV, substModel, doc);
	    	editor.m_isEstimatedBox.setVisible(false);
	    	m_listBox.getChildren().add(editor);
		} catch (NoSuchMethodException | SecurityException | ClassNotFoundException | InstantiationException
				| IllegalAccessException | IllegalArgumentException | InvocationTargetException e1) {
			e1.printStackTrace();
		}

    	
    	// doc.getInputEditorFactory().addInputs(m_listBox, substModel, this, null, doc);

    	BooleanInputEditor polyEditor = new BooleanInputEditor(doc);
    	polyEditor.init(treeLikelihood.m_usenNonPolymorphic, treeLikelihood, -1, ExpandOption.FALSE, false);
    	m_listBox.getChildren().add(polyEditor);

    	BooleanInputEditor muAtRootEditor = new BooleanInputEditor(doc);
    	muAtRootEditor.init(treeLikelihood.mutationOnlyAtRoot, treeLikelihood, -1, ExpandOption.FALSE, false);
    	m_listBox.getChildren().add(muAtRootEditor);

    	BooleanInputEditor showPatterLikelihoodEditor = new BooleanInputEditor(doc);
    	showPatterLikelihoodEditor.init(treeLikelihood.showPatternLikelihoodsAndQuit, treeLikelihood, -1, ExpandOption.FALSE, false);
    	m_listBox.getChildren().add(showPatterLikelihoodEditor);

    	BooleanInputEditor useLikelihoodCorrectionEditor = new BooleanInputEditor(doc);
    	useLikelihoodCorrectionEditor.init(treeLikelihood.useLogLikelihoodCorrection, treeLikelihood, -1, ExpandOption.FALSE, false);
    	m_listBox.getChildren().add(useLikelihoodCorrectionEditor);

    	getChildren().add(m_listBox);
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
