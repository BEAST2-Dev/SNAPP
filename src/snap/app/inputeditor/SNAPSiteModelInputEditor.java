package snap.app.inputeditor;



import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;


import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.inputeditor.IntegerInputEditor;
import beastfx.app.inputeditor.ParameterInputEditor;
import beastfx.app.inputeditor.SiteModelInputEditor;
import beastfx.app.inputeditor.SmallLabel;
import beastfx.app.util.FXUtils;
import javafx.scene.control.CheckBox;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import beast.base.core.BEASTInterface;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.CompoundDistribution;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.inference.operator.DeltaExchangeOperator;
import beast.base.evolution.sitemodel.SiteModel;
import snap.likelihood.SNAPSiteModel;

public class SNAPSiteModelInputEditor extends SiteModelInputEditor {

	public SNAPSiteModelInputEditor(BeautiDoc doc) {
		super(doc);
	}
	
    @Override
    public Class<?> type() {
        return SNAPSiteModel.class;
    }

    
    private static final long serialVersionUID = 1L;

    IntegerInputEditor categoryCountEditor;
    TextField categoryCountEntry;
    InputEditor gammaShapeEditor;

    // vars for dealing with mean-rate delta exchange operator
    CheckBox fixMeanRatesCheckBox;
    DeltaExchangeOperator operator;
    protected SmallLabel fixMeanRatesValidateLabel;

	    
    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
    		ExpandOption isExpandOption, boolean addButtons) {
    	fixMeanRatesCheckBox = new CheckBox("Fix mean substitution rate");
    	fixMeanRatesCheckBox.setId("FixMeanMutationRate");
    	fixMeanRatesCheckBox.setDisable(doc.autoUpdateFixMeanSubstRate);
    	super.init(input, beastObject, itemNr, isExpandOption, addButtons);

		List<Operator> operators = ((MCMC) doc.mcmc.get()).operatorsInput.get();
    	fixMeanRatesCheckBox.setOnAction(e -> {
				CheckBox averageRatesBox = (CheckBox) e.getSource();
				doFixMeanRates(averageRatesBox.isSelected());
				if (averageRatesBox.isSelected())
					// set up relative weights
					setUpOperator();
			});
    	operator = (DeltaExchangeOperator) doc.pluginmap.get("FixMeanMutationRatesOperator");
    	if (operator == null) {
    		operator = new DeltaExchangeOperator();
    		try {
    			operator.setID("FixMeanMutationRatesOperator");
				operator.initByName("weight", 2.0, "delta", 0.75);
			} catch (Throwable e1) {
				// ignore initAndValidate exception
			}
    		doc.addPlugin(operator);
    	}
		fixMeanRatesCheckBox.setSelected(operators.contains(operator));
		HBox box = FXUtils.newHBox();
		box.getChildren().add(fixMeanRatesCheckBox);
		// box.add(HBox.createHorizontalGlue());
		fixMeanRatesValidateLabel = new SmallLabel("x", "green");
		fixMeanRatesValidateLabel.setVisible(false);
		box.getChildren().add(fixMeanRatesValidateLabel);
		
//    	if (doc.alignments.size() >= 1 && operator != null) {
//        	JComponent component = (JComponent) getComponents()[0];
//    		component.add(box);
//    	}
		pane.getChildren().add(box);
		setUpOperator();
    }
    

	private void doFixMeanRates(boolean averageRates) {
		List<Operator> operators = ((MCMC) doc.mcmc.get()).operatorsInput.get();
		if (averageRates) {
			// connect DeltaExchangeOperator
			if (!operators.contains(operator)) {
				operators.add(operator);
			}
		} else {
			operators.remove(operator);
			fixMeanRatesValidateLabel.setVisible(false);
			repaint();
		}
	}

    public InputEditor createMutationRateEditor() {
    	SiteModel sitemodel = ((SiteModel) m_input.get()); 
        final Input<?> input = sitemodel.muParameterInput;
        ParameterInputEditor mutationRateEditor = new ParameterInputEditor(doc);
        mutationRateEditor.init(input, sitemodel, -1, ExpandOption.FALSE, true);
        mutationRateEditor.getEntry().setDisable(doc.autoUpdateFixMeanSubstRate);
        return mutationRateEditor;
    }
	
	public InputEditor createGammaCategoryCountEditor() {
    	SiteModel sitemodel = ((SiteModel) m_input.get()); 
        final Input<?> input = sitemodel.gammaCategoryCount;
        categoryCountEditor = new IntegerInputEditor(doc) {
			private static final long serialVersionUID = 1L;

			@Override
			public void validateInput() {
        		super.validateInput();
            	SiteModel sitemodel = (SiteModel) m_beastObject; 
                if (sitemodel.gammaCategoryCount.get() < 2 && sitemodel.shapeParameterInput.get().isEstimatedInput.get()) {
                	m_validateLabel.setColor("orange");
                	m_validateLabel.setTooltip(new Tooltip("shape parameter is estimated, but not used"));
                	m_validateLabel.setVisible(true);
                }
        	};
        };
        
        categoryCountEditor.init(input, sitemodel, -1, ExpandOption.FALSE, true);
        categoryCountEntry = categoryCountEditor.getEntry();
        categoryCountEntry.setOnKeyReleased(e->processEntry2());
//        categoryCountEntry.getDocument().addDocumentListener(new DocumentListener() {
//            @Override
//            public void removeUpdate(DocumentEvent e) {
//                processEntry2();
//            }
//
//            @Override
//            public void insertUpdate(DocumentEvent e) {
//                processEntry2();
//            }
//
//            @Override
//            public void changedUpdate(DocumentEvent e) {
//                processEntry2();
//            }
//        });
        
       	categoryCountEditor.validateInput();
        return categoryCountEditor;
    }

    void processEntry2() {
        String categories = categoryCountEntry.getText();
        try {
            int categoryCount = Integer.parseInt(categories);
        	RealParameter shapeParameter = ((SiteModel) m_input.get()).shapeParameterInput.get();
            if (!gammaShapeEditor.getComponent().isVisible() && categoryCount >= 2) {
            	// we are flipping from no gamma to gamma heterogeneity accross sites
            	// so set the estimate flag on the shape parameter
            	shapeParameter.isEstimatedInput.setValue(true, shapeParameter);            	
            } else if (gammaShapeEditor.getComponent().isVisible() && categoryCount < 2) {
            	// we are flipping from with gamma to no gamma heterogeneity accross sites
            	// so unset the estimate flag on the shape parameter
            	shapeParameter.isEstimatedInput.setValue(false, shapeParameter);            	
            }
            Object o = ((ParameterInputEditor)gammaShapeEditor).getComponent();
            if (o instanceof ParameterInputEditor) {
	            ParameterInputEditor e = (ParameterInputEditor) o;
	            e.m_isEstimatedBox.setSelected(shapeParameter.isEstimatedInput.get());
            }
            gammaShapeEditor.getComponent().setVisible(categoryCount >= 2);
            repaint();
        } catch (java.lang.NumberFormatException e) {
            // ignore.
        }
    }

    public InputEditor createShapeEditor() throws NoSuchMethodException, SecurityException, ClassNotFoundException, InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException {
        final Input<?> input = ((SiteModel) m_input.get()).shapeParameterInput;
        gammaShapeEditor = doc.getInputEditorFactory().createInputEditor(input, (BEASTInterface) m_input.get(), doc);
        gammaShapeEditor.getComponent().setVisible(((SiteModel) m_input.get()).gammaCategoryCount.get() >= 2);
        return gammaShapeEditor;
    }

    public static boolean customConnector(BeautiDoc doc) {
 		try {
 	        DeltaExchangeOperator operator = (DeltaExchangeOperator) doc.pluginmap.get("FixMeanMutationRatesOperator");
 	        if (operator == null) {
 	        	return false;
 	        }

 	       	List<RealParameter> parameters = operator.parameterInput.get();
 	    	parameters.clear();
		   	//String weights = "";
		    CompoundDistribution likelihood = (CompoundDistribution) doc.pluginmap.get("likelihood");
		    boolean hasOneEstimatedRate = false;
		    List<String> rateIDs = new ArrayList<>();
		    List<Integer> weights = new ArrayList<>();
			for (Distribution d : likelihood.pDistributions.get()) {
				GenericTreeLikelihood treelikelihood = (GenericTreeLikelihood) d;
	    		Alignment data = treelikelihood.dataInput.get(); 
	    		int weight = data.getSiteCount();
	    		if (data.isAscertained) {
	    			weight -= data.getExcludedPatternCount();
	    		}
	    		if (treelikelihood.siteModelInput.get() instanceof SiteModel) {
		    		SiteModel siteModel = (SiteModel) treelikelihood.siteModelInput.get();
		    		RealParameter mutationRate = siteModel.muParameterInput.get();
		    		//clockRate.m_bIsEstimated.setValue(true, clockRate);
		    		if (mutationRate.isEstimatedInput.get()) {
		    			hasOneEstimatedRate = true;
		    			if (rateIDs.indexOf(mutationRate.getID()) == -1) {
			    			parameters.add(mutationRate);
			    			weights.add(weight);
			    			rateIDs.add(mutationRate.getID());
		    			} else {
		    				int k = rateIDs.indexOf(mutationRate.getID());
			    			weights.set(k,  weights.get(k) + weight);
		    			}
		    		}
	    		}
	    	}
			
			
		    IntegerParameter weightParameter;
			if (weights.size() == 0) {
		    	weightParameter = new IntegerParameter();
			} else {
				String weightString = "";
				for (int k : weights) {
					weightString += k + " ";
				}
		    	weightParameter = new IntegerParameter(weightString);
				weightParameter.setID("weightparameter");
				
			}
			weightParameter.isEstimatedInput.setValue(false, weightParameter);
	    	operator.parameterWeightsInput.setValue(weightParameter, operator);
	    	return hasOneEstimatedRate;
		} catch (Exception e) {
			
		}
		return false;
    }
    
    /** set up relative weights and parameter input **/
    public void setUpOperator() {
    	boolean isAllClocksAreEqual = true;
    	try {
    		boolean hasOneEstimatedRate = customConnector(doc);
		    if (doc.autoUpdateFixMeanSubstRate) {
		    	fixMeanRatesCheckBox.setSelected(hasOneEstimatedRate);
		    	doFixMeanRates(hasOneEstimatedRate);
		    }


     		try {
     	    	double commonClockRate = -1;
    		    CompoundDistribution likelihood = (CompoundDistribution) doc.pluginmap.get("likelihood");
    			for (Distribution d : likelihood.pDistributions.get()) {
    				GenericTreeLikelihood treelikelihood = (GenericTreeLikelihood) d;
    	    		if (treelikelihood.siteModelInput.get() instanceof SiteModel) {
    		    		SiteModel siteModel = (SiteModel) treelikelihood.siteModelInput.get();
    		    		RealParameter mutationRate = siteModel.muParameterInput.get();
    		    		//clockRate.m_bIsEstimated.setValue(true, clockRate);
    		    		if (mutationRate.isEstimatedInput.get()) {
    		    			if (commonClockRate < 0) {
    		    				commonClockRate = mutationRate.valuesInput.get().get(0);
    		    			} else {
    		    				if (Math.abs(commonClockRate - mutationRate.valuesInput.get().get(0)) > 1e-10) {
    		    					isAllClocksAreEqual = false;
    		    				}
    		    			}
    		    		}
    	    		}
    	    	}

    		} catch (Exception e) {
    			
    		}
   		
    		List<RealParameter> parameters = operator.parameterInput.get();
	    	if (!fixMeanRatesCheckBox.isSelected()) {
	    		fixMeanRatesValidateLabel.setVisible(false);
				repaint();
	    		return;
	    	}
	    	if (parameters.size() == 0) {
	    		fixMeanRatesValidateLabel.setVisible(true);
	    		fixMeanRatesValidateLabel.setColor("red");
	    		fixMeanRatesValidateLabel.setTooltip(new Tooltip("The model is invalid: At least one substitution rate should be estimated."));
				repaint();
	    		return;
	    	}
	    	if (!isAllClocksAreEqual) {
	    		fixMeanRatesValidateLabel.setVisible(true);
	    		fixMeanRatesValidateLabel.setColor("orange");
	    		fixMeanRatesValidateLabel.setTooltip(new Tooltip("Not all substitution rates are equal. Are you sure this is what you want?"));
	    	} else if (parameters.size() == 1) {
	    		fixMeanRatesValidateLabel.setVisible(true);
	    		fixMeanRatesValidateLabel.setColor("orange");
	    		fixMeanRatesValidateLabel.setTooltip(new Tooltip("At least 2 clock models should have their rate estimated"));
	    	} else if (parameters.size() < doc.getPartitions("SiteModel").size()) {
	    		fixMeanRatesValidateLabel.setVisible(true);
	    		fixMeanRatesValidateLabel.setColor("orange");
	    		fixMeanRatesValidateLabel.setTooltip(new Tooltip("Not all partitions have their rate estimated"));
	    	} else {
	    		fixMeanRatesValidateLabel.setVisible(false);
	    	}
			repaint();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
