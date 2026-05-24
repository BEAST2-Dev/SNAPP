module SNAPP {
    requires transitive beast.pkgmgmt;
    requires transitive beast.base;
    requires transitive beast.fx;
    requires transitive javafx.controls;
	requires java.xml;
	requires java.desktop;
//	requires mtj;
	requires arpack.combined.all;

	requires org.apache.commons.statistics.distribution;

    exports snap;
    exports snap.spec;
    exports snap.datatype;
    exports snap.distribution;
    exports snap.matrix;
    exports snap.tree;
    exports snap.util;
    exports snap.likelihood;
    exports snap.spec.likelihood;
    exports snap.operators;
    exports snap.spec.operators;
    exports snap.app.inputeditor;

	
    provides beast.base.core.BEASTInterface with
			snap.AncestralTreeHeightLogger,
			snap.CoalescentUnitTreeLogger,
			snap.Data,
			snap.datatype.IntegerData2,
			snap.GammaParameter,
			snap.likelihood.BetaApproximationLikelihood,
			snap.likelihood.NormalisedDistributionLogger,
			snap.likelihood.RatePrior,
			snap.likelihood.SnAPPrior,
			snap.spec.likelihood.SnAPPrior,
			snap.likelihood.SNAPSiteModel,
			snap.likelihood.SnapSubstitutionModel,
			snap.spec.likelihood.SnapSubstitutionModel,
			snap.likelihood.SnAPTreeLikelihood,
			snap.spec.likelihood.SnAPTreeLikelihood,
			snap.likelihood.ThresholdTreeLikelihood,
			snap.MCMC,
			snap.ML,
			snap.NodeData,
			snap.operators.ApproximateDistanceBasedLikelihood,
			snap.operators.ApproximateSampledLikelihood,
			snap.operators.BudgerAndScaler,
			snap.operators.CompoundConstantSitesSampler,
			snap.operators.ConstantSitesSampler,
			snap.operators.DelayedAcceptanceOperator,
			snap.operators.DistDAOperator,
			snap.operators.GammaMoveAll,
			snap.operators.GammaMover,
			snap.spec.operators.GammaMover,
			snap.operators.MergeSplitSpeciesTree,
			snap.operators.MutationMover,
			snap.spec.operators.MutationMover,
			snap.operators.NodeBudger,
			snap.operators.NodeSwapper,
			snap.operators.RateMixer,
			snap.spec.operators.RateMixer,
			snap.operators.RootGammaMover,
			snap.operators.ScaleOperator,
			snap.spec.operators.ScaleOperator,
			snap.operators.SubtreeMoveTheta,
			snap.RateToTheta,
			snap.spec.RateToTheta,
			snap.SNPSequence,
			snap.SubSampledData,
			snap.ThetaLogger,
			snap.spec.ThetaLogger,
			snap.tree.YulePriorOneOnXBirthRatePrior,
			snap.TreeLengthLogger,
			snap.TreeNodeLogger,
			snap.util.AncestralTree,
			snap.util.SkylineAnalyser,
			snap.util.TreeSetAnalyser,
			snap.util.TreeSetAnalyser2,
			snap.util.TreeSetAnalyser3,
			snap.WeightedData;

     provides beast.base.evolution.datatype.DataType with
			snap.datatype.IntegerData2;

     provides beastfx.app.inputeditor.InputEditor with
			snap.app.inputeditor.DataInputEditor,
			snap.app.inputeditor.SNAPPPriorListEditor,
			snap.app.inputeditor.SNAPPSubstitutionModelInputEditor,
			snap.app.inputeditor.SNAPPTreeLikelihoodEditor,
			snap.app.inputeditor.SNAPSiteModelInputEditor;
}