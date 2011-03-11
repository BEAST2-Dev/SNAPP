#!/bin/bash


#Assumes that executable for simsnap and jar file for snap are in this directory.


SIMSNAP="./simsnap"

HARDTREE4="(((A[theta=0.006]:0.0057,B[theta=0.006]:0.0057)[theta=0.008]:0.0005,C[theta=0.006]:0.0062)[theta=0.005]:0.0078,D[theta=0.006]:0.014)[theta=0.006]:0.0;"
EASYTREE4="(((A[theta=0.006]:0.0057,B[theta=0.006]:0.0057)[theta=0.005]:0.0045,C[theta=0.006]:0.0102)[theta=0.005]:0.0138,D[theta=0.006]:0.024)[theta=0.006]:0.0;"

HARDTREE8="((((A[theta=0.006]:0.005,B[theta=0.006]:0.005)[theta=0.008]:0.0026,C[theta=0.006]:0.0076)[theta=0.009]:0.0004, D[theta=0.006]:0.008)[theta=0.005]:.010,((E[theta=0.006]:0.003,F[theta=0.006]:0.003)[theta=0.001]:0.004,(G[theta=0.006]:0.0068,H[theta=0.006]:0.0068)[theta=0.014]:0.0002)[theta=0.001]:0.011)[theta=0.011]:0.0;"
EASYTREE8="((((A[theta=0.006]:0.003,B[theta=0.006]:0.003)[theta=0.005]:0.0056,C[theta=0.006]:0.0086)[theta=0.004]:0.0054, D[theta=0.004]:0.014)[theta=0.002]:.004,((E[theta=0.006]:0.003,F[theta=0.006]:0.003)[theta=0.001]:0.009,(G[theta=0.006]:0.0048,H[theta=0.006]:0.0048)[theta=0.004]:0.0072)[theta=0.001]:0.006)[theta=0.021]:0.0;"

TAXANAMES=(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z)


##############################################################################################################################

#Makes input for simsnap. Species names assumed to be A B C D .... 
#Arguments are:
#	filename
#	nspecies
#	nsamples per species
#	Number of replicates of tree
#	tree (as string)

function makeSimSnapInput {
	 
	echo $2 > $1
	for (( i=0; i<$2; i++ ))
	do
		echo ${TAXANAMES[$i]} "   " $3 >> $1
	done
	echo 1 1 >> $1     #Assume unit mutation rates
	echo $4 >> $1   
	for (( i=1; i<=$4 ; i++ ))
	do
		echo $5 >> $1
	done

}


#Modifies the prior values in a SNAP xml input file
#Arguments are:
#	filename
#	alpha value
#	beta value
#	lambda value

function modifySnapXML {
	sed -ie 's/ALPHA/'$2'/g' $1
	sed -ie 's/BETA/'$3'/g' $1
	sed -ie 's/LAMBDA/'$4'/g' $1
	sed -ie 's/LENGTH/'$5'/g' $1
	rm -rf $1'e'
}

function runSnapSingleSim {
	#Creates a temporary directory. Copies the data file to it, runs snap, runs post-analysis, then returns 
	# results of post-analysis, tree file, and log file.
	# if   foo.xml is the name of the data file, then results file is foo_seed_results.txt.
        # Creates foo_seed.tgz which unpacks to the tree file and log file.
	# These are created in the current director.
	# Assumes that there is a copy of snap.jar in the current directory, and uses this.
	#First argument:  xml file to analyse
	#Second argument: Random seed
	#Third argument: true tree, for postanalysis
	len=${#1}  #length of filename
	let len-=4
	FILEROOT=${1:0:len}
	echo $FILEROOT
	
	HERE=$(pwd)
	TMP=$(mktemp -d)
	
	
	cp $1 $TMP
	cd $TMP
	
	echo 'Running SNAP on '$1' in tmp directory '$TMP
	
	#Run snap
	java -jar $HERE/snap.jar  -overwrite -seed $2 $1 &>snapjunk
	
	#Rename output files
	mv *.log $FILEROOT'_seed'$2'.log'
	mv *.trees $FILEROOT'_seed'$2'.trees'
	
	java -cp $HERE/snap.jar  snap.util.TreeSetAnalyser3 -tree $3 *.trees 2> $FILEROOT'_seed'$2'_results.txt' 1>&-
	
	#Back to the original directory
	tar -czf $HERE/$FILEROOT'_seed'$2.tgz $FILEROOT'_seed'$2.log $FILEROOT'_seed'$2'.trees' $1
        
	cp $FILEROOT'_seed'$2'_results.txt' $HERE
	cd $HERE

}

function runSNAPandDelete {
    runSnapSingleSim $1 $2 $3
    rm $1
}

##############################################################################################################################




#Simulates data from 'easy' and 'hard' 4 and 8 taxa trees.
#Runs snap (using 4 combinations of prior parameters)
#For each one, extracts the size of the credible set
# and whether this set includes the true tree.
#Output to file TruePosteriorOut.txt
#Columns are:
# Name of tree
# Gamma prior
# Lambda prior
# Replicate number
# Size of credibility set
# Y or N indicating if true is in credibility set.



makeSimSnapInput thetree.txt 4 1 1 $EASYTREE4 
for nsites in 100 200 300 400 500 600
do
	$SIMSNAP -s $nsites thetree.txt 1> original_$nsites'.txt' 2>&-

	rm thetree.txt
	
#2x2 combinations of prior.

	cp thetree_tree_1.xml 'tree_'$nsites'_aa.xml'
	modifySnapXML 'tree_'$nsites'_aa.xml' 1 200 50 200000
	cp thetree_tree_1.xml 'tree_'$nsites'_ab.xml'
	modifySnapXML 'tree_'$nsites'_ab.xml' 1 200 25 200000
	cp thetree_tree_1.xml 'tree_'$nsites'_ba.xml'
	modifySnapXML 'tree_'$nsites'_ba.xml' 1 2000 50 200000
	cp thetree_tree_1.xml 'tree_'$nsites'_bb.xml'
	modifySnapXML 'tree_'$nsites'_bb.xml' 1 2000 25 200000
	
	rm thetree_tree_1.xml

	for filename in 'tree_'$nsites'_aa.xml' 'tree_'$nsites'_ab.xml' 'tree_'$nsites'_ba.xml' 'tree_'$nsites'_bb.xml'
	do
		 	 runSNAPandDelete $filename 200 $EASYTREE4  &
      	done


done





