#!/bin/bash


cd ../
svn update
ant
mv dist/snap.jar BestSimulations/
cd SimSnap
make
mv simsnap ../BestSimulations
cd ../BestSimulations


