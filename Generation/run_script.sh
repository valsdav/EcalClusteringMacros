#!/bin/sh -e

export X509_USER_PROXY=/afs/cern.ch/user/d/dvalsecc/.proxy
voms-proxy-info

cp -r /afs/cern.ch/work/d/dvalsecc/private/CMSSW_10_6_0 .
cd CMSSW_10_6_0/src

echo -e "evaluate"
eval `scramv1 ru -sh`
export HOME='/afs/cern.ch/user/d/dvalsecc'

JOBID=$1; shift; 
OUTPUTFILE=$1;
NEVENTS=$2;
EMIN=$3
EMAX=$4
ZMIN=$5
ZMAX=$6
RMIN=$7
RMAX=$8
SEED1=$9
SEED2=${10}
SEED3=${11}
SEED4=${12}
DR=${13}
DPHI=${14}
DZ=${15}
EFIX=${16}

cd RecoSimStudies/Dumpers/test

echo -e "cmsRun..";
echo -e ">>> STEP1";
cmsRun step1_CloseEcal_Overlap_cfi_GEN_SIM.py jobid=$JOBID  maxEvents=$NEVENTS     seed1=$SEED1 seed2=$SEED2 seed3=$SEED3 seed4=$SEED4     emin=$EMIN emax=$EMAX zmin=$ZMIN zmax=$ZMAX rmin=$RMIN rmax=$RMAX     deltaR=$DR deltaPhi=$DPHI deltaZ=$DZ efix=$EFIX;

echo -e ">>> STEP2";
cmsRun step2_DIGI_L1_DIGI2RAW_HLT.py
 
#xrdcp --nopbar step2.root root://eosuser.cern.ch/${OUTPUTFILE}_step2.root;

echo -e ">>> STEP3";
cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py

xrdcp --nopbar step3.root root://eosuser.cern.ch/${OUTPUTFILE}_step3_${JOBID}.root;
#xrdcp --nopbar root://eosuser.cern.ch/${OUTPUTFILE}_step3_${JOBID}.root step3.root;

echo -e "Running dumper.."

cd ..
cmsRun python/RecoSimDumper_cfg.py inputFile=test/step3.root outputFile=output_${JOBID}.root


echo -e "Copying result to: $OUTPUTFILE";
xrdcp --nopbar  output_${JOBID}.root root://eosuser.cern.ch/${OUTPUTFILE}_jobid${JOBID}.root;

echo -e "DONE";
