import sys
import os
import argparse
import random
from math import *

'''
This scripts runs hadd on single crystal files to 
group them in strips reading a DOF file
'''
parser = argparse.ArgumentParser()

#parser.add_argument("-f", "--files", type=str, help="input file", required=True)
parser.add_argument("--energy", type=float,nargs='+', help="energies", required=True)
parser.add_argument("--eta", type=float,nargs='+', help="etas", required=True)
parser.add_argument("-n", "--nevents", type=int, help="n events", required=True)
parser.add_argument("-o", "--outputdir", type=str, help="Outputdir", required=True)
parser.add_argument("-c", "--cmssw", type=str, help="CMSSW tar", required=True)
parser.add_argument("-q", "--queue", type=str, help="Condor queue", default="longlunch", required=True)
parser.add_argument("-e", "--eos", type=str, default="user", help="EOS instance user/cms", required=False)
args = parser.parse_args()


# Prepare condor jobs
condor = '''executable              = run_script.sh
output                  = output/strips.$(ClusterId).$(ProcId).out
error                   = error/strips.$(ClusterId).$(ProcId).err
log                     = log/strips.$(ClusterId).log
transfer_input_files    = run_script.sh

+JobFlavour             = "{queue}"
queue arguments from arguments.txt
'''

condor = condor.replace("{queue}", args.queue)
user = os.environ["USER"]

script = '''#!/bin/sh -e

export X509_USER_PROXY=/afs/cern.ch/user/{user1}/{user}/.proxy
voms-proxy-info

cp -r  {cmssw_loc} .
tar zxf {cmssw_file}.tar.gz {cmssw_file}
cd {cmssw_file}/src

echo -e "evaluate"
eval `scramv1 ru -sh`

JOBID=$1; shift; 
OUTPUTFILE=$1;
NEVENTS=$2;
EMIN=$3
EMAX=$4
ZMIN=$5
ZMAX=$6
SEED1=$7
SEED2=$8
SEED3=$9
SEED4=${10}

cd RecoSimStudies/Dumpers/test

echo -e "cmsRun.."
cmsRun step1_CloseEcal_cfi_GEN_SIM.py jobid=$JOBID  maxEvents=$NEVENTS \
    emin=$EMIN emax=$EMAX zmin=$ZMIN zmax=$ZMAX

cmsRun step2_DIGI_L1_DIGI2RAW_HLT.py
cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py

echo -e "Running dumper.."

cd ..
cmsRun python/RecoSimDumper_cfg.py inputFile=test/step3.root outputFile=output_${JOBID}.root


echo -e "Copying result to: $OUTPUTFILE";
xrdcp --nopbar  output_${JOBID}.root root://eos{eosinstance}.cern.ch/${OUTPUTFILE};

echo -e "DONE";
'''

script = script.replace("{eosinstance}", args.eos)
script = script.replace("{user1}", user[:1])
script = script.replace("{user}", user)
cmssw_file = args.cmssw.split("/")[-1]
script = script.replace("{cmssw_loc}", args.cmssw)
script = script.replace("{cmssw_file}", cmssw_file)

arguments= []
if not os.path.exists(args.outputdir):
    os.makedirs(args.outputdir)

def getZ(eta):
    R = 130
    st = 1 / cosh(eta)
    return (R * sqrt(1- st**2)) / st

jobid = 0
for en in args.energy:
    for eta in args.eta:
        z = getZ(eta)
        print("Eta: {} | Z: {}".format(eta, z))

        jobid +=1
        outputfile = args.outputdir + "cluster_{:.1f}_{:.1f}.root".format(en,eta)
        arguments.append("{} {} {} {} {} {} {} {} {} {} {}".format(
            jobid,outputfile,args.nevents, en-0.0001,en+0.0001,z-0.0001,z+0.0001,
            random.randint(1,10000),random.randint(1,10000),
            random.randint(1,10000),random.randint(1,10000)))

print("Njobs: ", len(arguments))
    
with open("condor_job.txt", "w") as cnd_out:
    cnd_out.write(condor)

with open("arguments.txt", "w") as args:
    args.write("\n".join(arguments))

with open("run_script.sh", "w") as rs:
    rs.write(script)

#os.system("condor_submit condor_job.txt")




