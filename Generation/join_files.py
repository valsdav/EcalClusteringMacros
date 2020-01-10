import os
from multiprocessing import Pool
import sys

inputdir = sys.argv[1]
files = os.listdir(inputdir)
print(files)

nfgroup = int(sys.argv[2])
output = sys.argv[3]
ngroups = len(files) // nfgroup

def join(i):
    fs = files[i*nfgroup: (i+1)*nfgroup ]
    os.system("hadd -f {0}/raw_data_{1}.root {2}".format(output, i," ".join([ inputdir+ "/"+ f for f in fs])))

p = Pool()

p.map(join, range(ngroups))
