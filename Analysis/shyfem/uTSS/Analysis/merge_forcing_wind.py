import os

for i in xrange(1,13):
    sdate="%2.2d" % (i)
    cmd = 'cp '+'2018'+sdate+'01_bsfs/wp.dat 2018year/dnmw'+sdate
    print cmd
    os.system(cmd)



pathname = '/work/mi19918/Projects/uTSS/preproc/runs/2018year'
os.chdir(pathname)

cmd = 'cat dnmw01 dnmw02 dnmw03 dnmw04 dnmw05 dnmw06 dnmw07 dnmw08 dnmw09 dnmw10 dnmw11 dnmw12 > wp.dat'
os.system(cmd)
