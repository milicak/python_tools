import os

for i in xrange(1,13):
    sdate="%2.2d" % (i)
    cmd = 'cp '+'2018'+sdate+'01_bsfs/tc.dat 2018year/dnmtc'+sdate
    print cmd
    os.system(cmd)



pathname = '/work/mi19918/Projects/uTSS/preproc/runs/2018year'
os.chdir(pathname)

cmd = 'cat dnmtc01 dnmtc02 dnmtc03 dnmtc04 dnmtc05 dnmtc06 dnmtc07 dnmtc08 dnmtc09 dnmtc10 dnmtc11 dnmtc12 > tc.dat'
os.system(cmd)
