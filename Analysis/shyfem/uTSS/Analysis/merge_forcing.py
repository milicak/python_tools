import os

year = '2018'

for i in xrange(1,13):
    sdate="%2.2d" % (i)
    cmd = 'cp '+year+sdate+'01_mfs/boundn_L1_TSS.dat '+year+'year/dnmeL1'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_mfs/boundn_L2_TSS.dat '+year+'year/dnmeL2'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_bsfs/boundn_L3_TSS.dat '+year+'year/dnmeL3'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_mfs/saltn_L1_TSS.dat '+year+'year/dnms'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_mfs/saltn_L2_TSS.dat '+year+'year/dnmsL2'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_bsfs/saltn_L3_TSS.dat '+year+'year/dnmsL3'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_mfs/tempn_L1_TSS.dat '+year+'year/dnmt'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_mfs/tempn_L2_TSS.dat '+year+'year/dnmtL2'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_bsfs/tempn_L3_TSS.dat '+year+'year/dnmtL3'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_mfs/uv3d_L1_TSS.dat '+year+'year/dnmu'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_mfs/uv3d_L2_TSS.dat '+year+'year/dnmuL2'+sdate
    print cmd
    os.system(cmd)
    cmd = 'cp '+year+sdate+'01_bsfs/uv3d_L3_TSS.dat '+year+'year/dnmuL3'+sdate
    print cmd
    os.system(cmd)
    # cmd = 'cp '+year+sdate+'01/tc.dat '+year+'year/dnmtc'+sdate
    # print cmd
    # os.system(cmd)
    # cmd = 'cp '+year+sdate+'01/wp.dat '+year+'year/dnmw'+sdate
    # print cmd
    # os.system(cmd)



pathname = '/work/mi19918/Projects/uTSS/preproc/runs/'+year+'year'
os.chdir(pathname)

cmd = 'cat dnmeL101 dnmeL102 dnmeL103 dnmeL104 dnmeL105 dnmeL106 dnmeL107 dnmeL108 dnmeL109 dnmeL110 dnmeL111 dnmeL112 > boundn_L1_TSS.dat'
os.system(cmd)
cmd = 'cat dnmeL201 dnmeL202 dnmeL203 dnmeL204 dnmeL205 dnmeL206 dnmeL207 dnmeL208 dnmeL209 dnmeL210 dnmeL211 dnmeL212 > boundn_L2_TSS.dat'
os.system(cmd)
cmd = 'cat dnmeL301 dnmeL302 dnmeL303 dnmeL304 dnmeL305 dnmeL306 dnmeL307 dnmeL308 dnmeL309 dnmeL310 dnmeL311 dnmeL312 > boundn_L3_TSS.dat'
os.system(cmd)
cmd = 'cat dnms01 dnms02 dnms03 dnms04 dnms05 dnms06 dnms07 dnms08 dnms09 dnms10 dnms11 dnms12 > saltn_L1_TSS.dat'
os.system(cmd)
cmd = 'cat dnmsL201 dnmsL202 dnmsL203 dnmsL204 dnmsL205 dnmsL206 dnmsL207 dnmsL208 dnmsL209 dnmsL210 dnmsL211 dnmsL212 > saltn_L2_TSS.dat'
os.system(cmd)
cmd = 'cat dnmsL301 dnmsL302 dnmsL303 dnmsL304 dnmsL305 dnmsL306 dnmsL307 dnmsL308 dnmsL309 dnmsL310 dnmsL311 dnmsL312 > saltn_L3_TSS.dat'
os.system(cmd)
cmd = 'cat dnmt01 dnmt02 dnmt03 dnmt04 dnmt05 dnmt06 dnmt07 dnmt08 dnmt09 dnmt10 dnmt11 dnmt12 > tempn_L1_TSS.dat'
os.system(cmd)
cmd = 'cat dnmtL201 dnmtL202 dnmtL203 dnmtL204 dnmtL205 dnmtL206 dnmtL207 dnmtL208 dnmtL209 dnmtL210 dnmtL211 dnmtL212 > tempn_L2_TSS.dat'
os.system(cmd)
cmd = 'cat dnmtL301 dnmtL302 dnmtL303 dnmtL304 dnmtL305 dnmtL306 dnmtL307 dnmtL308 dnmtL309 dnmtL310 dnmtL311 dnmtL312 > tempn_L3_TSS.dat'
os.system(cmd)
cmd = 'cat dnmu01 dnmu02 dnmu03 dnmu04 dnmu05 dnmu06 dnmu07 dnmu08 dnmu09 dnmu10 dnmu11 dnmu12 > uv3d_L1_TSS.dat'
os.system(cmd)
cmd = 'cat dnmuL201 dnmuL202 dnmuL203 dnmuL204 dnmuL205 dnmuL206 dnmuL207 dnmuL208 dnmuL209 dnmuL210 dnmuL211 dnmuL212 > uv3d_L2_TSS.dat'
os.system(cmd)
cmd = 'cat dnmuL301 dnmuL302 dnmuL303 dnmuL304 dnmuL305 dnmuL306 dnmuL307 dnmuL308 dnmuL309 dnmuL310 dnmuL311 dnmuL312 > uv3d_L3_TSS.dat'
os.system(cmd)
# cmd = 'cat dnmtc01 dnmtc02 dnmtc03 dnmtc04 dnmtc05 dnmtc06 dnmtc07 dnmtc08 dnmtc09 dnmtc10 dnmtc11 dnmtc12 > tc.dat'
# os.system(cmd)
# cmd = 'cat dnmw01 dnmw02 dnmw03 dnmw04 dnmw05 dnmw06 dnmw07 dnmw08 dnmw09 dnmw10 dnmw11 dnmw12 > wp.dat'
# os.system(cmd)

