date=20160101
# days=[366,732]
# days=[0,732]
days=[0,1096]
# days=[0,366]
dtt=86400
title='uTSS'
bas='Marmara_basbathy_ser.bas'
simul='uTSS_lobc_chunk_'
param='param_lobc_chunk_'
#itmrst=43200
#idtrst=43200
itmrst=86400
idtrst=86400
itend=86400
itendConst=0          # itendConst adds constant value to itend
idtout=3600
idtcon=3600
boundfile='bound1.txt'
boundfile2='bound2.txt'
boundfile3='bound3.txt'

"""
levels=[5.0,5.866306,6.977921,8.171897,9.461089,10.86027,12.38635,14.05869,15.89927,17.93301,20.18801,
        22.6957,25.49102,28.61248,32.10208,36.00515,40.36995,45.24715,50.68907,56.7488,63.47917,70.93165,
        79.15538,88.19627,98.09648,108.8941,120.6236,133.3163,147.0015,161.7082,177.4665,194.3095,212.2754,
        231.4089,251.7629,273.4002,296.394,320.8295,346.804,374.4275,403.8227,435.1251,468.482,504.0522,
        542.0048,582.5175,625.775,671.9667,721.2841,773.9175,830.0533,889.8698,953.5339,1021.197,1092.993,
        1169.031,1249.397,1334.148,1423.314,1516.893,1614.854,1717.136,1823.65,1934.282,2048.897,2167.338,
        2289.434,2415.0,2543.845,2675.771,2810.579,2948.069,3088.046,3230.32,3374.708,3521.033,3669.131,
        3818.845,3970.029,4122.546,4276.271,4431.087,4586.889,4743.577,4901.063,5059.266,5218.113,5377.536,6071.000]


levels = [1.0, 2.5, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 25.0, 40.0,60.0,
	 80.0, 100.0, 150.0, 300.0, 500.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6071.0]

"""

levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
          14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
          25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 
          36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46,
          47, 48, 49, 50, 55, 60, 65, 70, 75, 80, 85,
          90, 95, 100, 120, 140, 160, 180, 200, 250, 300,
          350, 400, 450, 500, 550, 600, 650, 700, 750, 800,
          850, 900, 950, 1000, 1100, 1200, 1300, 1400, 1500,
          1600, 1700, 1800, 1900, 2000, 2100, 2200]

f = open(boundfile,'r')
bound1 = f.readlines()
bound1 = [int(i) for i in bound1 ]
f.close()
f = open(boundfile2,'r')
bound2 = f.readlines()
bound2 = [int(i) for i in bound2 ]
f.close()
f = open(boundfile3,'r')
bound3 = f.readlines()
bound3 = [int(i) for i in bound3 ]
f.close()

def setParam(title, day, date, bas, param, simul, dtt, itmrst, idtrst,
             itend,itendConst, idtout, idtcon, levels, bound1, bound2, bound3):

    f = open('%s%s.str'%(param,str(day).zfill(4)), 'w')
    f.write('$title\n\t%s\n'%title)
    #f.write('\t%s_%s\n' % (date, day))
    f.write('\t%s%s\n' % (simul,str(day).zfill(4)))
    f.write('\t%s\n' % bas)
    f.write('$end\n\n$para\n')
    f.write('\tdate = %s\n' % date)
    f.write('\ttime = 000000\n\n')
    #f.write('\trtide = 1\n')
    #f.write('\tltidec = 0.00007\n\n')
    f.write('\tconst = 0\n\n')
    f.write('\tcoumax = 0.9 \n')
    f.write('\titsplt = 2 \n')
    f.write('\tidtsyn = 3600 \n')
    f.write('\tidtmin = 0.5 \n\n')
    f.write('\titanf = %s itend = %s idt = %s\n' % (day*dtt, itend + itendConst+(dtt*day), 40) )
    f.write('\tidtout = %s itmout = %s \n'% (idtout,day*dtt+3600))
    f.write('\tidtcon = %s itmcon = %s \n\n'% (idtcon,day*dtt+3600))
    if day>0:
        f.write('\tityrst = 1\n')
        f.write('\titrst = %s\n'%(day*dtt))

    f.write('\tidtrst = %s\n' % idtrst)
    f.write('\titmrst = %s\n\n' % (itmrst+(dtt*day)))
    f.write('\tilin = 0\n\t'
	    'itlin = 0\n\t'
	    'iclin = 0\n\n\t'
	    'eostype = 2\n\n\t'
	    'rlin = 1\n\n\t'
            'ampar = 0.6\n\t'
            'azpar = 0.6\n\t'
            'aapar = 0\n\n\t'
            'ievap = 1\n\t'
            'isolp = 1\n\t'
            'hdecay = 0\n\t'
            'botabs = 1\n\t'
            'ihtype = 3\n\t'
            'iheat  = 8\n\t'
            'iwtype = 1\n\t'
            'itdrag = 4\n\t'
            'dragco = 0.0025\n\n\t'
            'isphe = 1\n\n\t'
            'noslip = 1\n\n\t'
            'iturb = 1\n\n\t'
            'ireib = 6\n\t'
            'ihwadv = 2\n\t'
            'czdef = 0.01\n\n\t'
            'itvd = 2\n\t'
            'itvdv = 1\n\n\t'
            'ibarcl = 1\n\t'
            'idtau = 1\n\t'
            'itemp = 1\n\t'
            'isalt = 1\n\n\t'
            'idhtyp = 3\n\t'
            'ahpar = 2.2\n\n\t'
            'vismol = 1.e-6\n\t'
            'difmol = 1.e-7\n\t'
            'ilytyp = 3\n\t'
            'hlvmin = 0.5\n\t'
            'vreps  = 1.e-3\n\n\t'
            'chpar = 0.2\n\t'
            'shpar = 0.0\n\t'
            'thpar = 0.0\n\t'
            'dhpar = 0.2\n'
            '$end\n\n$levels\n')
    for level in levels:
        f.write('\t%s\n' % level)
    f.write('$end\n\n')
    f.write('$bound1\n')
    f.write('\tkbound = \n')
    for nknb in bound1:
	f.write('\t%s\n'% nknb)
    f.write('\tibtyp = 1\n')		
    f.write('\tintpol = 2\n')		
    f.write('\ttnudge = 300\n')		
    f.write("\tboundn   = 'input/boundn_L1anew_TSS.dat'\n")		
    f.write("\tsaltn  	= 'input/saltn_L1_TSS.dat'\n")		
    f.write("\ttempn  	= 'input/tempn_L1_TSS.dat'\n")		
    f.write("\tvel3dn 	= 'input/uv3d_L1_TSS.dat'\n")
    f.write('$end\n\n')		
    f.write('$bound2\n')
    f.write('\tkbound = \n')
    for nknb in bound2:
	f.write('\t%s\n'% nknb)
    f.write('\tibtyp = 1\n')		
    f.write('\tintpol = 2\n')		
    f.write('\ttnudge = 300\n')		
    f.write("\tboundn   = 'input/boundn_L2anew_TSS.dat'\n")		
    f.write("\tsaltn  	= 'input/saltn_L2_TSS.dat'\n")		
    f.write("\ttempn  	= 'input/tempn_L2_TSS.dat'\n")		
    f.write("\tvel3dn 	= 'input/uv3d_L2_TSS.dat'\n")
    f.write('$end\n\n')		
    f.write('$bound3\n')
    f.write('\tkbound = \n')
    for nknb in bound3:
	f.write('\t%s\n'% nknb)
    f.write('\tibtyp = 1\n')		
    f.write('\tintpol = 2\n')		
    f.write('\ttnudge = 300\n')		
    f.write("\tboundn   = 'input/boundn_L3anew_TSS.dat'\n")		
    f.write("\tsaltn  	= 'input/saltn_L3_TSS.dat'\n")		
    f.write("\ttempn  	= 'input/tempn_L3_TSS.dat'\n")		
    f.write("\tvel3dn 	= 'input/uv3d_L3_TSS.dat'\n")
    f.write('$end\n\n')		
    f.write('$name\n')
    if day==0:
        f.write("\ttempin = 'input/tempin_%s.dat'\n" % date)
        f.write("\tsaltin = 'input/saltin_%s.dat'\n" % date)
        #f.write("\tuvinit = 'input/uvin_%s.dat'\n" % date)
    else:
        f.write("\trestrt = 'restart_files/%s%s.rst'\n" % (simul,str(day-1).zfill(4)))
    f.write("\twind   = 'input/wp.dat'\n\tqflux  = 'input/tc.dat'\n\train   = 'input/rain.dat'\n\tgotmpa = 'input/gotmturb.nml'\n$end")
    #f.write("\twind   = 'input/wp1.dat'\n\tqflux  = 'input/tc1.dat'\n\tgotmpa = 'input/gotmturb.nml'\n\tqsol = 'input/qsol.dat'\n$end")
    f.close()

# loop on days of simulation to 
# create param files. For each day
# 1 param file is created
for day in range(days[0],days[1]):
    setParam(title, day, date, bas, param, simul, dtt, itmrst, idtrst,
             itend,itendConst, idtout, idtcon, levels, bound1, bound2, bound3)
