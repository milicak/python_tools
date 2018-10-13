import os

pathname = "/work/mi19918/Projects/uTSS/preproc/runs"
os.chdir(pathname)

mnthnames = ['20160101', '20160201', '20160301', '20160401',
             '20160501', '20160601', '20160701', '20160801',
             '20160901', '20161001', '20161101', '20161201']
# mnthnames = ['20170101', '20170201', '20170301', '20170401',
             # '20170501', '20170601', '20170701', '20170801',
             # '20170901', '20171001', '20171101', '20171201']

for foldername in mnthnames:
        command = "mkdir "+foldername
        os.system(command)
        filenames=['saltn_L1_TSS.dat','saltn_L2_TSS.dat',
                   'boundn_L1_TSS.dat','boundn_L2_TSS.dat',
                   'tempn_L1_TSS.dat','tempn_L2_TSS.dat',
                   'uv3d_L1_TSS.dat','uv3d_L2_TSS.dat',
                   'tc.dat','wp.dat']
        for fnames in filenames:
                command = "cp "+foldername+"/"+fnames+" "+foldername+"/."
                os.system(command)
        # for fnames in filenames:
                # command = "cp "+foldername+"_mfs/"+fnames+" "+foldername+"/."
                # os.system(command)



        filenames=['saltn_L3_TSS.dat',
                   'boundn_L3_TSS.dat',
                    'tempn_L3_TSS.dat',
                    'uv3d_L3_TSS.dat']                      
        for fnames in filenames:
                command = "cp "+foldername+"_bsfs/"+fnames+" "+foldername+"/."
                os.system(command)




