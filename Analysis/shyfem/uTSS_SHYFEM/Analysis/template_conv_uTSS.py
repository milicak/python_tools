bas='Marmara_basbathy_ser.bas'
simul='uTSS_lobc_chunk_'
param='param_lobc_chunk_'

days=[0,366]
days=[366,732]
days=[0,1096]
# days=[366,732]
# days=[0,732]
# days=[0,1096]
## set the resolution (0 for non-structured output)
resol=0
## set box for zoom (Optional)
#box = [26.1,39.8,30.3,41.4]
#box = [-6.,35.0,-1.,37.0]
#box = [-5.6,35.8,-5.25,36.2]
#box = [28.9,41.0,29.26,41.24]
#[lonmin,latmin ]
#box = [26.0,39.9,26.85,40.55]
#boxname = '.marmara' 
#boxname = '.gibraltar2' 
#boxname = '.bosphorus' 
#boxname = '.dardanelles' 
try:
	box
except NameError:
	box = None
	boxname = ''


## outfolder depends on regular/unstructured case
out_folder = 'OUT'
regroot = ''
## extension for outfile
if resol !=0:
	out_folder = 'nc_reg'
	if resol < 1:		
		regroot = '.reg'+str(resol).replace('.','')
	else:
		regroot = '.reg'+str(resol)


def setCONV_OUS(bas,day,simul):

	f = open('inp_files/ousnc_conv%s.sh'%(str(day).zfill(3)), 'w')
	#f = open('ousnc_conv%s.sh'%day, 'w')
	f.write('#!/bin/bash\n\n')
	#f.write('#BSUB -q serial_30min\n')
	f.write('#BSUB -q s_short\n')
	f.write('#BSUB -P R000\n')
	f.write('#BSUB -J ous2nc%s\n' % day)
	f.write('#BSUB -o logout.%J.out\n' )
	f.write('#BSUB -e logerr.%J.err\n' )
	f.write('#BSUB -w \"done(nos2nc%s)\"\n\n' % (day))
#	f.write('#BSUB -w \"done(amedbs%s)\"\n\n' % day)
	f.write('sim_ous_name=\"%s%s.ous\"\n'% (simul,str(day).zfill(4)))
	f.write('sim_bas_name=\"../%s\"\n\n'% bas)
	f.write('output_name=`echo $sim_ous_name | rev | cut -c 5- | rev`\n\n')
	f.write('echo \"output name=\" $output_name\n\n')
	f.write('module purge\n')
	f.write('module purge\n')
	f.write('module load impi19.5/netcdf-threadsafe/C_4.7.2-F_4.5.2_CXX_4.3.1\n')
	f.write('module load intel19.5/19.5.281\n')
	f.write('module load impi19.5/19.5.281\n')
	f.write('module load impi19.5/petsc/3.7.5\n\n')
	f.write('module load intel19.5/cdo/1.9.8\n\n')
	f.write('module load intel19.5/magics/3.3.1\n\n')
	f.write('module load intel19.5/eccodes/2.12.5\n\n')
	# f.write('cd out_%s\n\n'% (str(day).zfill(4)))
	f.write('/users_home/opa/mi19918/models/shympi_r8_4/fembin/memory -s $sim_ous_name\n')
	f.write('/users_home/opa/mi19918/models/shympi_r8_4/fembin/memory -b $sim_bas_name\n')
	f.write('apn2nc=\"apnous2nc.str\"\n')
	f.write('rm $apn2nc\n\n')
	#f.write('echo -e " " >> "$apn2nc"\n')
	f.write('echo -e " " >> "$apn2nc"\n')
	if resol !=0:
		f.write('echo -e "%.6f " >> "$apn2nc" #window tot\n'% resol)
		if box != None:
			f.write('echo -e "%.2f,%.2f,%.2f,%.2f" >> "$apn2nc" #AFS\n' % (box[0],box[1],box[2],box[3]))
		else:
			f.write('echo -e " " >> "$apn2nc"\n')
	else: 
		f.write('echo -e " " >> "$apn2nc"\n')
		f.write('echo -e " " >> "$apn2nc"\n')
	f.write('echo -e " " >> "$apn2nc"\n')
	f.write('echo -e " " >> "$apn2nc"\n')
	f.write('echo -e " " >> "$apn2nc"\n')
	f.write('time /users_home/opa/mi19918/models/SHYFEM/bin/ous2nc < "$apn2nc" > ous2nc.log 2>&1\n\n')
	f.write('sleep 10\n\n')
	f.write('cdo -f nc4 -z zip_4 copy out.nc out2.nc\n')
	f.write('sleep 10\n\n')
	f.write('rm -f out.nc \n')
	f.write('sleep 10\n\n')
	f.write('cdo timavg out2.nc out3.nc \n')
	f.write('sleep 10\n\n')
	f.write('rm -f out2.nc \n')
	f.write('sleep 10\n\n')
	f.write('mv out3.nc ../%s/%s%s%s%s.ous.nc\n'% (out_folder,simul,str(day).zfill(4),regroot,boxname))
	f.write('echo \"Job completed at: \" `date`\n\n')
	f.close()

def setCONV_NOS(bas,day,simul):

	#f = open('nosnc_conv%s.sh'%day, 'w')
	f = open('inp_files/nosnc_conv%s.sh'%(str(day).zfill(3)), 'w')
	f.write('#!/bin/bash\n\n')
	#f.write('#BSUB -q serial_30min\n')
	f.write('#BSUB -q s_short\n')
	f.write('#BSUB -P R000\n')
	f.write('#BSUB -J nos2nc%s\n' % day)
	f.write('#BSUB -o logout.%J.out\n' )
	f.write('#BSUB -e logerr.%J.err\n' )
#	f.write('#BSUB -w \"done(amedbs%s)\"\n\n' % day)
	f.write('sim_nos_name=\"%s%s.nos\"\n'% (simul,str(day).zfill(4)))
	f.write('sim_bas_name=\"../%s\"\n\n'% bas)
	f.write('output_name=`echo $sim_nos_name | rev | cut -c 5- | rev`\n\n')
	f.write('echo \"output name=\" $output_name\n\n')
	f.write('module purge\n')
	f.write('module purge\n')
	f.write('module load impi19.5/netcdf-threadsafe/C_4.7.2-F_4.5.2_CXX_4.3.1\n')
	f.write('module load intel19.5/19.5.281\n')
	f.write('module load impi19.5/19.5.281\n')
	f.write('module load impi19.5/petsc/3.7.5\n\n')
	f.write('module load intel19.5/cdo/1.9.8\n\n')
	f.write('module load intel19.5/magics/3.3.1\n\n')
	f.write('module load intel19.5/eccodes/2.12.5\n\n')
	f.write('/users_home/opa/mi19918/models/shympi_r8_4/fembin/memory -s $sim_nos_name\n')
	f.write('/users_home/opa/mi19918/models/shympi_r8_4/fembin/memory -b $sim_bas_name\n\n')
	f.write('apn2nc=\"apnnos2nc.str\"\n')
	f.write('rm $apn2nc\n\n')
	f.write('echo -e " " >> "$apn2nc"\n')
	if resol !=0:
		f.write('echo -e "%.6f " >> "$apn2nc" #window tot\n'% resol)
		if box != None:
			f.write('echo -e "%.2f,%.2f,%.2f,%.2f" >> "$apn2nc" #AFS\n' % (box[0],box[1],box[2],box[3]))
	f.write('echo -e " " >> "$apn2nc"\n')
	f.write('echo -e " " >> "$apn2nc"\n')
	f.write('echo -e " " >> "$apn2nc"\n')
#	f.write('time /users/home/ib31217/shyfem-develop/fem3d/nos2nc < "$apn2nc" > nos2nc.log 2>&1\n\n')
	#f.write('time /work/tessa_gpfs2/sanifs/shyfemDeploy/umedbs_op/data/model/bin/nos2nc-enlarged-7_5_13 < "$apn2nc" > nos2nc.log 2>&1\n\n')
	f.write('time /users_home/opa/mi19918/models/SHYFEM/bin/nos2nc < "$apn2nc" > nos2nc.log 2>&1\n\n')
	f.write('sleep 10\n\n')
	f.write('cdo -f nc4 -z zip_4 copy out.nc out2.nc\n')
	f.write('sleep 10\n\n')
	f.write('rm -f out.nc \n')
	f.write('sleep 10\n\n')
	f.write('cdo timavg out2.nc out3.nc \n')
	f.write('sleep 10\n\n')
	f.write('rm -f out2.nc \n')
	f.write('sleep 10\n\n')
	f.write('mv out3.nc ../%s/%s%s%s%s.nos.nc\n'% (out_folder,simul,str(day).zfill(4),regroot,boxname))
	f.write('echo \"Job completed at: \" `date`\n\n')
	f.close()

def setCONV_WND(bas,day,simul):

	f = open('inp_files/wndnc_conv%s.sh'%(str(day).zfill(3)), 'w')
	f.write('#!/bin/bash\n\n')
	f.write('#BSUB -q serial_6h\n')
	f.write('#BSUB -J wnd2nc%s\n' % day)
	f.write('#BSUB -o logout.%J.out\n' )
	f.write('#BSUB -e logerr.%J.err\n' )
	f.write('#BSUB -w \"done(wnd2nc%s)\"\n\n' % (day))
#	f.write('#BSUB -w \"done(amedbs%s)\"\n\n' % day)
	f.write('sim_wnd_name=\"%s%s.wnd\"\n'% (simul,str(day).zfill(4)))
	f.write('sim_bas_name=\"../%s\"\n\n'% bas)
	f.write('output_name=`echo $sim_wnd_name | rev | cut -c 5- | rev`\n\n')
	f.write('echo \"output name=\" $output_name\n\n')
	f.write('module purge\n')
	f.write('module purge\n')
	f.write('module load intel19.5/19.5.281\n')
	f.write('module load impi19.5/netcdf-threadsafe/C_4.7.2-F_4.5.2_CXX_4.3.1\n')
	f.write('module load intel19.5/szip/2.1.1\n')
	f.write('module load impi19.5/hdf5-threadsafe/1.10.5\n')
	f.write('module load impi19.5/petsc/3.7.5\n\n')
	f.write('module load intel19.5/cdo/1.9.8\n\n')
	f.write('module load intel19.5/nco/4.8.1\n\n')
		# f.write('cd out_%s\n\n'% (str(day).zfill(4)))
	f.write('/users_home/opa//mi19918/models/shympi_r8_4/fembin/memory -s $sim_wnd_name\n')
	f.write('/users_home/opa/mi19918/models/shympi_r8_4/fembin/memory -b $sim_bas_name\n\n')
	f.write('apn2nc=\"apnnos2nc.str\"\n')
	f.write('rm $apn2nc\n\n')
	f.write('echo -e " " >> "$apn2nc"\n')
	if resol !=0:
		f.write('echo -e "%.6f " >> "$apn2nc" #window tot\n'% resol)
		if box != None:
			f.write('echo -e "%.2f,%.2f,%.2f,%.2f" >> "$apn2nc" #AFS\n' % (box[0],box[1],box[2],box[3]))
	f.write('echo -e " " >> "$apn2nc"\n')
	f.write('echo -e " " >> "$apn2nc"\n')
	f.write('echo -e " " >> "$apn2nc"\n')
	f.write('time /users_home/opa//mi19918/models/shympiv1_1/2nc/nos2nc < "$apn2nc" > nos2nc.log 2>&1\n\n')
	f.write('sleep 10\n\n')
	# f.write('cdo sellevel,1 out.nc out2.nc\n')
	f.write('ncwa -d level,0 -a level out.nc out2.nc\n')
	f.write('sleep 10\n\n')
	f.write('mv out2.nc out.nc \n')
	f.write('sleep 10\n\n')
	f.write('cdo -f nc4 -z zip_4 copy out.nc out2.nc\n')
	f.write('sleep 10\n\n')
	f.write('rm -f out.nc \n')
	f.write('sleep 10\n\n')
	f.write('cdo timavg out2.nc out3.nc \n')
	f.write('sleep 10\n\n')
	f.write('rm -f out2.nc \n')
	f.write('sleep 10\n\n')
	f.write('mv out3.nc ../%s/%s%s%s%s.wnd.nc\n'% (out_folder,simul,str(day).zfill(4),regroot,boxname))
	f.write('echo \"Job completed at: \" `date`\n\n')
	f.close()

for day in range(days[0],days[1]):
	# setCONV_WND(bas,day,simul)	
	# setCONV_NOS(bas,day,simul)	
	setCONV_OUS(bas,day,simul)	
