date=20160101
days=[0,366]
#days=[0,31]
dtt=86400
title='uTSS'
bas='Marmara_basbathy_ser.bas'
simul='uTSS_lobc_chunk_'
param='param_lobc_chunk_'

def setBSUBS(paramfile,day,day0,simul, Nprocs):

	f = open('shympi_chunk_%s.sh'% str(day).zfill(4), 'w')
	f.write('#!/bin/bash\n\n')
	f.write('#BSUB -x\n')
	f.write('#BSUB -q poe_long\n')
	f.write('#BSUB -a poe\n')
	f.write('#BSUB -J uTSS%s\n' % day)
	f.write('#BSUB -n %s\n' % Nprocs)
	f.write('#BSUB -o logout.%J.out\n' )
	f.write('#BSUB -e logerr.%J.err\n' )
	#f.write('#BSUB -sla SC_sanifs_longrun_PS4\n')
	if day>day0:
		f.write('#BSUB -w \"done(uTSS%s)\"\n\n' % (day-1))
	else:
		f.write('\n')
	f.write('module purge\n')
	f.write('module purge\n')
	f.write('module load INTEL/intel_xe_2015.3.187\n')
	f.write('module load IMPI/intel_mpi_5.0.3.048\n')
	f.write('module load NETCDF/netcdf-C_4.3.3.1-F_4.4.2_C++_4.2.1_parallel\n')
	f.write('module load INTEL/intel_xe_2013\n')
	f.write('module load NETCDF/netcdf-4.3\n\n')
	f.write('export LD_LIBRARY_PATH=/users/home/sco116/PETSC/2015/petsc-3.7.4/linux-gnu-intel/lib:${LD_LIBRARY_PATH}\n')
	f.write('export RUNTIME_OPTS=\"-ksp_type bcgs -ksp_rtol 1e-8 -ksp_atol 1e-8\"\n\n')
	#f.write('time mpirun_Impi5 -l /work/tessa_gpfs2/sanifs/mpi_codes/oper_saniv2fixsp/fem3d/shympi %s%s.str $RUNTIME_OPTS\n\n' % (paramfile,str(day).zfill(4)) )
	f.write('time mpirun_Impi5 -l /users/home/mi19918/models/shympi_def_start/fem3d/shympi %s%s.str $RUNTIME_OPTS > logmsS.log\n\n' % (paramfile,str(day).zfill(4)) )
	f.write('echo \"Job completed at: \" `date`\n\n')
	f.write('sleep 100 \n\n') 
	f.write('mv %s%s*.rst restart_files \n\n' %(simul,str(day).zfill(4))) 
	f.write('sleep 100 \n\n') 
	f.write('mkdir out_%s\n\n' %(str(day).zfill(4))) 
	f.write('sleep 100 \n\n') 
	f.write('mv nosnc_conv%s.sh  out_%s\n' %(str(day).zfill(3),str(day).zfill(4))) 
	f.write('mv ousnc_conv%s.sh  out_%s\n' %(str(day).zfill(3),str(day).zfill(4))) 
	f.write('sleep 100 \n\n') 
	f.write('mv %s%s*.nos out_%s\n\n' %(simul,str(day).zfill(4),str(day).zfill(4))) 
	f.write('sleep 100 \n\n') 
	f.write('mv %s%s*.ous out_%s' %(simul,str(day).zfill(4),str(day).zfill(4))) 
	f.close()


for day in range(days[0],days[1]):
    setBSUBS(param,day,days[0],simul, Nprocs=768)
