date=20160101
# days=[0,732]
#days=[366,732]
days=[0,1]
dtt=86400
title='uTSS'
bas='Marmara_basbathy_ser.bas'
simul='uTSS_lobc_chunk_'
param='param_lobc_chunk_'

def setBSUBS(paramfile,day,day0,simul, Nprocs):
    f = open('shympi_chunk_%s.sh'% str(day).zfill(4), 'w')
    f.write('#!/bin/bash\n\n')
    f.write('#SBATCH -A googr\n')
    f.write('#SBATCH -p b224q\n')
    f.write('#SBATCH --time=06:00:00\n')
    f.write('#SBATCH -J uTSS%s\n' % day)
    f.write('#SBATCH -n %s\n' % Nprocs)
    f.write('#SBATCH --output=%j.out\n')
    f.write('#SBATCH --error=%j.err\n')
    if day>day0:
    	f.write('#BSUB -w \"done(uTSS%s)\"\n\n' % (day-1))
    else:
    	f.write('\n')
    f.write('source /okyanus/users/milicak/loadmodule_shyfem \n')
    # f.write('export LD_LIBRARY_PATH=/users/home/sco116/PETSC/2015/petsc-3.7.4/linux-gnu-intel/lib:${LD_LIBRARY_PATH}\n')
    # f.write('export BG_MAXALIGNEXP=-1\n')
    # f.write('export KMP_STACKSIZE=32M\n')
    f.write('export RUNTIME_OPTS=\"-ksp_type bcgs -ksp_rtol 1e-8 -ksp_atol 1e-8\"\n\n')
    # f.write('export I_MPI_FABRICS=shm:ofa\n')
    f.write('time mpirun -np 228 /okyanus/users/milicak/models/shyfem-share/fem3d/shympi %s%s.str $RUNTIME_OPTS > logmsS.log\n\n' % (paramfile,str(day).zfill(4)) )
    f.write('echo \"Job completed at: \" `date`\n\n')
    f.write('sleep 100 \n\n')
    f.write('mv %s%s*.rst restart_files \n\n' %(simul,str(day).zfill(4)))
    f.write('sleep 100 \n\n')
    f.write('mkdir out_%s\n\n' %(str(day).zfill(4)))
    f.write('sleep 100 \n\n')
    f.write('cp nosnc_conv%s.sh  out_%s\n' %(str(day).zfill(3),str(day).zfill(4)))
    f.write('cp ousnc_conv%s.sh  out_%s\n' %(str(day).zfill(3),str(day).zfill(4)))
    f.write('cp wndnc_conv%s.sh  out_%s\n' %(str(day).zfill(3),str(day).zfill(4)))
    f.write('sleep 100 \n\n')
    f.write('mv %s%s*.nos out_%s\n\n' %(simul,str(day).zfill(4),str(day).zfill(4)))
    f.write('sleep 100 \n\n')
    f.write('mv %s%s*.ous out_%s\n\n' %(simul,str(day).zfill(4),str(day).zfill(4)))
    f.write('sleep 100 \n\n')
    f.write('mv %s%s*.wnd out_%s\n\n' %(simul,str(day).zfill(4),str(day).zfill(4)))
    f.write('cd out_%s\n\n' %(str(day).zfill(4)))
    f.write('sleep 100 \n\n')
    f.write('sbatch nosnc_conv%s.sh \n\n' %(str(day).zfill(3)))
    f.write('sbatch ousnc_conv%s.sh \n\n' %(str(day).zfill(3)))
    # f.write(' && ')
    f.write('sbatch wndnc_conv%s.sh' %(str(day).zfill(3)))
    f.close()


for day in range(days[0],days[1]):
    setBSUBS(param,day,days[0],simul, Nprocs=228)
