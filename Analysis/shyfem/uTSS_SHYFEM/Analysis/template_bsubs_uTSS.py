date=20200101
days=[0,366]
# days=[0,1096]
#days=[366,732]
#days=[0,31]
dtt=86400
title='uTSS'
bas='Marmara_basbathy_ser.bas'
simul='uTSS_lobc_chunk_'
param='param_lobc_chunk_'

def setBSUBS(paramfile,day,day0,simul, Nprocs):

    f = open('inp_files/shympi_chunk_%s.sh'% str(day).zfill(4), 'w')
    f.write('#!/bin/bash\n\n')
    f.write('#BSUB -x\n')
    f.write('#BSUB -q p_short\n')
    f.write('#BSUB -P 0285\n')
    # f.write('#BSUB -q p_medium\n')
    f.write('#BSUB -J uTSS_SHYFEM%s\n' % day)
    f.write('#BSUB -n %s\n' % Nprocs)
    f.write('#BSUB -o logout.%J.out\n' )
    f.write('#BSUB -e logerr.%J.err\n' )
    #f.write('#BSUB -sla SC_sanifs_longrun_PS4\n')
    if day>day0:
        f.write('#BSUB -w \"done(uTSS_SHYFEM%s)\"\n\n' % (day-1))
    else:
        f.write('\n')
    f.write('module purge\n')
    f.write('module purge\n')
    f.write('module load intel19.5/19.5.281\n')
    f.write('module load impi19.5/19.5.281\n')
    f.write('module load impi19.5/petsc/3.7.5\n\n')
    f.write('module load intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1\n\n')
    f.write('export RUNTIME_OPTS=\"-ksp_type bcgs -ksp_rtol 1e-12 -ksp_atol 1e-12\"\n\n')
    f.write('export FI_PROVIDER="verbs;ofi_rxm"\n\n')
    f.write('echo "FI_PROVIDER: ${FI_PROVIDER}"\n\n')
    f.write('if [ "$I_MPI_HYDRA_BOOTSTRAP" == "" ]; then\n')
    f.write('  export export I_MPI_HYDRA_BOOTSTRAP=lsf\n')
    f.write('fi\n')
    f.write('export I_MPI_HYDRA_BRANCH_COUNT=7\n')
    f.write('if [ "$I_MPI_HYDRA_BRANCH_COUNT" == "" ]; then\n')
    f.write('  export I_MPI_HYDRA_BRANCH_COUNT=`cat $LSB_DJOB_HOSTFILE | uniq | wc -l`\n')
    f.write('fi\n\n')
    f.write('if [ "$I_MPI_LSF_USE_COLLECTIVE_LAUNCH" == "" ]; then\n')
    f.write('  export I_MPI_LSF_USE_COLLECTIVE_LAUNCH=1\n')
    f.write('fi\n\n')
    f.write('cd /work/opa/mi19918/Projects/uTSS_SHYFEM/work\n\n')
    f.write('time mpiexec.hydra -l /users_home/opa/mi19918/models/SHYFEM/bin/shympi %s%s.str $RUNTIME_OPTS > logmsS.log\n\n' % (paramfile,str(day).zfill(4)) )
    f.write('echo \"Job completed at: \" `date`\n\n')
    f.write('sleep 100 \n\n') 
    f.write('mv %s%s*.rst restart_files \n\n' %(simul,str(day).zfill(4))) 
    f.write('sleep 100 \n\n') 
    f.write('mkdir out_%s\n\n' %(str(day).zfill(4))) 
    f.write('sleep 100 \n\n') 
    f.write('cp nosnc_conv%s.sh  out_%s\n' %(str(day).zfill(3),str(day).zfill(4))) 
    f.write('cp ousnc_conv%s.sh  out_%s\n' %(str(day).zfill(3),str(day).zfill(4))) 
    # f.write('cp wndnc_conv%s.sh  out_%s\n' %(str(day).zfill(3),str(day).zfill(4))) 
    f.write('sleep 100 \n\n') 
    f.write('mv %s%s*.nos out_%s\n\n' %(simul,str(day).zfill(4),str(day).zfill(4))) 
    f.write('sleep 100 \n\n') 
    f.write('mv %s%s*.ous out_%s\n\n' %(simul,str(day).zfill(4),str(day).zfill(4))) 
    f.write('sleep 100 \n\n') 
    # f.write('mv %s%s*.wnd out_%s\n\n' %(simul,str(day).zfill(4),str(day).zfill(4))) 
    f.write('cd out_%s\n\n' %(str(day).zfill(4))) 
    f.write('sleep 100 \n\n') 
    f.write('bsub < nosnc_conv%s.sh \n\n' %(str(day).zfill(3)))
    f.write('bsub < ousnc_conv%s.sh \n\n' %(str(day).zfill(3)))
    # f.write('bsub < wndnc_conv%s.sh' %(str(day).zfill(3)))
    f.close()


for day in range(days[0],days[1]):
    setBSUBS(param,day,days[0],simul, Nprocs=228)
