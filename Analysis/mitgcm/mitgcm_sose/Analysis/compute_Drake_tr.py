import xmitgcm


df = xmitgcm.open_mdsdataset(
                             '/cluster/work/users/milicak/RUNS/mitgcm/mitgcm_sose/Exp01_0s',
                             prefix='UVELMASS',geometry='sphericalpolar')



voltr = df.UVELMASS*df.dyG*df.drF

# new section
Tr = voltr[:,:,160:278,3550]
Trtime = Tr.sum(('Z','YC'))
