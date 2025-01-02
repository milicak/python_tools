import numpy as np
import scipy.io

# FNC1:string - FNC1:dz_min,H_total,power,precision

# fortran routine follows as
#  !> Parses a string and generates a dz(:) profile that goes like k**power.
#  subroutine dz_function1( string, dz )
#    character(len=*),   intent(in)    :: string !< String with list of parameters in form
#                                                !! dz_min, H_total, power, precision
#    real, dimension(:), intent(inout) :: dz     !< Profile of nominal thicknesses
#    ! Local variables
#    integer :: nk, k
#    real    :: dz_min, power, prec, H_total
#
#    nk = size(dz) ! Number of cells
#    prec = -1024.
#    read( string, *) dz_min, H_total, power, prec
#    if (prec == -1024.) call MOM_error(FATAL,"dz_function1: "// &
#            "Problem reading FNC1: string  ="//trim(string))
#    ! Create profile of ( dz - dz_min )
#    do k = 1, nk
#      dz(k) = (real(k-1)/real(nk-1))**power
#    enddo
#    dz(:) = ( H_total - real(nk) * dz_min ) * ( dz(:) / sum(dz) ) ! Rescale to so total is H_total
#    dz(:) = anint( dz(:) / prec ) * prec ! Rounds to precision prec
#    dz(:) = ( H_total - real(nk) * dz_min ) * ( dz(:) / sum(dz) ) ! Rescale to so total is H_total
#    dz(:) = anint( dz(:) / prec ) * prec ! Rounds to precision prec
#    dz(nk) = dz(nk) + ( H_total - sum( dz(:) + dz_min ) ) ! Adjust bottommost layer
#    dz(:) = anint( dz(:) / prec ) * prec ! Rounds to precision prec
#    dz(:) = dz(:) + dz_min ! Finally add in the constant dz_min
#
#  end subroutine dz_function1

dz_min = 2
H_total = 4000
power = 4.5
precision = 0.01

dz_min = 5
H_total = 8000
power = 1
precision = 0.01

# number of cells
nk = 75
prec = precision

# Create profile of ( dz - dz_min )
dz = np.zeros(nk)
for k in range(0,nk):
    dz[k] = (np.real(k)/np.real(nk-1))**power

#   Rescale to so total is H_total
dz = ( H_total - np.real(nk) * dz_min ) * ( dz / np.sum(dz) )
dz = np.round(dz/prec)*prec
dz = ( H_total - np.real(nk) * dz_min ) * ( dz / np.sum(dz) )
dz = np.round(dz/prec)*prec
dz = ( H_total - np.real(nk) * dz_min ) * ( dz / np.sum(dz) )
dz = np.round(dz/prec)*prec

dz += dz_min

deltaz = H_total-dz.sum()
dz[-1] += deltaz

# Create a mosaic file
fout  = 'vgrid.nc'
rg = scipy.io.netcdf_file(fout,'w')
rg.createDimension('nz',nk)
hz = rg.createVariable('dz','double',('nz',))
hz.units = 'm'
hz.long_name = 'z coordinate level thickness'
hz[:] = np.copy(dz)
rg.close()

