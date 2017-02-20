% 
clear all
resolution = 'quarterdegree'

switch resolution
    case 'onedegree'
        lonmicom
        mask = ncgetvar('/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc','pmask');
        filename = 'noresm_ESMF_grid_tnx1v1.nc';
    case 'quarterdegree'
        lonmicom2
        mask = ncgetvar('/hexagon/work/shared/noresm/inputdata/ocn/micom/tnx0.25v1/20130930/grid.nc','pmask');
        filename = 'noresm_ESMF_grid_tnx0.25v1.nc';
end

grid_size = size(pclon,2);
grid_corners = 4;
grid_rank = 2;

% create the netcdf file
ncid=netcdf.create(filename,'NC_CLOBBER');
% Define dimensions.
ni_dimid=netcdf.defDim(ncid,'grid_size',grid_size);
nj_dimid=netcdf.defDim(ncid,'grid_corners',grid_corners);
nz_dimid=netcdf.defDim(ncid,'grid_rank',grid_rank);

tdim_varid=netcdf.defVar(ncid,'grid_dims','nc_int',[grid_rank]);
netcdf.putAtt(ncid,tdim_varid,'long_name','grid_dims');

tlat_varid=netcdf.defVar(ncid,'grid_center_lat','double',[ni_dimid]);
netcdf.putAtt(ncid,tlat_varid,'long_name','grid_center_lat');
netcdf.putAtt(ncid,tlat_varid,'units','degrees');

tlon_varid=netcdf.defVar(ncid,'grid_center_lon','double',[ni_dimid]);
netcdf.putAtt(ncid,tlon_varid,'long_name','grid_center_lon');
netcdf.putAtt(ncid,tlon_varid,'units','degrees');

tmask_varid=netcdf.defVar(ncid,'grid_imask','nc_int',[ni_dimid]);
netcdf.putAtt(ncid,tmask_varid,'long_name','grid_imask');
netcdf.putAtt(ncid,tmask_varid,'units','unitless');

tclat_varid=netcdf.defVar(ncid,'grid_corner_lat','double',[ni_dimid nj_dimid]);
netcdf.putAtt(ncid,tclat_varid,'long_name','grid_corner_lat');
netcdf.putAtt(ncid,tclat_varid,'units','degrees');

tclon_varid=netcdf.defVar(ncid,'grid_corner_lon','double',[ni_dimid nj_dimid]);
netcdf.putAtt(ncid,tclon_varid,'long_name','grid_corner_lon');
netcdf.putAtt(ncid,tclon_varid,'units','degrees');

%lont_bounds_varid=netcdf.defVar(ncid,'lon_bnds','float',[nz_dimid ni_dimid]);
%netcdf.putAtt(ncid,lont_bounds_varid,'long_name','longitude boundaries of T cells');
%netcdf.putAtt(ncid,lont_bounds_varid,'units','degrees_east');

%latt_bounds_varid=netcdf.defVar(ncid,'lat_bds','float',[nz_dimid nj_dimid]);
%netcdf.putAtt(ncid,latt_bounds_varid,'long_name','latitude boundaries of T cells');
%netcdf.putAtt(ncid,latt_bounds_varid,'units','degrees_north');

% End definitions and leave define mode.
netcdf.endDef(ncid)

nx = size(lon,1);
ny = size(lon,2);

grank = [nx,ny];

% Provide values for time invariant variables.
netcdf.putVar(ncid,tlon_varid,lon(:));
netcdf.putVar(ncid,tlat_varid,lat(:));
netcdf.putVar(ncid,tmask_varid,mask(:));
netcdf.putVar(ncid,tdim_varid,grank);
netcdf.putVar(ncid,tclat_varid,pclat);
netcdf.putVar(ncid,tclon_varid,pclon);

% Close netcdf file
netcdf.close(ncid)
