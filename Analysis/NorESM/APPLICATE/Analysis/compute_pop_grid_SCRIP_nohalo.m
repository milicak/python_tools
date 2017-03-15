% 
clear all
resolution = 'onedegree'
%resolution = 'quarterdegree'

switch resolution
    case 'onedegree'
        gridfile = '/export/grunchfs/unibjerknes/milicak/bckup/noresm/CORE2/Arctic/DATA/ncar-pop/areacello_fx_CCSM4_piControl_r0i0p0.nc';
        lon = ncread(gridfile,'lon');
        lat = ncread(gridfile,'lat');
        pclon = ncread(gridfile,'lon_vertices');
        pclat = ncread(gridfile,'lat_vertices');
        maskfile = '/export/grunchfs/unibjerknes/milicak/bckup/noresm/CORE2/Arctic/DATA/ncar-pop/CESM_Arctic_mask.nc';
        mask = ncread(maskfile,'Arctic_mask');
        mask(mask~=0)=1;
        filename = 'ncar_pop_ESMF_grid_gx1v6_nohalo.nc';
end
pclon = reshape(pclon,[size(pclon,1) size(pclon,2)*size(pclon,3)]);
pclat = reshape(pclat,[size(pclon,1) size(pclon,2)*size(pclon,3)]);

nx = size(lon,1);
ny = size(lon,2);
grank = [nx,ny];
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


% Provide values for time invariant variables.
netcdf.putVar(ncid,tlon_varid,lon(:));
netcdf.putVar(ncid,tlat_varid,lat(:));
netcdf.putVar(ncid,tmask_varid,mask(:));
netcdf.putVar(ncid,tdim_varid,grank);
netcdf.putVar(ncid,tclat_varid,pclat);
netcdf.putVar(ncid,tclon_varid,pclon);

% Close netcdf file
netcdf.close(ncid)
