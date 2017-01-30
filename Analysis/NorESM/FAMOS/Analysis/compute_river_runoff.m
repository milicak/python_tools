clear all
rho0 = 1e3 ; % fresh water density
area = ncread('runoff.daitren.iaf.20120419.nc','area');
load('Arctic_mask_woa_grid.mat')
runoff = ncread('runoff.daitren.iaf.20120419.nc','runoff');
lon = ncread('runoff.daitren.iaf.20120419.nc','xc');
lat = ncread('runoff.daitren.iaf.20120419.nc','yc');

area = repmat(area, [1 1 size(runoff,3)]);
maskriver = repmat(maskriver, [1 1 size(runoff,3)]);

% convert from kg/sm2 to Sv [m3/s]
runoffArctic = runoff.*area.*maskriver./rho0;
runoffArctic = squeeze(nansum(runoffArctic,1));
runoffArctic = squeeze(nansum(runoffArctic,1));
runoffArctic = runoffArctic*1e-6; % Sv conversion
runoffArctic = reshape(runoffArctic, [12 64]);

