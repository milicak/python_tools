% this subroutine computes masks for Beaufort Gyre and Greenland Sea
% provided by Jeff Scott
% Beaufort Gyre: 130W-170W, 70.5-80.5N
% Greenland Sea: 345-5E, 72.0-77.0N
clear all

% load 1 degree tripolar grid lon lat
lonmicom

% Beaufort Gyre mask
lon1 = [-130 -130 -170 -170];
lat1 = [70.5 80.5 80.5 70.5];
lon1(end+1) = lon1(1);
lat1(end+1) = lat1(1);
in = insphpoly(lon,lat,lon1,lat1,0.,90.);
in = double(in);

% let's plot the lon1 lat1
figure
m_projpolar
hold on
m_pcolor(lon,lat,in);shfn
m_plot(lon1,lat1,'r-')
m_coast('patch',[.7 .7 .7])

BG_mask = in;

% Greenland Sea mask
lon1 = [345.0 345.0 5.0 5.0];
lat1 = [72.0 77.0 77.0 72.0];
lon1(end+1) = lon1(1);
lat1(end+1) = lat1(1);
in = insphpoly(lon,lat,lon1,lat1,0.,90.);
in = double(in);
% let's plot the lon1 lat1
figure
m_projpolar
hold on
m_pcolor(lon,lat,in);shfn
m_plot(lon1,lat1,'r-')
m_coast('patch',[.7 .7 .7])

GS_mask = in;

save('BG_GS_masks.mat','BG_mask','GS_mask');
