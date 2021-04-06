
# Web pages that are useful
# https://github.com/nikizadehgfdl/grid_generation/wiki/Generating-the-Gridspec-file-(Mosaic-Grid)
# https://github.com/nikizadehgfdl/grid_generation/wiki/Regrid-runoff-and-salt-restore-files


/okyanus/users/milicak/models/FRE-NCtools/bin/make_solo_mosaic --num_tiles 1 --dir . --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc --periodx 360.

/okyanus/users/milicak/models/FRE-NCtools//tools/make_quick_mosaic/make_quick_mosaic --input_mosaic ocean_mosaic.nc --mosaic_name grid_spec --ocean_topog ocean_topog.nc

ncks -C -O -x -v gridtiles ocean_mosaic.nc dnm.nc
ncks -v gridtiles ~/dataset/MOM6/OM4_05/mosaic_ocean.v20180227.unpacked/ocean_mosaic.nc dnm2.nc
ncks -A dnm2.nc dnm.nc
mv dnm.nc ocean_mosaic.nc 

cd WHERE_MOSAIC_GRID_IS
# FOR RIVER RUNOFF
../../../python_tools/Analysis/mom6/Arctic/Analysis/regrid_runoff/regrid_runoff.py --fast_pickle ocean_hgrid.nc ocean_mask.nc ~/dataset/CORE2/NYF_v2.0/runoff.daitren.clim.v2011.02.10.nc --fms runoff.daitren.clim.nc

../../../python_tools/Analysis/mom6/Arctic/Analysis/regrid_runoff/regrid_runoff.py --fast_pickle ocean_hgrid.nc ocean_mask.nc ~/dataset/CORE2/NYF_v2.0/runoff.daitren.iaf.20120419.nc --fms runoff.daitren.iaf.nc

# for SEAWIFS
../../../python_tools/Analysis/mom6/Arctic/Analysis/OM4_025_seawifs/interp_and_fill/interp_and_fill.py ocean_hgrid.nc ocean_mask.nc ~/dataset/CORE2/NYF_v2.0/seawifs-clim-1997-2010.nc chlor_a --fms seawifs-clim-1997-2010.smoothed.nc
/okyanus/users/milicak/python_tools/Analysis/mom6/Arctic/Analysis/OM4_025_seawifs/interp_and_fill/interp_and_fill.py ocean_hgrid.nc ocean_mask.nc ~/dataset/CORE2/NYF_v2.0/seawifs-clim-1997-2010.nc chlor_a --fms seawifs-clim-1997-2010.smoothed.nc

# for SALT RESTORE
../../../python_tools/Analysis/mom6/Arctic/Analysis/OM4_025_seawifs/interp_and_fill/interp_and_fill.py ocean_hgrid.nc ocean_mask.nc ~/dataset/CORE2/NYF_v2.0/PHC2_salx.2004_08_03.corrected.nc SALT --fms --closest salt_restore_PHC2.nc

# for GEOTHERMAL HEATING
python ../../../python_tools/Analysis/mom6/Arctic/Analysis/OM4_025_preprocessing_geothermal/regrid_geothermal.py
