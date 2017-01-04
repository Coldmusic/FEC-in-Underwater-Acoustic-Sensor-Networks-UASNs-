function runtests
% Comparison of range-independent, adiabatic, and coupled mode options for
% the Gulf Stream eddy

% Note that in the new version of KRAKEN, this can all be done in one input
% file.

% Possible compiler problem showed up on this test case
% A segmentation fault in evalad.f (adiabatic field evaluation) occured when the number of receiver depths was
% increased
global units jkpsflag
units = 'km';

runkraken = which( 'kraken.exe' );
runfield  = which( 'field.exe' );

% run kraken for each of the profiles
eval( [ '! "' runkraken '" gulf_rd' ] );

% do the field calculations
copyfile( 'gulf_rd.mod', 'gulf.mod' )
eval( [ '! "' runfield '" gulf < gulf_ri.flp > field.prt' ] );
movefile( 'gulf.shd', 'gulf_ri.shd' )
plotshd( 'gulf_ri.shd', 3, 1, 1 )
caxisrev( [ 70 100 ] )

eval( [ '! "' runfield '" gulf < gulf_cm.flp > field.prt' ] );
movefile( 'gulf.shd', 'gulf_cm.shd' )
plotshd( 'gulf_cm.shd', 3, 1, 2 )
caxisrev( [ 70 100 ] )

eval( [ '! "' runfield '" gulf < gulf_ad.flp > field.prt' ] );
movefile( 'gulf.shd', 'gulf_ad.shd' )
plotshd( 'gulf_ad.shd', 3, 1, 3 )
caxisrev( [ 70 100 ] )

if ( jkpsflag )
    print -deps2 Fig5_18.eps
end
%%
figure
plotssp 'Gulf_ray_ri'

bellhop 'Gulf_ray_ri'
figure
plotray 'Gulf_ray_ri'

figure
plotssp2D 'Gulf_ray_rd'

bellhop 'Gulf_ray_rd'
figure
plotray 'Gulf_ray_rd'
%%

bellhop 'GulfB_ri_geo'
plotshd( 'GulfB_ri_geo.shd', 2, 2, 1 )
caxisrev( [ 70 100 ] )

bellhop 'GulfB_ri_gb'
plotshd( 'GulfB_ri_gb.shd', 2, 2, 2 )
caxisrev( [ 70 100 ] )

bellhop 'GulfB_rd_geo'
plotshd( 'GulfB_rd_geo.shd', 2, 2, 3 )
caxisrev( [ 70 100 ] )

bellhop 'GulfB_rd_gb'
plotshd( 'GulfB_rd_gb.shd', 2, 2, 4 )
caxisrev( [ 70 100 ] )

