% Munk profile test cases
% mbp

figure
plotssp( 'MunkB_ray' )

bellhop( 'MunkB_ray' )
figure
plotray( 'MunkB_ray' )

bellhop( 'MunkB_eigenray' )
figure
plotray( 'MunkB_eigenray' )

bellhop( 'MunkB_Coh' )
plotshd( 'MunkB_Coh.shd', 2, 2, 1 )
caxis( [ 50 100 ] )

bellhop( 'MunkB_gb' )
plotshd( 'MunkB_gb.shd', 2, 2, 2 )
caxis( [ 50 100 ] )

kraken( 'MunkK' )
plotshd( 'MunkK.shd', 2, 2, 3 )
caxis( [ 50 100 ] )

scooter( 'MunkS' )
plotshd( 'MunkS.shd', 2, 2, 4 )
caxis( [ 50 100 ] )

figure; plotmode( 'MunkK.mod' )

%%

% tests of different Gaussian beam methods

bellhop( 'MunkB_Coh_gb' )
plotshd( 'MunkB_Coh_gb.shd', 2, 2, 1 )
caxis( [ 50 100 ] )

bellhop( 'MunkB_Coh_CervenyR' )
plotshd( 'MunkB_Coh_CervenyR.shd', 2, 2, 2 )
caxis( [ 50 100 ] )

bellhop( 'MunkB_Coh_CervenyC' )
plotshd( 'MunkB_Coh_CervenyC.shd', 2, 2, 3 )
caxis( [ 50 100 ] )

bellhop( 'MunkB_Coh_SGB' )
plotshd( 'MunkB_Coh_SGB.shd', 2, 2, 4 )
caxis( [ 50 100 ] )
%%

% tests of incoherent, semi-coherent options

bellhop( 'MunkB_Coh' )
plotshd( 'MunkB_Coh.shd', 3, 2, 1 )
caxis( [ 50 100 ] )

bellhop( 'MunkB_Coh_gb' )
plotshd( 'MunkB_Coh_gb.shd', 3, 2, 2 )
caxis( [ 50 100 ] )

bellhop( 'MunkB_Semi' )
plotshd( 'MunkB_Semi.shd', 3, 2, 3 )
caxis( [ 50 100 ] )

bellhop( 'MunkB_Semi_gb' )
plotshd( 'MunkB_Semi_gb.shd', 3, 2, 4 )
caxis( [ 50 100 ] )

bellhop( 'MunkB_Inc' )
plotshd( 'MunkB_Inc.shd', 3, 2, 5 )
caxis( [ 50 100 ] )

bellhop( 'MunkB_Inc_gb' )
plotshd( 'MunkB_Inc_gb.shd', 3, 2, 6 )
caxis( [ 50 100 ] )

% test of Green's function plotting
%figure;
%plotgrn( 'MunkS.grn' )
