% Munk profile test cases
% mbp
global units
units = 'km'

figure
plotssp( 'MunkB_ray' )

bellhop( 'MunkB_ray' )
figure
plotray( 'MunkB_ray' )

bellhop( 'MunkB_eigenray' )
figure
plotray( 'MunkB_eigenray' )

bellhop( 'MunkB_Coh' )
bellhop( 'MunkB_gb' )
kraken(  'MunkK' )
scooter( 'MunkS' )

if ( exist( 'MunkB_Coh.shd', 'file' ) )
    plotshd( 'MunkB_Coh.shd', 2, 2, 1 )
    caxisrev( [ 50 100 ] )
end

if ( exist( 'MunkB_gb.shd', 'file' ) )
    plotshd( 'MunkB_gb.shd', 2, 2, 2 )
    caxisrev( [ 50 100 ] )
end

if ( exist( 'MunkK.shd', 'file' ) )
    plotshd( 'MunkK.shd', 2, 2, 3 )
    caxisrev( [ 50 100 ] )
end

if ( exist( 'MunkS.shd', 'file' ) )
    plotshd( 'MunkS.shd', 2, 2, 4 )
    caxisrev( [ 50 100 ] )
end

plotmode( 'MunkK.mod', [ 1 10 30 60 ] )   % mode selection does not work; fix read_modes_bin

%%

% test of Green's function plotting
figure
plotgrn( 'MunkS.grn' )

% TL slices
figure
plottld( 'MunkS.shd', 50 )
axis( [ 70 120  0 5000 ] )

figure
plottlr( 'MunkS.shd', 800 )
axis( [ 0 100 50 120 ] )
%%

% tests of different Gaussian beam methods

bellhop( 'MunkB_Coh_gb' )
plotshd( 'MunkB_Coh_gb.shd', 2, 2, 1 )
caxisrev( [ 50 100 ] )

bellhop( 'MunkB_Coh_CervenyR' )
plotshd( 'MunkB_Coh_CervenyR.shd', 2, 2, 2 )
caxisrev( [ 50 100 ] )

bellhop( 'MunkB_Coh_CervenyC' )
plotshd( 'MunkB_Coh_CervenyC.shd', 2, 2, 3 )
caxisrev( [ 50 100 ] )

bellhop( 'MunkB_Coh_SGB' )
plotshd( 'MunkB_Coh_SGB.shd', 2, 2, 4 )
caxisrev( [ 50 100 ] )
%%

% tests of incoherent, semi-coherent options

bellhop( 'MunkB_Coh' )
plotshd( 'MunkB_Coh.shd', 3, 2, 1 )
caxisrev( [ 50 100 ] )

bellhop( 'MunkB_Coh_gb' )
plotshd( 'MunkB_Coh_gb.shd', 3, 2, 2 )
caxisrev( [ 50 100 ] )

bellhop( 'MunkB_Semi' )
plotshd( 'MunkB_Semi.shd', 3, 2, 3 )
caxisrev( [ 50 100 ] )

bellhop( 'MunkB_Semi_gb' )
plotshd( 'MunkB_Semi_gb.shd', 3, 2, 4 )
caxisrev( [ 50 100 ] )

bellhop( 'MunkB_Inc' )
plotshd( 'MunkB_Inc.shd', 3, 2, 5 )
caxisrev( [ 50 100 ] )

bellhop( 'MunkB_Inc_gb' )
plotshd( 'MunkB_Inc_gb.shd', 3, 2, 6 )
caxisrev( [ 50 100 ] )

