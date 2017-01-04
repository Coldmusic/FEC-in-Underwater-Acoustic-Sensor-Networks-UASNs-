% Dickins seamount test
% mbp
global units jkpsflag

units = 'km';


%%
% Coupled mode runs for range-dependent case

% run Dickins case, using the wedge ocean

D = linspace( 3000, 500, 301 );
preCalcAll(   'DickinsK', D )
%pause( 200 )
fieldLoadAll( 'DickinsK', D )

figure
plotshd( 'DickinsK.mat' )
caxisrev( [ 70 120 ] )
