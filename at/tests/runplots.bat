REM test plotting routines

call plotssp  framiv

copy fieldarc.flp field.flp
call kraken   framiv
call plottlr  framiv
call plottld  framiv
call plotmode framiv
call plotgrn  framiv
REM call field    framiv

call bounce   refl
call plotrth  refl

call sparc    iso
call plotts   iso

call bellhop  munkB
call plotray  munkB
call bellhop  calibB

call plottri med

REM test conversion routines

call toasc framiv
call tobin framiv

modasc framiv
modbin framiv


