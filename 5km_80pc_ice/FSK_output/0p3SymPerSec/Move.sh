#!/bin/bash
for ((c=1; c<=30; c++))
 do
#mkdir goff_random_40pc_5km_E_no$c/filters62POLY
# mv goff_random_40pc_5km_E_no$c/output* goff_random_40pc_5km_E_no$c/filters62POLY
# mv goff_random_40pc_5km_E_no$c/BER goff_random_40pc_5km_E_no$c/filters62Poly
 cp goff_random_80pc_5km_E_no$c/resultsfilters64POLY/* /Volumes/'Seagate Backup Plus Drive'/DataCollected/BPSK_80pc_5km/0.333SymPerSecSenderUsing0p35PulseShapingFilterbw/goff_random_80pc_5km_E_no$c/resultsfilters64POLY/

done
