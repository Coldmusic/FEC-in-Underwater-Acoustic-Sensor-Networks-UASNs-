#!/bin/bash

arrivalFilename="ice_0pc_5km_no"
receiver_type="resultsfilters64POLY"
NoF=1
#y is the number of arrival files you are trying to process
for ((y=1; y<=$NoF; y++))do

#do the below for each arrival file
for ((SNR=-50; SNR<=30; SNR=SNR+5))
do
sp="_"
outputfile="$arrivalFilename$y"'/'"$receiver_type"'/'"outputBinary$y$sp$SNR"
inputfile="$arrivalFilename$y/$arrivalFilename$y"'SNR'"$SNR.wav"

#extract the new delay
#a=$(grep 'lagDiff=\+' myresult$y$SNR.out | grep  -o "[^=|.]\d\+")
a=$(grep 'lagDiff=\+' "$arrivalFilename$y"'/'"$receiver_type"'/'myresult$y$SNR.out | grep  -o "[-]*[0-9]\+")
echo "THIS IS THE RESULT OF CORRELATION $a"
#calculating BER with the delay
wait ${!}
python2.7 BER.py inputBinary $outputfile $a >>logs/logfile2 &
wait ${!}
echo 'This is the value of '"$SNR"

done
done


