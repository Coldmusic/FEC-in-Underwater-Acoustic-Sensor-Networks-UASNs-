%recordTime in sec 
function PSKprocessing(i,recordTime,TypeOfModulation)
format long e
%mkdir OFDM_Output;

impulseResponse = strcat('goff_random_80pc_5km_E_no',int2str(i),'.arr');

  
%Reads The PSK.wav file into s1
%[s1,Fs]=audioread('QPSK_output.wav')

%Upsamples s1 and create s2  because delay and sum will downsample this by 2
%s2=upsample(s1,2)

%Writing the upsampled signal to a .wav file
%audiowrite('inputSig.wav',s2,96000,'BitsPerSample',16)
Distance=5
travelDistTime=Distance/1.5
AdjustedrecordTime=recordTime+travelDistTime
%calling the delayandsum- impulseResponse is a string with the name of .arr
%file  for example impulseResponse='goff_random_40pc_2km_E_no1.arr'
[rts,sample_rate]=delayandsum(impulseResponse,AdjustedrecordTime,TypeOfModulation);

direct= strcat('FHFSK_output/0p3SymPerSec/',impulseResponse(1:end-4));
mkdir(direct);
   
for SNR= (-50:5:30)
y=awgn(rts,SNR,'measured');
%normalize rts 
%normalized=(rts-min(rts))/(max(rts)-min(rts));
%time=[0:size(rts)-1]*1/sample_rate;
%plot(time,rts);
%figure;
%plot(time,y);

%savefig('normalized.fig');p

%plot(time,rts);
%savefig('original.fig');
%writing the output signal of delayandsum into a wav file to be processed
%by gnuradio

output=strcat(direct,'/','goff_random_80pc_5km_E_no',int2str(i),'SNR',int2str(SNR),'.wav');

audiowrite(output,y*1000,sample_rate,'BitsPerSample',16)

end


disp('Done!')