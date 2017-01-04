
function x= CorrelateGnuRadio(inputBinaryFilefromRIN,outputBinaryFileFromROUT)
s1=csvread(inputBinaryFilefromRIN)
s2=csvread(outputBinaryFileFromROUT)

Fs= 48000

[acor,lag] = xcorr(s2,s1);

[~,I] = max(abs(acor));
lagDiff = lag(I)
timeDiff = lagDiff/Fs

%figure
%plot(lag,acor)
%a3 = gca;
%a3.XTick = sort([-3000:1000:3000 lagDiff]);

x=lagDiff;
y = sprintf('lagDiff=%d',x);
disp(y)
end
