
%Calculating the RMS of the delay spread 
% formula of mean delay spread is m= (sum{(amp(t)^2)(t)}/sum{amp(t)^2}
%formula of rms delay spread is rms= sqrt(n-m^2)
%n= (sum{(amp(t)^2)(t)^2}/sum{amp(t)^2}
%calculate and write the delay spread to a file
%we calculate the rms of each file then we get the mean of those.
size=100;
meanDsp=zeros(size,1);
n=zeros(size,1);
rmsDsp=zeros(size,1);
for i=1:size
 ARRFIL{i}=strcat('goff_random_80pc_5km_E_no',int2str(i),'.arr');
 [ Arr(i), Pos(i) ] = read_arrivals_asc( char(ARRFIL(i)), 100 );  % read the arrivals file
end
for i=1:size

 sum1=0;
 sum2=0;
 sum3=0;
 for j=1:Arr(i).Narr
  
     sum1=sum1+(abs(Arr(i).A(j))^2)*(Arr(i).delay(j));
     sum2=sum2+(abs(Arr(i).A(j))^2);
     sum3=sum3+(abs(Arr(i).A(j))^2)*((Arr(i).delay(j))^2);
 end
 meanDsp(i)=sum1/sum2;
 n(i)=sum3/sum2;
 rmsDsp(i)=sqrt((n(i)-((meanDsp(i))^2)));
 
end

stem(rmsDsp);


