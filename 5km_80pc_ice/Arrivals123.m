%Arrivals

size=100;
ARRFIL=cell(size:1);
for i=1:size
 ARRFIL{i}=strcat('goff_random_80pc_5km_E_no',int2str(i),'.arr');
 %disp(ARRFIL)
end
max=0;
for i=1:size
[ Arr(i), Pos(i) ] = read_arrivals_asc( char(ARRFIL(i)), 100 );  % read the arrivals file
% Finds the max number of arrivals in the arrival files read
nmax=Arr(i).Narr;
  if max<nmax
    max= nmax;
  end 
end 

for i=1:size
    
A(i)=Arr(i).Narr
    
end
