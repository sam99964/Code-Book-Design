clear all;clc;
tic
load('LDPC_BPSK_simple.mat', 'Deintlv');
Deintlv=Deintlv(1,:);

% use Deintlv to generate Edgeintlv
Edgeintlv=zeros(1,length(Deintlv),'int32');
for m=1:length(Edgeintlv)
    Edgeintlv(Deintlv(m))=m; 
end

% check interleaver
Edgesorted=sort(Edgeintlv);
tmp=1:length(Deintlv);
if(Edgesorted ~= tmp )
    sprintf('error interleaver ? \n')  
end    

% transform these data into txt file
fid = fopen('Cq_Deintlv.txt','wt');
fprintf(fid,'%d\n',Deintlv-1);   % -1 for index start from 0
fclose(fid);

fid = fopen('Cq_intlv.txt','wt');
fprintf(fid,'%d\n',Edgeintlv-1);   % -1 for index start from 0
fclose(fid);
toc