function [usdt]=ffrt_parfor_OneEp_csearch(filenameepc,modelno,ep,Qa,DOPs,totalNsamp,Nblocks)%Qa,Qb,Qab,Ps
% 0. General Initialization
%%%%%%%%%%%%%%%%Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
prtstr=strcat('entering ffrt_parfor_OneEp...ep=',num2str(ep));
disp(prtstr);
% disp('entering new ffrt_parfor_OneEp...');
%Nblocks=500;%debug
Nsamp=ceil(totalNsamp/Nblocks);%make sure Nsamp is an integer
totalNsamp=Nsamp*Nblocks;% totalNsamp changed slightly due to the ceil function
na      = size(Qa,1);
atrue=randi([-200 200],na,1);% generate random integers
% ncands=2;
[~,~,~,D,~,~] = decorrel(Qa,atrue);
% Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');
counts = struct('na',0,'Psb',0,'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],'falsealarmrt',[]);


% 1. FFRT simulation
% 1.1 Initialization for FFRT
%mus=0.00001:0.00001:1;
mus=1:3;
mus=1./mus';
muslen=length(mus);

ffrtclnmodelno=1; ffrtclnPDOP=2; ffrtclnGDOP=3; ffrtclnADOP=4; 
ffrtclnep=5;ffrtclnna=6;ffrtclnPsb=7;ffrtclnPsILS=8;
ffrtclnmu=9;ffrtclnPfixrt=10;ffrtclnPsrt=11; ffrtclnPfrt=12; ffrtclnPfart=13;

resffrt=zeros(muslen,13);


% 1.2 initialization of the struct array
    counts.Psb=prod(2 * normcdf(1./(2*sqrt(D))) -1 );%scalar
    counts.sucnumILS=0;%scalar
    
    counts.fixnumrt=zeros(muslen,1);%vector (muslen x 1)
    counts.sucnumrt=counts.fixnumrt;%vector (muslen x 1)
    counts.failnumrt=counts.fixnumrt;%vector (muslen x 1)
    counts.falsealarmrt=counts.fixnumrt;%vector (muslen x 1)

filenameblocks=cell(Nblocks,1);
t_be_par=toc;
display(t_be_par);

parfor iblock=1:Nblocks
    filenameblocks{iblock}=strcat(filenameepc,'block',num2str(iblock),'.mat');
    ffrt_OneBlock_csearch(filenameblocks{iblock},Qa,Nsamp);%write the simulation result within the function
end
t_af_par=toc;
display(t_af_par);
for iblock=1:Nblocks
    if ~exist(filenameblocks{iblock},'file')%blockinfos(iblock)==0
        disp(strcat(filenameblocks{iblock},'does not exist!'));
        continue;
    else
        load(filenameblocks{iblock});
        counts=addcounts(counts,blockcounts);
        delete(filenameblocks{iblock});
    end
end
t_af_com=toc;
display(t_af_com);
% parameter for FFRT
    counts.sucnumrt=counts.sucnumrt./counts.fixnumrt; % Pscon
    counts.failnumrt=counts.failnumrt/totalNsamp; % Pf
    counts.fixnumrt=counts.fixnumrt/totalNsamp; % Pfix
    counts.falsealarmrt=counts.falsealarmrt/totalNsamp; % Pfa
    
    resffrt(:,ffrtclnmodelno)=modelno;
    resffrt(:,ffrtclnPDOP)=DOPs(2);
    resffrt(:,ffrtclnGDOP)=DOPs(3);
    resffrt(:,ffrtclnADOP)=DOPs(4);
    resffrt(:,ffrtclnep)=ep;
    resffrt(:,ffrtclnna)=na;
    
    resffrt(:,ffrtclnPsb)=counts.Psb;
    resffrt(:,ffrtclnPsILS)=counts.sucnumILS/totalNsamp;
    resffrt(:,ffrtclnmu)=mus;
    resffrt(:,ffrtclnPfixrt)=counts.fixnumrt;
    resffrt(:,ffrtclnPsrt)=counts.sucnumrt;
    resffrt(:,ffrtclnPfrt)=counts.failnumrt;
    resffrt(:,ffrtclnPfart)=counts.falsealarmrt;
    

% write res for each epoch
% filenameepc=strcat(filenamepre,'_epc',num2str(ep));
printformatstr=makeSformat(size(resffrt,2));
% printformatstrs{epi}=printformatstr;

fid=fopen(filenameepc,'w');
fprintf(fid,printformatstr,resffrt');
fclose(fid);
t_af_writefile=toc;
display(t_af_writefile);
disp('leaving ffrt_parfor_OneEp...');
usdt=toc;
