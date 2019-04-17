function [t_lv_block]=ffrt_OneBlock_csearch(filenameblock,Qa,Nsamp)%Qa,Qb,Qab,Ps
% 0. General Initialization
%%%%%%%%%%%%%%%%Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
disp('entering new ffrt_OneBlock...');
% info=0;
na      = size(Qa,1);
atrue=randi([-200 200],na,1);% generate random integers
ncands=2;
[Qzhat,~,L,D,ztrue,~] = decorrel(Qa,atrue);
Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');

% blockcounts = repmat(struct('Psb',[],'ratio',[],'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],'falsealarmrt',[],...
%     'themus',[],'thePfix',[],'thePfs',[],'thePscons',[],'theFalseAlarm',[]), na, 1 );
% 1. FFRT simulation
% 1.1 Initialization for FFRT
% mus=0.00001:0.00001:1;
mus=1:3;
mus=1./mus';
muslen=length(mus);

Psb=0;
SucnumILSs=Psb;
fixnumrts=zeros(muslen,1);
sucnumrts=fixnumrts;
falsealarmrts=fixnumrts;

% 1.2 initialization of the struct array
    Psb=prod(2 * normcdf(1./(2*sqrt(D))) -1 );

% fileprogressname=strcat(filenameblock,'progress.txt');
% fid=fopen(fileprogressname,'w');
% 1.3 Simulate all the samples and do LAMBDA, ratio test.
%     zhats=mvnrnd(ztrue,Qzhat,Nsamp)';%10^7 samples, used up 1.28GB memory
% load('myzhats.mat');
for i=1:Nsamp
%     if ~mod(i,1E5)
%         str=strcat('In ',filenameblock,' Progress ',num2str(i));
%         disp(str);
%     end
    %     zhat=zhats(:,i);
    zhat=mvnrnd(ztrue,Qzhat,1)';
        [zpar,sqnorm] = csearchL(zhat,L,D,ncands);
        zpar=zpar(:,1);
        
        sucflag=sum(~(ztrue==zpar))==0;
        
        SucnumILSs = SucnumILSs+sucflag;         % Correctly fixed
        ratio=sqnorm(1)/sqnorm(2);
        rtpass=ratio<mus;
        
        fixnumrts=fixnumrts+rtpass;          % fixed, but correctness unknown
        sucnumrts=sucnumrts+ double(rtpass&sucflag); % fixed and correct
        falsealarmrts=falsealarmrts+ double((~rtpass)&sucflag); % rejected but correct
    
end
failnumrts=fixnumrts-sucnumrts;

blockcounts = struct('na',0,'Psb',0,'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],'falsealarmrt',[] );
    blockcounts.na=na;
    blockcounts.Psb=Psb;
    blockcounts.sucnumILS=SucnumILSs;
    blockcounts.fixnumrt=fixnumrts;
    blockcounts.sucnumrt=sucnumrts;
    blockcounts.falsealarmrt=falsealarmrts;
    blockcounts.failnumrt=failnumrts;
save(filenameblock,'blockcounts');
% fclose(fid);
% delete(fileprogressname);
t_lv_block=toc;
display(t_lv_block);
disp('leaving new ffrt_OneBlock...');
% info=1;
