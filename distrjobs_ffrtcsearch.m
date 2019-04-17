function usdt=distrjobs_ffrtcsearch(alljobs,epcs,ntds,hours,mins)
tic
cluster1 = parcluster('PBS');
argum=strcat('-l nodes=',num2str(ntds),':ppn=16,walltime=',num2str(hours),':',num2str(mins),':00');
cluster1.SubmitArguments=argum;%'-l nodes=6:ppn=16,walltime=4:00:00';
%cluster1.SubmitArguments='-l nodes=6:ppn=16,walltime=4:00:00';
S=load('./Haigesimudata/Haige_option.mat');
opts=S.opts;
clear S;

ppnc=16;
% lenepochs=48;
totalNsamp=1E5;
Nblocks=ntds*ppnc-1;
for i=1:length(alljobs)
    num=alljobs(i);
    filename=opts(num).filename;
    loadfile=strcat('./Haigesimudata/',filename,'.mat');
    if ~exist(loadfile,'file')
        display(strcat('File doesnt exist: ',loadfile));
        usdt=0;
        return;
    end
    S=load(loadfile);
    res=S.res; clear S;
    %change the dir name, and the filenamepre
    filenamepre=strcat('./results/','FFRTModel',num2str(num));
    
    for ep=epcs %1:lenepochs
        filenameepc=strcat(filenamepre,'Epoch',num2str(ep),'.txt');
        if ~exist(filenameepc,'file')           
          job = batch(cluster1,@ffrt_parfor_OneEp_csearch,1,{filenameepc,num,ep,res(ep).Qa,res(ep).DOPs,totalNsamp,Nblocks},'Pool',Nblocks);
        end
    end
end
toc
usdt=toc;
