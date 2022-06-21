% We developed a Matlab package, to rapidly obtain the apparent EArthquake 
% SOurce time function using the Average of near-source P-wave content of 
% the observations (EASOt-AP) and estimate the earthquake source parameters, 
% such as the seismic moment, the rupture radius, and the average static stress drop. 
% The algorithm implemented in this package is based on a rapid and straightforward 
% methodology that is recently developed by Zollo et al. (2021). 
% To this purpose, EASOt-AP retrieves and models the time evolution of the 
% average P-wave displacement in the logarithm scale named as LPDT curve. 
% LPDT Method to charactrize the earthquake source M0/Mw, Stress Drop and Source radius
% Sahar Nazeri, 2020, Naples, Italy
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clc;clear;
path0=pwd;
% =============================================================Don't change
%%
warning off
tic
global path_Val path_Reg
global  Vp Vs 

addpath(genpath([path0,'\CODES']))
cd (path0)
% ================================================= Input Parameters
A=textread([path0 '\INPUT\','Input.txt'],'%s','headerlines',5,'delimiter','\n');
Region=A{2};
wave_type=A{4};
dirparam=A{6};
clear A0; A0=split(A{10}, ',');
Unit=A0{1};PU=str2num(A0{2});
minSta=str2num(A{12});
clear A0; A0=split(A{15}, ',');
Rm1=str2num(A0{1});Rm3=str2num(A0{2});Rm5=str2num(A0{3});RmM=str2num(A0{4});
clear A0; A0=split(A{17}, ',');
SNR=str2num(A0{1});
clear A0; A0=split(A{21}, ',');
nPol=str2num(A0{1});cor1=str2num(A0{2});cor2=str2num(A0{3});
clear A0; A0=split(A{24}, ',');
Vfilter=A0{1};nPolV=str2num(A0{2});cor1V=str2num(A0{3});cor2V=str2num(A0{4});
clear A0; A0=split(A{29}, ',');
WinFix=A0{1};WinMax=str2num(A0{2});S_P_Coef=str2num(A0{3});
clear A0; A0=split(A{31}, ',');
WinFix_fit=A0{1};WinMax_fit=str2num(A0{2});
clear A0; A0=split(A{35}, ',');
rho=str2num(A0{1});Vp=str2num(A0{2});Vs=str2num(A0{3});SD=str2num(A0{4});
RadP=str2num(A{37});
clear A0; A0=split(A{41}, ',');
QFilter=A0{1};Q=str2num(A0{2});
FitFunc=A{43};
clear A0; A0=split(A{47}, ',');
SD1 =str2num(A0{1}); SD2 =str2num(A0{2}); SDN = str2num(A0{3});
clear A A0
% =========================================================================
% Read the magnitude information from a list
% [eve, M] = textread([path0 '\INPUT\','maginfo-' Region '.txt'],'%s %f');
% M=round(M.*10)/10;
...........................................................................
% Path Definiations
path_sac=[path0,'\INPUT\',Region];
mkdir(path0,'OUTPUT');

path_out=[path0,'\OUTPUT'];
mkdir([path0,'\OUTPUT\',Region]);
path_out=[path0,'\OUTPUT\',Region];

mkdir(path_out,'SignalProcessing');

path_Fig=[path_out,'\SignalProcessing'];
path_Val=[path_out,'\Values'];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PeaksAndPlots;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (size(coefE,2)~=1)
    mkdir(path_out,'Regression');
    path_Reg=[path_out,'\Regression'];
end
...........................................................................
...........................................................................
output = fullfile(path_out,['Result_',Region,'.xlsx']);
...........................................................................
...........................................................................
[~,~,logM0c,logM0t_M,d_logM0t_M,DeltaS_M,dDeltaS_M,loga_M,dloga_M,TX_M,dTX_M]=App2('m',Mweve,PL.MainAmp(:,4),...
PL.MainTime(:,4)-PL.MainTime(:,1),Depth,coefE(9,:),coefE(10,:));
App2M.event = EveL;
App2M.logM0_Cata=round(100.*logM0c)./100;
App2M.Mw_Cata=(App2M.logM0_Cata-9.1)./1.5;
App2M.logM0_LPDT=round(100.*logM0t_M)./100;
App2M.d_logM0_LPDT=round(1000.*d_logM0t_M)./1000;
App2M.Mw_LPDT=(App2M.logM0_LPDT-9.1)./1.5;
App2M.d_Mw_LPDT=(App2M.d_logM0_LPDT)./1.5;
App2M.StressDrop_MPa=round(10000.*DeltaS_M)./10000;
App2M.d_StressDrop=round(1000.*dDeltaS_M)./1000;
App2M.Radius_Km=round(10.*((10.^loga_M)/1000))./10;
App2M.d_Radius=round(1000.*((10.^dloga_M)/1000))./1000;
App2M.Duration_s=round(1000.*TX_M)./1000;
App2M.d_Duration=round(1000.*dTX_M)./1000;

clear T2; T2 = struct2table(App2M);
% writetable(T2,fullfile(path_Val,'App2.txt'),'Delimiter','\t','WriteRowNames',false);
writetable(T2,output,'WriteRowNames',false,'Sheet','App2_Mode');
...........................................................................
...........................................................................
if (strcmp(QFilter,'on')==1)
    AppQ.Mode.event = EveL;
    AppQ.Mode.logM0_Cata=round(100.*logM0c)./100;
    AppQ.Mode.Mw_Cata=(App2M.logM0_Cata-9.1)./1.5;
    AppQ.Mode.logM0_LPDT=round(100.*(1.5*QCorI_M(:,1)+9.1))./100;
    AppQ.Mode.d_logM0_LPDT=round(1000.*d_logM0t_M)./1000;
    AppQ.Mode.Mw_LPDT=round(100.*QCorI_M(:,1))./100;
    AppQ.Mode.d_Mw_LPDT=(AppQ.Mode.d_logM0_LPDT)./1.5;
    AppQ.Mode.StressDrop_MPa=round(10000.*QCorI_M(:,2))./10000;
    AppQ.Mode.d_StressDrop=round(10000.*dDeltaS_M)./10000;
    AppQ.Mode.Radius_Km=round(10.*((10.^QCorI_M(:,6))/1000))./10;
    AppQ.Mode.d_Radius=round(1000.*((10.^dloga_M)/1000))./1000;
    AppQ.Mode.Duration_s=round(100.*2.*QCorI_M(:,3))./100;
    AppQ.Mode.d_Duration=round(1000.*dTX_M)./1000;        
    T3 = struct2table(AppQ.Mode);
%     writetable(T3,fullfile(path_Val,'Qfilter.txt'),'Delimiter','\t','WriteRowNames',false);
    writetable(T3,output,'WriteRowNames',false,'Sheet','App2_Qfilter_Mode');
    ...................................................................PLOT
    if (size(coefE,2)~=1)
    subplot 221; 
    hold on
    scatter(logM0c,AppQ.Mode.logM0_LPDT,30,[0.7 0.7 0.7],'fill','MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
    subplot 222; 
    hold on
    scatter(AppQ.Mode.logM0_LPDT,log10(AppQ.Mode.Duration_s),30,[0.7 0.7 0.7],'fill','MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
    subplot 223; 
    hold on;
    scatter(AppQ.Mode.logM0_LPDT,log10(AppQ.Mode.Radius_Km.*1000),30,[0.7 0.7 0.7],'fill','MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
    subplot 224; 
    hold on
    h=histfit(log10(AppQ.Mode.StressDrop_MPa));   
    h(1).FaceColor = [0.7 0.7 0.7]; 
    h(2).Color = [.6 .6 .6];
    h(1).FaceAlpha=0.7;
    fsd=find((h(2).YData>= max(h(2).YData)));
    maxnorm =10.^(h(2).XData(fsd(1)));
    xline(h(2).XData(fsd(1)),'Linewidth',2,'color',[.6 .6 .6])
%     xt = xticklabels; 
%     xnum=10.^( cellfun(@str2num,xt) );
%     xticklabels(round(1000.*xnum)./1000)
    title(['Stress Drop= ',num2str(round(1000.*maxnorm)./1000),' MPa'],'FontSize',9)    
    end
end
...........................................................................
if (size(coefE,2)~=1)
clear Mo_test M0_test_min M0_test_max; 
if (strcmp(QFilter,'on')==1)
M0_test_min=min([AppQ.Mode.logM0_LPDT App2M.logM0_LPDT],[],'all');
M0_test_max=max([AppQ.Mode.logM0_LPDT App2M.logM0_LPDT],[],'all');
else
M0_test_min=min(App2M.logM0_LPDT);
M0_test_max=max(App2M.logM0_LPDT);
end
M0_test=M0_test_min:((M0_test_max-M0_test_min)/10):M0_test_max;
a1_01=log10((7.*(10.^M0_test)./(16*0.01*1e6)).^(1./3));
a1_1=log10((7.*(10.^M0_test)./(16*0.1*1e6)).^(1./3));
a1=log10((7.*(10.^M0_test)./(16*1*1e6)).^(1./3));
a10=log10((7.*(10.^M0_test)./(16*10*1e6)).^(1./3));
subplot 221;hold on;box on;grid on
subplot 222;hold on;box on;grid on
subplot 224;hold on;box on;grid on
subplot 223;hold on;box on;grid on
plot(M0_test,a1_01,':k','Linewidth',1);plot(M0_test,a1_1,'--k','Linewidth',1);
plot(M0_test,a1,'.-k','Linewidth',1);plot(M0_test,a10,'k','Linewidth',1)
xlabel('log M0_L_P_D_T');ylabel('log(Source Radius, m)');box on
saveas(gca,fullfile(path_Reg,'App2_Mode.png'));
close
end
......
...........................................................................
[Depth2,Pl,logM0c,logM0t_C,d_logM0t_C,DeltaS_C,dDeltaS_C,loga_C,dloga_C,TX_C,dTX_C]=App2('c',Mweve,PL.MainAmp(:,3),...
PL.MainTime(:,3)-PL.MainTime(:,1),Depth,coefE(9,:),coefE(10,:));
App2C.event = EveL;
App2C.logM0_Cata=round(100.*logM0c)./100;
App2C.Mw_Cata=(App2C.logM0_Cata-9.1)./1.5;
App2C.logM0_LPDT=round(100.*logM0t_C)./100;
App2C.d_logM0_LPDT=round(1000.*d_logM0t_C)./1000;
App2C.Mw_LPDT=(App2C.logM0_LPDT-9.1)./1.5;
App2C.d_Mw_LPDT=(App2C.d_logM0_LPDT)./1.5;
App2C.StressDrop_MPa=round(10000.*DeltaS_C)./10000;
App2C.d_StressDrop=round(1000.*dDeltaS_C)./1000;
App2C.Radius_Km=round(10.*((10.^loga_C)/1000))./10;
App2C.d_Radius=round(1000.*((10.^dloga_C)/1000))./1000;
App2C.Duration_s=round(1000.*TX_C)./1000;
App2C.d_Duration=round(1000.*dTX_C)./1000;

clear T2; T2 = struct2table(App2C);
writetable(T2,output,'WriteRowNames',false,'Sheet','App2_Carvature');
...........................................................................
    ...........................................................................
if (strcmp(QFilter,'on')==1)
    AppQ.Curvature.event = EveL;
    AppQ.Curvature.logM0_Cata=round(100.*logM0c)./100;
    AppQ.Curvature.Mw_Cata=(App2M.logM0_Cata-9.1)./1.5;
    AppQ.Curvature.logM0_LPDT=round(100.*(1.5*QCorI_C(:,1)+9.1))./100;
    AppQ.Curvature.d_logM0_LPDT=round(1000.*d_logM0t_M)./1000;
    AppQ.Curvature.Mw_LPDT=round(100.*QCorI_C(:,1))./100;
    AppQ.Curvature.d_Mw_LPDT=(AppQ.Curvature.d_logM0_LPDT)./1.5;
    AppQ.Curvature.StressDrop_MPa=round(10000.*QCorI_C(:,2))./10000;
    AppQ.Curvature.d_StressDrop=round(10000.*dDeltaS_M)./10000;
    AppQ.Curvature.Radius_Km=round(10.*((10.^QCorI_C(:,6))/1000))./10;
    AppQ.Curvature.d_Radius=round(1000.*((10.^dloga_M)/1000))./1000;
    AppQ.Curvature.Duration_s=round(100.*2.*QCorI_C(:,3))./100;
    AppQ.Curvature.d_Duration=round(1000.*dTX_M)./1000;        
    T3 = struct2table(AppQ.Curvature);
%     writetable(T3,fullfile(path_Val,'Qfilter.txt'),'Delimiter','\t','WriteRowNames',false);
    writetable(T3,output,'WriteRowNames',false,'Sheet','App2_Qfilter_Curvature');
    ...................................................................PLOT
    if (size(coefE,2)~=1)
    subplot 221; 
    hold on
    scatter(logM0c,AppQ.Curvature.logM0_LPDT,30,[0.7 0.7 0.7],'fill','MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
    subplot 222; 
    hold on
    scatter(AppQ.Curvature.logM0_LPDT,log10(AppQ.Curvature.Duration_s),30,[0.7 0.7 0.7],'fill','MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
    subplot 223; 
    hold on;
    scatter(AppQ.Curvature.logM0_LPDT,log10(AppQ.Curvature.Radius_Km.*1000),30,[0.7 0.7 0.7],'fill','MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
    subplot 224; 
    hold on
    h=histfit(log10(AppQ.Curvature.StressDrop_MPa));   
    h(1).FaceColor = [0.7 0.7 0.7]; 
    h(2).Color = [.6 .6 .6];
    h(1).FaceAlpha=0.7;
    fsd=find((h(2).YData>= max(h(2).YData)));
    maxnorm =10.^(h(2).XData(fsd(1))); 
    xline(h(2).XData(fsd(1)),'Linewidth',2,'color',[.6 .6 .6])
%     xt = xticklabels; 
%     xnum=10.^( cellfun(@str2num,xt) );
%     xticklabels(round(1000.*xnum)./1000)
    title(['Stress Drop= ',num2str(round(1000.*maxnorm)./1000),' MPa'],'FontSize',9)
    end
end
...........................................................................
if (size(coefE,2)~=1)
clear Mo_test M0_test_min M0_test_max; 
if (strcmp(QFilter,'on')==1)
M0_test_min=min([AppQ.Curvature.logM0_LPDT App2C.logM0_LPDT],[],'all');
M0_test_max=max([AppQ.Curvature.logM0_LPDT App2C.logM0_LPDT],[],'all');
else
M0_test_min=min(App2C.logM0_LPDT);
M0_test_max=max(App2C.logM0_LPDT);
end

M0_test=M0_test_min:((M0_test_max-M0_test_min)/10):M0_test_max;
a1_01=log10((7.*(10.^M0_test)./(16*0.01*1e6)).^(1./3));
a1_1=log10((7.*(10.^M0_test)./(16*0.1*1e6)).^(1./3));
a1=log10((7.*(10.^M0_test)./(16*1*1e6)).^(1./3));
a10=log10((7.*(10.^M0_test)./(16*10*1e6)).^(1./3));
subplot 221;hold on;box on;grid on
subplot 222;hold on;box on;grid on
subplot 224;hold on;box on;grid on
subplot 223;hold on;box on;grid on
plot(M0_test,a1_01,':k','Linewidth',1);plot(M0_test,a1_1,'--k','Linewidth',1);
plot(M0_test,a1,'.-k','Linewidth',1);plot(M0_test,a10,'k','Linewidth',1)
xlabel('log M0_L_P_D_T');ylabel('log(Source Radius, m)');box on
saveas(gca,fullfile(path_Reg,'App2_Carvature.png'));
close
end
...........................................................................
if (strcmp(QFilter,'on')==1)
%     QPLS=AppQ.Curvature.logM0_LPDT-log10(QCorI_C(:,3)./Ap); 
    QPLS=log10(QCorI_C(:,4));  
    [logM0App1]=App1Cal(Mweve,QPLS,SD);
else
    [logM0App1]=App1Cal(Mweve,PL.MainAmp(:,4),SD);
end
...........................................................................
d_logM0t_App1=d_logM0t_M;
App1S.event = EveL;
App1S.logM0_Cata=round(100.*logM0c)./100;
App1S.Mw_Cata=round(100*(App1S.logM0_Cata-9.1)./1.5)./100;
App1S.logM0_LPDT=round(100.*logM0App1)./100;
App1S.d_logM0_LPDT=round(1000.*d_logM0t_App1)./1000;
App1S.Mw_LPDT=round(100*(App1S.logM0_LPDT-9.1)./1.5)./100;
App1S.d_Mw_LPDT=(App1S.d_logM0_LPDT)./1.5;
T1 = struct2table(App1S);
% writetable(T1,fullfile(path_Val,'App1.txt'),'Delimiter','\t','WriteRowNames',false);
writetable(T1,output,'WriteRowNames',false,'Sheet','App1');
...........................................................................
...........................................................................
[p] = App1Val(Mweve,PL.MainAmp(:,4), SD1, SD2, SDN);
App1ValS.alpha = p(:,1);
App1ValS.beta = p(:,2);
App1ValS.stress_drop = p(:,3);
App1ValS.RMS = p(:,4);
T4 = struct2table(App1ValS);
writetable(T4,output,'WriteRowNames',false,'Sheet','App1_FindMinStress');
clear fsd; fsd=find(T4.RMS==min(T4.RMS));
minSD=T4.stress_drop(fsd);
figure(1); set(gca,'FontSize',20)
semilogx(T4.stress_drop,T4.RMS,'Linewidth', 2)
xlabel('Stress Drop (MPa)'); ylabel('RMS')
hold on; scatter(minSD,T4.RMS(fsd),80,'b','fill');   
grid on
title(['Min Stress Drop: ', num2str(round(100*minSD)/100),' MPa']);
saveas(figure (1),fullfile(path_Reg,'App1-MinSD.png'));
close
...........................................................................
cd (path_out);ginfo=dir('*xls');
result = []; 

for i = 1:size(ginfo,1)
    [~,~, filecontent]= xlsread(fullfile([path_out,'\Values'], ginfo(i).name)); 
    xlswrite(output,filecontent,ginfo(i).name,i)
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
clc; close all; cd(path0)
save([path_out,'\',Region,'.mat'])

