function [Depth2,Pl,logM0c,logM0t,d_logM0t,DeltaS,d_DeltaS,loga,d_loga...
    ,TX,dTX]=App2(p,Mw,Pl,TC,Depth,dp,dTC)
global Vp Ap Vr
if size(Mw,1)==1;Mag=Mw';else; Mag=Mw; end
if size(Pl,1)==1;Pl=Pl';end    
if size(TC,1)==1;TC=TC';end 
if size(Depth,1)==1;Depth2=Depth';else; Depth2=Depth; end
if size(dTC,1)==1;dTC=dTC';end 
if size(dp,1)==1;dp=dp';end 

% [Mag , I]=sort(Mag);
% Pl=Pl(I);TC=TC(I);Depth2=Depth2(I);
% dp=dp(I);dTC=dTC(I);
dTX=2.*dTC;

% clear fmT; fmT=find(dTC>=(0.1.*TC));
% dTX(fmT,:)=0.1;dTC(fmT,:)=0.05;

% Pl(fmT,:)=[];TC(fmT,:)=[];Depth2(fmT,:)=[];
% dp(fmT,:)=[];
% Mag(fmT,:)=[];

logM0c=(1.5*Mag)+9.1; %logM0
Vp=Vp.*ones(length(Pl),1); 
Vr=Vr.*ones(length(Pl),1);
Ap=Ap.*ones(length(Pl),1);
TX=2.*TC; 
% =========================================log M0
logM0t=(Pl-log10(Ap./TC));
d_logM0t=dp-((dTC.*TC)./(ones(length(dp),1).*log(10)));
% d_M0t=(10.^logM0t).*d_logM0t.*(ones(length(dp),1).*log(10));
% ========================================
% ========================================log a
loga=log10(TC./((1./Vr)-(2./(pi().*Vp))));
d_loga=((dTC.*TC)./(ones(length(dp),1).*log(10)));
% d_a=(10.^loga).*d_loga.*(ones(length(dp),1).*log(10));
M0t=10.^logM0t; 
at=(10.^loga);
% =========================================
DeltaS=(7.*M0t)./(16.*(at.^3)); %Pa
DeltaS=DeltaS./(1e6);
% d_DeltaS=( (7.*d_M0t.*(16.*(at.^3)))-(16.*3.*(at.^2).*d_a.*7.*M0t) )./((16.*(at.^3)).^2);
d_DeltaS=d_logM0t-3.*d_loga;
% d_DeltaS=d_DeltaS./(1e6);
% ============================================
if (size(Mw,1)~=1)
%     figure (3);
if(strcmp(p,'m')==1)
    subplot 221; 
    hold on
    errorbar(logM0c,logM0t,d_logM0t,'kd');% 'LineStyle','none');
    scatter(logM0c,logM0t,60,'c','fill','MarkerEdgeColor','k');
    plot(logM0c,logM0c,'k','Linewidth',2)
    xlabel('log M0_c_a_t_a');ylabel('log M0_L_P_D_T');box on
    .....................................
    subplot 222; 
    hold on
    errorbar(logM0t,log10(TC),((dTC.*TC)./(ones(length(dp),1).*log(10))),'kd');% 'LineStyle','none');
    scatter(logM0t,log10(TC),60,'c','fill','MarkerEdgeColor','k');
    % set(gca, 'YScale', 'log')
    xlabel('log M0_L_P_D_T');ylabel('log(HD,s)');box on    
    subplot 223; 
    hold on;
    errorbar(logM0t,loga, d_loga,'kd');% 'LineStyle','none');
    scatter(logM0t,loga,60,'c','fill','MarkerEdgeColor','k');

    subplot 224; 
    hold on
    h=histfit(log10(real(DeltaS)));   
    h(1).FaceColor = 'c'; 
    h(2).Color = [.2 .2 .2];
    fsd=find((h(2).YData>= max(h(2).YData)));
    maxnorm =10.^(h(2).XData(fsd(1)));
    xline(h(2).XData(fsd(1)),'Linewidth',2,'color',[.2 .2 .2])
    xt = xticklabels; 
    xnum=10.^( cellfun(@str2num,xt) );
    xticklabels(round(1000.*xnum)./1000)
    title(['Stress Drop= ',num2str(round(1000.*maxnorm)./1000),' MPa'],'FontSize',9) 
%     errorbar(logM0t,(DeltaS),10.^(d_DeltaS),'kd');% 'LineStyle','none');
%     scatter(logM0t,(DeltaS),60,'c','fill','MarkerEdgeColor','k');
%     set(gca, 'YScale', 'log')
%     xlabel('log M0_L_P_D_T');
    xlabel('Stress Drop (MPa)');box on
%     sgtitle('(LPDT)');
elseif(strcmp(p,'c')==1)
    subplot 221; 
    hold on
    errorbar(logM0c,logM0t,d_logM0t,'kd');% 'LineStyle','none');
    scatter(logM0c,logM0t,60,'dm','fill','MarkerEdgeColor','k');
    plot(logM0c,logM0c,'k','Linewidth',2)
    xlabel('log M0_c_a_t_a');ylabel('log M0_L_P_D_T');box on
    .....................................
    subplot 222; 
    hold on
    errorbar(logM0t,log10(TC),((dTC.*TC)./(ones(length(dp),1).*log(10))),'kd');% 'LineStyle','none');
    scatter(logM0t,log10(TC),60,'dm','fill','MarkerEdgeColor','k');
    % set(gca, 'YScale', 'log')
    xlabel('log M0_L_P_D_T');ylabel('log(HD,s)');box on    
    subplot 223; 
    hold on;
    errorbar(logM0t,loga, d_loga,'kd');% 'LineStyle','none');
    scatter(logM0t,loga,60,'dm','fill','MarkerEdgeColor','k');

    subplot 224; 
    hold on
    h=histfit(log10(real(DeltaS)));   
    h(1).FaceColor = 'm'; 
    h(2).Color = [.2 .2 .2];
    fsd=find((h(2).YData>= max(h(2).YData)));
    maxnorm =10.^(h(2).XData(fsd(1)));
    xline(h(2).XData(fsd(1)),'Linewidth',2,'color',[.2 .2 .2])
    xt = xticklabels; 
    xnum=10.^( cellfun(@str2num,xt) );
    xticklabels(round(1000.*xnum)./1000)
    title(['Stress Drop= ',num2str(round(1000.*maxnorm)./1000),' MPa'],'FontSize',9) 
%     errorbar(logM0t,(DeltaS),10.^(d_DeltaS),'kd');% 'LineStyle','none');
%     scatter(logM0t,(DeltaS),60,'dm','fill','MarkerEdgeColor','k');
%     set(gca, 'YScale', 'log')
%     xlabel('log M0_L_P_D_T');
    xlabel('Stress Drop (MPa)');box on
else
    subplot 221; 
    hold on
    errorbar(logM0c,logM0t,d_logM0t,'kd');% 'LineStyle','none');
    scatter(logM0c,logM0t,70,'gd','fill','MarkerEdgeColor','k');
    plot(logM0c,logM0c,'k','Linewidth',2)
    xlabel('log M0_c_a_t_a');ylabel('log M0_L_P_D_T');box on
    .....................................
    subplot 222; 
    hold on
    errorbar(logM0t,log10(TC),((dTC.*TC)./(ones(length(dp),1).*log(10))),'kd');% 'LineStyle','none');
    scatter(logM0t,log10(TC),70,'gd','fill','MarkerEdgeColor','k');
    % set(gca, 'YScale', 'log')
    xlabel('log M0_L_P_D_T');ylabel('log(HD,s)');box on    
    subplot 223; 
    hold on;
    errorbar(logM0t,loga, d_loga,'kd');% 'LineStyle','none');
    scatter(logM0t,loga,70,'gd','fill','MarkerEdgeColor','k');
    subplot 224; 
    hold on
    h=histfit(log10(DeltaS));   
    h(1).FaceColor = 'g'; 
    h(2).Color = [.2 .2 .2];
    fsd=find((h(2).YData>= max(h(2).YData)));
    maxnorm =10.^(h(2).XData(fsd(1)));
    xline(h(2).XData(fsd(1)),'Linewidth',2,'color',[.2 .2 .2])
    xt = xticklabels; 
    xnum=10.^( cellfun(@str2num,xt) );
    xticklabels(round(1000.*xnum)./1000)
    title(['Stress Drop= ',num2str(round(1000.*maxnorm)./1000),' MPa'],'FontSize',9) 
%     errorbar(logM0t,(DeltaS),10.^(d_DeltaS),'kd');% 'LineStyle','none');
%     scatter(logM0t,(DeltaS),70,'gd','fill','MarkerEdgeColor','k');
%     set(gca, 'YScale', 'log')
    xlabel('log M0_L_P_D_T');ylabel('Stress Drop (MPa)');box on
end
end

end