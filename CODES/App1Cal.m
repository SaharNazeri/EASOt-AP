function [logM0]=App1Cal(Mw,LPXTin,SD)
global path_Reg
global Vr Vp Ap
if size(Mw,1)==1;Mag=Mw';else; Mag=Mw; end
if size(LPXTin,1)==1;LPXT=LPXTin';else; LPXT=LPXTin; end  


StressDrop=SD*1e6;  %Pa
K=((7/16)^(1/3)) .* ((1./StressDrop).^(1/3)) .* ((1./Vr)-(2./(pi().*Vp)));
logM0=(3/2).*(LPXT-log10(Ap./K));

dobs=1.5*Mag+9.1;
% MagND=(logM0-9.1)./1.5; 
MagND=logM0; 
if (size(Mag,1)~=1)
    figure; %set(gcf,'units','centimeter','position',[8 4 15 20]); 
    hold on;box on
%     set(gca,'fontsize', 15); 
    [p,S] = polyfit(dobs,MagND,1);
    [y_fit,delta] = polyval(p,dobs,S);
%     caption1=['M (LPDT) = ' sprintf('%.2f',p(1)) '*M_w + ' sprintf('%.2f',p(2))];
    caption1=['log M0 (LPDT) = ' sprintf('%.2f',p(1)) '*log M0 + ' sprintf('%.2f',p(2))];
    caption2=['SE  = ' sprintf('%.2f',mean(2*delta))];
    caption3=['SD (MPa)= ' sprintf('%.2f',StressDrop/1e6)];

    scatter(dobs,MagND,100,'b','fill','MarkerEdgeColor','k'); 
    plot(dobs,y_fit,'b','LineWidth',3)
    [dobsP,I]=sort(dobs); 
    clear P1; clear P2;
%     P1=y_fit+(2*delta); P1=P1(I);
%     P2=y_fit-(2*delta); P2=P2(I);
    P1=y_fit+std(dobsP); P1=P1(I);
    P2=y_fit-std(dobsP); P2=P2(I);
    plot(dobsP,P1,':b',dobsP,P2,':b','LineWidth',2.5)
    % plot(dobs,P1,':b',dobs,P2,':b','LineWidth',2.5)
    hold on;
    plot(dobs,dobs,'k','LineWidth',3); %legend('Estimate','True')
%     xlabel('Mw (cata)','Fontsize',15);ylabel(['M (LPDT)'],'Fontsize',15);
    xlabel('log M0_c_a_t_a');ylabel('log M0_L_P_D_T');
    ax=gca;
%     ax.XAxis.FontSize = 15;
%     ax.YAxis.FontSize = 15;
    grid on
%     title ({[caption1 '    ;  ' caption2],[caption3]});
      title (caption1);
    ax1 = gca; % current axes
    ax1.XColor = 'k';
    ax1.YColor = 'k';
%     lim1=(min([dobs MagND],[],'all'));
%     lim2=(max([dobs MagND],[],'all'));
%     xlim([lim1 lim2])
%     ylim([lim1 lim2])
    saveas(gca,fullfile(path_Reg,'App1.png'));

    close all
end

end