cd (path_sac); g=dir('*');
disp('LPDT Package: Charactrization of the Source Time Function in time domain')
disp('Sahar Nazeri, Naples, Italy')
disp('Reference, Zollo et al., 2021')
disp('************************************************************************')
disp(['Region: ',Region,',  Total Number of events: ', num2str(length(g)-2)])

nuM=1; % it is counting the used events
nuMSta=1;
Depth=[];Time_PL=[]; Amp_PL=[];Events=[];
for l=3:length(g)   
    cd (g(l).name); 
    s=dir(convertCharsToStrings(dirparam));
    if(isempty(s)==1);disp(['** Warning: folder of "',g(l).name,'" is empty']);
    cd ..;continue; end
    a=s(1).name; tr=readsacs(a);
    ....................................................Event and magnitude
    Magco=tr.h1(40);
%     if (Magco==-12345) || (isnan(Magco)==1) || (isinf(Magco)==1)
%     clear fM; fM=find(strcmp(g(l).name,eve)==1); Magco=M(fM);
%     end    
    disp('======================================================')
    disp(['#  ',num2str(l-2),'   Event ',num2str(l-2),': ',g(l).name,'    , Mw: ', num2str(Magco)])
    fprintf('\n');
    ........................... using the stations whitin the range of Rmax
    if (Magco>=5); Rmax=RmM; elseif(Magco>=3 && Magco<5); Rmax=Rm5; 
    elseif(Magco>=1 && Magco<3); Rmax=Rm3; else; Rmax=Rm1; end
       
    if (strcmp(WinFix,'Yes')==1)
        win_stopMax=WinMax;
    else
        win_stopMax=S_P_Coef*Rmax;
    end

    LPDT=[];LPDT_mean=[];Dis=[];DisR=[];Travel=[];LPDT_Max=[];%Pl_error=[];
    num=1; % it is counting the traces   
    
    figure (1); set(gcf,'WindowState','maximized');
    set(gca,'FontSize',20)
    sgtitle([num2str(g(l).name),'     Mw =',num2str(mean(Magco))]);
   
    subplot 222; hold on; box on; grid on;
    set(gca,'Linewidth',1.5,'FontSize',12)
    ylabel('Displacement, (m)','FontSize',16); 
    xlabel('Time (s)','FontSize',16);   
    xline(0,'--');
    subplot 224; hold on; box on; grid on;
    set(gca,'Linewidth',1.5,'FontSize',12)
    ylabel('LPDT (SI unit)','FontSize',16);   
    xlabel('Time (s)','FontSize',16);   
    xline(0,'--');
    subplot (2,2,[1,3]);
    geoscatter(0,0,50,'fill','p','LineWidth',0.05,...
    'MarkerEdgeColor','k');geobasemap topographic
   ........................................................................ 
    for i=1:length(s) 
            a=s(i).name; tr=readsacs(a);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Excluding Non-Picked waveforms
            if (tr.h1(9)==-12345); continue; end
            Ppick=tr.h1(9);
            if(Ppick<tr.bt);Ppick=tr.bt;end
            T0=Ppick-tr.bt;
            T0=round(T0*100)/100;          
            tau=tr.tau; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Excluding waveforms outside the given distance
            if (isnan(tr.salt) == 0);depCo=tr.edep+(tr.salt/1000);
            else depCo=tr.edep;end
            Rhypo=sqrt(tr.distkm^2+depCo^2);
            if (Rhypo>Rmax); continue; end           
            win_stop=S_P_Coef*Rhypo;            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %read the waveforms
            if(strcmp(Unit,'Yes'))==1 
                accf=(tr.trace).*PU; %m/s
            else
%                 accf=(tr.trace).*(tr.h1(4)); %m/s
                accf=(tr.trace).*(tr.par2/tr.par3); %m/s
            end
            ...............................................................           
            time=0:tau:(length(accf)-1)*tau;
            time=time-T0;
            accf=detrend(accf,'constant');
            accf=detrend(accf,'linear');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % acc, vel, dis 
            if (time(1,1)<-5)    
                acc0=accf(time>=-5 & time<=(win_stop+9)); 
                time0=time(time>=-5 & time<=(win_stop+9));                
            elseif (time(1,1)>=-5 && time(1,1)<-2)    
                acc0=accf(time>=-2 & time<=(win_stop+9));
                time0=time(time>=-2 & time<=(win_stop+9));
            else 
                fprintf(2,['** Warning: Data Excluded:',a,' less that 2 seconds before P-wave\n']);
                continue;
            end
            acc0 = acc0.*tukeywin(length(acc0),0.25);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Noise & SNR
            SigN=acc0(time0>=(-1.5) & time0<=(-0.5));
            SigS=acc0(time0>=(0) & time0<=(win_stop/2));
            SNRSig=20.*log10(max(abs(SigS))./max(abs(SigN)));
            if (SNRSig<=SNR)
                fprintf(2,['** Warning: Station of "',tr.staname,'" is excluded, because of low SNR\n']);
                continue;
            end                
%             ...................................
            if (strcmp(wave_type,'A'))==1 
                vel=real(cumtrapz(acc0).*tau); 
                vel=detrend(vel,'constant');
                vel=detrend(vel,'linear'); 
            else
                vel=acc0;
            end          
            
            [bhV,ahV] = butter(nPolV,[cor1V, cor2V]/(1/(2*tau)),'bandpass');
            if(strcmp(Vfilter,'Yes'))==1 
                vel=filter(bhV,ahV,vel);  
            end
%             [bh,ah] = butter(nPol,[cor1, cor2]/(1/(2*tau)),'bandpass');
            [bh,ah] = butter(nPol,cor1/(1/(2*tau)),'high');
            diss=real(cumtrapz(vel).*tau); 
            diss=detrend(diss,'constant');
            diss=detrend(diss,'linear');                        
            ...................................
            dissh=filter(bh,ah,diss);   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Noise correction
            accN=acc0(time0>=(-0.2) & time0<=(-0.1));
            velN=vel(time0>=(-0.2) & time0<=(-0.1));
            disshN=dissh(time0>=(-0.2) & time0<=(-0.1));
            acc0=acc0-mean(accN);vel=vel-mean(velN);dissh=dissh-mean(disshN);
            %....... For Plot
            if (time(1,1)<-2)  
                disshp=dissh(time0>=-0.5 & time0<=win_stop);
                timep=time0(time0>=-0.5 & time0<=win_stop);
            else
                disshp=dissh(time0>=0 & time0<=win_stop);
                timep=time0(time0>=0 & time0<=win_stop);  
            end
            %.......
            acc0=acc0(time0>=0 & time0<=win_stop);
            vel=vel(time0>=0 & time0<=win_stop);
            dissh=dissh(time0>=0 & time0<=win_stop);
            
            timeN=time(time>=0 & time<=win_stop);
            timeMax=time(time>=0 & time<=win_stopMax);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % LPDT
            Cd=1;
            dissh_abs=abs(dissh);
            clear LPDT_t; LPDT_t=log10(dissh_abs)+Cd.*log10(1000*Rhypo);
            if(win_stopMax>=win_stop)
              LPDT_t( (length(timeN)+1):length(timeMax))=1e-20;
            else
              LPDT_t=LPDT_t(1:length(timeMax));
            end
            LPDT_t(isinf(LPDT_t)==1)=[];
            if(isempty(LPDT)==1)
                LPDT(:,num)=LPDT_t;
            else
                if(size(LPDT,1)<=size(LPDT_t,1))
                    LPDT(:,num)=LPDT_t(1:size(LPDT,1),1);
                else
                    LPDT_t(size(LPDT_t,1):size(LPDT,1),1)=LPDT_t(end,1);
                    LPDT(:,num)=LPDT_t;
                end
            end

            Dis(1,num)=log10(1000*Rhypo);
            DisR(1,num)=Rhypo;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
            if(num==1)
            subplot (2,2,[1,3]);cla
            geoscatter(tr.elat,tr.elon,400,'fill','rp','LineWidth',0.05,...
            'MarkerEdgeColor','k');                 
            end
            subplot (2,2,[1,3]);hold on
            geoscatter(tr.slat,tr.slon,250,'fill','v','LineWidth',0.05,...
                'MarkerFaceColor','b','MarkerEdgeColor','k');
            text(tr.slat-0.0022,tr.slon-0.0003, tr.staname,'FontSize',9)
            subplot 222; 
            hold on;
            plot(timep,disshp,'Linewidth',1.2,'color',[0.7 0.7 0.7]);      
            text(timep(end),disshp(end), tr.staname,'FontSize',9)
            xline(0,'--');
            subplot 224; 
            hold on; 
            plot(timeN,log10(dissh_abs)+Cd.*log10(Rhypo*1000),'Linewidth',1.2,'color',[0.7 0.7 0.7]);            
            xline(0,'--');   
            ...............................................................
            num=num+1;   
            nuMSta=nuMSta+1;
%             pause
    end
    subplot (2,2,[1,3]); 
    geoscatter(tr.elat,tr.elon,400,'fill','rp','LineWidth',0.05,...
    'MarkerEdgeColor','k');                 
    geobasemap topographic
    subplot 222; xlim([-0.2 inf]);
    if(length(LPDT)==1) 
         cd .. ; 
         disp (['LPDT is NOT generated for this event', g(l).name]);
         disp (['Less than' ,minSta ,'P-wave picks, Or less than 2 seconds pre P-wave window']);
         continue; 
     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mean LPDT
    num_c=1;nuMStaLPDT=[];   
    for i=1:size(LPDT,1)
        clear LPDT_t f0; LPDT_t=LPDT(i,LPDT(i,:)~=1e-20);        
        if (length(LPDT_t)>=minSta)
            LPDT_mean(num_c,1)=median(LPDT_t);
            LPDT_mean(num_c,2)=min([(median(LPDT_t)-min(LPDT_t));...
                (max(LPDT_t)-median(LPDT_t));std(LPDT_t)]);
            Dis_mean(num_c,1)=mean(Dis(1,1:length(LPDT_t)));
            nuMStaLPDT(num_c,1)=length(LPDT_t);
            num_c=num_c+1;
        end        
    end
    if(length(LPDT_mean)<=1)
        close figure 1; cd .. ;
        disp('Number of stations is not enough');
        continue; 
    end 
    
    clear tY; tY=find(isnan(LPDT_mean(:,1))==1 | isinf(LPDT_mean(:,1))==1);
    if isempty(tY)~=1;LPDT_mean(tY,:)=[]; Dis_mean(tY)=[]; end
    clear tY; tY=find(isnan(LPDT_mean(:,2))==1 | isinf(LPDT_mean(:,2))==1);
    if isempty(tY)~=1;LPDT_mean(tY,:)=[]; Dis_mean(tY)=[]; end
    timeN=timeMax(1,1:length(LPDT_mean));
    
    % Maximum of Mean LPDT-New method
    clear YA; YA=FMax(LPDT_mean(:,1));
    subplot 224; hold on; 
    plot(timeN,LPDT_mean(:,1),'k','Linewidth',2.5)
    plot(timeN,YA,'r','Linewidth',3)
    xlim([0 inf]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pause(2)
    saveas(figure (1),fullfile(path_Fig,[num2str(g(l).name)  '-1.png']));
    close figure 1 
    Xn=timeN(timeN>=0 ); Yn=YA(timeN>=0); 
    clear tY; tY=find(isnan(Yn)==1 | isinf(Yn)==1); 
    Yn(tY)=[]; Xn(tY)=[];
    .......................................................................
    %======================================================================
    if (strcmp(WinFix_fit,'Yes'))==1
        clear fFW; fFW=find(Xn>WinMax_fit);
        Xn(:,fFW)=[]; Yn(:,fFW)=[];
        timeN(:,fFW)=[]; LPDT_mean(fFW,:)=[]; nuMStaLPDT(fFW,:)=[];
    end
    figure (2);
    title([num2str(g(l).name),'     Mw =',num2str(mean(Magco))]);
    box on;grid on;
    hold on; plot(timeN,LPDT_mean(:,1),'color',[0.7 0.7 0.7],'Linewidth',1.5)   
   
    
    [pks,locs,Yn_temp]=SMLPXT(Magco,timeN,LPDT_mean(:,1));   
    if(isempty(pks)==1);pks(1,1)=LPDT_mean(end,1);locs(1,1)=timeN(end);end
    coefE(7,nuM)=pks(1,1);coefE(8,nuM)=locs(1,1); 
    clear YA; YA=FMax(Yn_temp);
    plot(Xn,YA,'r','Linewidth',2);
    .......................................................................
    .......................................................................
    [Plg,ygrid]=fit(Xn,YA);  
    [PlgHB,ygridHB]=fitHB(Xn,YA);    
    chekFit{nuM,1}=FitFunc; chekFit{nuM,2}=FitFunc; chekFit{nuM,3}=FitFunc;
    if (strcmp(chekFit{nuM,1},'HB'))==1
        ygrid=ygridHB; Plg=PlgHB;
    end
    coefE(1:3,nuM)=Plg;
    if (strcmp(chekFit{nuM,1},'HB'))==1
        coefE(2,nuM)=coefE(2,nuM)-Xn(1);
    else
        coefE([2 3],nuM)=coefE([2 3],nuM)-Xn(1);
    end
    coefE(4,nuM)=Yn(1);
    
    plot(Xn,ygrid,'k','LineWidth',1.2); 
    text(Xn(end),ygrid(end), chekFit{nuM,1},'FontSize',9)
    ygridD=diff((ygrid))./tau;
    tD0=find(ygridD==max(ygridD));
    yD=ygrid(tD0(1)+1);
    tD=tau.*(tD0(1)+1);
    tPl=find(YA==mode(YA));
    yPl=YA(tPl(1));
    tPl=tau.*tPl(1);
    [CNum, PLbin]=hist(YA,10);
    if(tPl~=tau)
        clear fM fm PL_level PL_Time; fM=find(CNum==max(CNum));         
        PL_level=PLbin(fM(1));
        fm=Xn(YA>=PL_level);
        PL_Time=fm(1);
        if ((tPl-PL_Time)>tau)
            tPl=PL_Time; yPl=PL_level;
        end
    end
    .......................................................................
    .......................................................................
    if (strcmp(chekFit{nuM,1},'exp'))==1
        clear fm K2L K2 
        [K2,K2L]=demo2d(Xn,ygrid);
        fm=Xn(K2L==max(K2L)); fm_min=Xn(K2L<=0.2*max(K2L));
        if (Magco>=5.5)
            fm=Xn(K2L==max(K2L)); fm_min=Xn(K2L<=0.1*max(K2L));
        end
        fm_nextMax=find(fm_min>fm);
        if (isempty(fm_nextMax)==1)
            disp(['** Warning: LPDT is not well generated for event of "',g(l).name,'"']);
            close; cd ..
            continue
        end
        fm_sol=find(Xn>=fm_min(fm_nextMax(1)));
        coefE(6,nuM)=Xn(fm_sol(1));
        coefE(5,nuM)=ygrid(Xn==coefE(6,nuM)); 
        if (Magco<5)
            clear ftemp; ftemp=find(Xn>=fm(1));
            coefE(6,nuM)=Xn(ftemp(1));    %max curvature
            coefE(5,nuM)=ygrid(Xn==coefE(6,nuM)); 
        end
    else
        if(tau*(tD0(1)+7)<Xn(1,end))
            Xnc=Xn(Xn>=tau*(tD0(1)+7));ygridc=ygrid(Xn>=tau*(tD0(1)+7));
        else
            Xnc=Xn;ygridc=ygrid;
        end
        clear fm K2L K2 
        [K2,K2L]=demo2d(Xnc,ygridc);
        fm=Xnc(K2L==max(K2L)); fm_min=Xnc(K2L<=0.2*max(K2L));
        if (Magco>=5.5)
            fm=Xnc(K2L==max(K2L)); fm_min=Xnc(K2L<=0.1*max(K2L));
        end
        fm_nextMax=find(fm_min>fm(1)); fm_sol=find(Xnc>=fm_min(fm_nextMax(1)));
        coefE(5,nuM)=max(ygrid); 
        coefE(6,nuM)=Xnc(fm_sol(1));    
        if (Magco<5)
            clear ftemp; ftemp=find(Xnc>=fm(1));
            coefE(6,nuM)=Xnc(ftemp(1));    %max curvature
            coefE(5,nuM)=ygrid(Xn==coefE(6,nuM));
        end
        clear fm ym K2L K2 
        if((tD0(1)+7)~=1)            
            [K2,K2L]=demo2d(Xn(1,1:(tD0(1)+7)),smooth(ygrid(1:(tD0(1)+7),1)));
            fm=Xn(K2L==max(K2L))+2.*tau;
            ym=ygrid(Xn==fm(end));            
        else
%             tD0(1)=tD0(1)-1; 
            fm=tau*(tD0(1));
            ym=ygrid(1);
        end

    end
    PL.MainTime(nuM,:)=[fm(end) tD coefE(6,nuM) tPl];
    PL.MainAmp(nuM,:)=[ym yD coefE(5,nuM) yPl];
    scatter(PL.MainTime(nuM,4),PL.MainAmp(nuM,4),60,'c','fill','MarkerEdgeColor','k')
    scatter(PL.MainTime(nuM,3),PL.MainAmp(nuM,3),60,'dm','fill','MarkerEdgeColor','k')
    plot([PL.MainTime(nuM,1) PL.MainTime(nuM,1)],[PL.MainAmp(nuM,1) PL.MainAmp(nuM,4)],':k','linewidth',1.5)
%%
    %======================================================================
    LPDT_M1=LPDT_mean(:,1)-LPDT_mean(:,2);
    LPDT_M2=LPDT_mean(:,1)+LPDT_mean(:,2);
    plot(timeN,LPDT_M1,'color',[0.7 0.7 0.7],'Linewidth',1.5)
    plot(timeN,LPDT_M2,'color',[0.7 0.7 0.7],'Linewidth',1.5)
    clear YAM1 YAM2 Yn_tempM1 Yn_tempM2;
    [~,~,Yn_tempM1]=SMLPXT(Magco,timeN,LPDT_M1);
    [~,~,Yn_tempM2]=SMLPXT(Magco,timeN,LPDT_M2);
    YAM1=FMax(Yn_tempM1);YAM2=FMax(Yn_tempM2);
    plot(Xn,YAM1,'r','Linewidth',2); plot(Xn,YAM2,'r','Linewidth',2);
    ..................................................
    [PlgM1,ygridM1]=fit(Xn,YAM1);
    [PlgM1HB,ygridM1HB]=fitHB(Xn,YAM1);
   if (strcmp(chekFit{nuM,2},'HB'))==1
        ygridM1=ygridM1HB; PlgM1=PlgM1HB;
    end
    plot(Xn,ygridM1,'-k','LineWidth',1);
    text(Xn(end),ygridM1(end), chekFit{nuM,2},'FontSize',9)
        ygridM1D=diff((ygridM1))./tau;
    tD0M1=find(ygridM1D==max(ygridM1D));
    yDM1=ygridM1(tD0M1(1)+1);
    tDM1=tau.*(tD0M1(1)+1);
    tPlM1=find(YAM1==mode(YAM1));
    yPlM1=YAM1(tPlM1(1));
    tPlM1=tau.*tPlM1(1);
    if(tPlM1~=tau)
        clear CNum PLbin;[CNum, PLbin]=hist(YAM1,10);
        clear fM fm PL_level PL_Time; fM=find(CNum==max(CNum));         
        PL_level=PLbin(fM(1));
        fm=find(YAM1>=PL_level);
        PL_Time=fm(1);
        if ((tPlM1-PL_Time)>tau)
            tPlM1=PL_Time; yPlM1=PL_level;
        end
    end
    .......................................................................
    if (strcmp(chekFit{nuM,2},'exp'))==1    
        clear fm1 K2M1 K2LM1 TC1
        [K2M1,K2LM1]=demo2d(Xn,ygridM1);
        if(isnan(K2M1)==1)
            TC1=Xn(1);
        else
        fm1=Xn(K2LM1==max(K2LM1)); fm_min1=Xn(K2LM1<=0.2*max(K2LM1));
        if (Magco>=5.5)
            fm1=Xn(K2LM1==max(K2LM1)); fm_min1=Xn(K2LM1<=0.1*max(K2LM1));
        end
        fm_nextMax1=find(fm_min1>fm1(1)); fm_sol1=find(Xn>=fm_min1(fm_nextMax1(1)));
        TC1=Xn(fm_sol1(1));    
        if (Magco<5)
            clear ftemp; ftemp=find(Xn>=fm1(1));
            TC1=Xn(ftemp(1));    %max curvature
        end
    end
    else
            clear K2M1 K2LM1 Xnc ygridM1c
        if(tau*(tD0M1(1))<Xn(1,end))
            Xnc=Xn(Xn>=tau*(tD0M1(1)-1));ygridM1c=ygridM1(Xn>=tau*(tD0M1(1)-1));
        else
            Xnc=Xn;ygridM1c=ygridM1;
        end
                       
        [K2M1,K2LM1]=demo2d(Xnc,ygridM1c);
        if(isnan(K2M1)==1)
            TC1=Xnc(1);
        else
            fm2=Xnc(K2LM1==max(K2LM1)); fm_min2=Xnc(K2LM1<=0.2*max(K2LM1));
            if (Magco>=5.5)
                fm2=Xnc(K2LM1==max(K2LM1)); fm_min2=Xnc(K2LM1<=0.1*max(K2LM1));
            end
            fm_nextMax2=find(fm_min2>fm2(1)); fm_sol2=find(Xnc>=fm_min2(fm_nextMax2(1)));
            TC1=Xnc(fm_sol2(1));    
            if (Magco<5)
                clear ftemp; ftemp=find(Xnc>=fm2(1));
                TC1=Xnc(ftemp(1));    %max curvature
            end
        end
        clear fm ym K2L K2 
        if(tD0M1(1)~=1)            
            [K2,K2L]=demo2d(Xn(1,1:tD0M1(1)),smooth(ygridM1(1:tD0M1(1),1)));
            fm=Xn(K2L==max(K2L));
            ymM1=ygridM1(Xn==fm(end));
        else
%             tD0M1(1)=tD0M1(1)-1; 
            fm=tau*tD0M1(1);
            ymM1=ygridM1(1);
        end

        
    end
    clear ftemp; ftemp=find(timeMax>=TC1);
    PL.DownTime(nuM,:)=[fm(end) tDM1 TC1 tPlM1];
    PL.DownAmp(nuM,:)=[ymM1 yDM1 ygridM1(ftemp(1)) yPlM1];
    scatter(PL.DownTime(nuM,4),PL.DownAmp(nuM,4),60,'c','fill','MarkerEdgeColor','k')
    scatter(PL.DownTime(nuM,3),PL.DownAmp(nuM,3),60,'dm','fill','MarkerEdgeColor','k')
 %%
    %======================================================================
    [PlgM2,ygridM2]=fitHB(Xn,YAM2);
    [PlgM2HB,ygridM2HB]=fitHB(Xn,YAM2); 
    if (strcmp(chekFit{nuM,3},'HB'))==1
        ygridM2=ygridM2HB; PlgM2=PlgM2HB;
    end
    plot(Xn,ygridM2,'-k','LineWidth',1); 
    text(Xn(end),ygridM2(end), chekFit{nuM,3},'FontSize',9)
    ygridM2D=diff((ygridM2))./tau;
    tD0M2=find(ygridM2D==max(ygridM2D));
    yDM2=ygridM2(tD0M2(1)+1);
    tDM2=tau.*(tD0M2(1)+1);
    tPlM2=find(YAM2==mode(YAM2));
    yPlM2=YAM2(tPlM2(1));
    tPlM2=tau.*tPlM2(1);
    if(tPlM2~=tau)
        clear CNum PLbin;[CNum, PLbin]=hist(YAM2,10);
        clear fM fm PL_level PL_Time; fM=find(CNum==max(CNum));         
        PL_level=PLbin(fM(1));
        fm=Xn(YAM2>=PL_level);
        PL_Time=fm(1);
        if ((tPlM2-PL_Time)>tau)
            tPlM2=PL_Time; yPlM2=PL_level;
        end
    end
   
%     NCfitM2=corrcoef(YAM2,ygridM2);
%     NCfitM2HB=corrcoef(YAM2,ygridM2HB);
% %     chekFit{l-2,3}='exp';
% %     if(NCfitM2HB(1,2)>NCfitM2(1,2))
% %         clear ygridM2 PlgM2
% %         ygridM2=ygridM2HB;PlgM2=PlgM2HB;
% %         chekFit{l-2,3}='HB';
% %     end
    clear fm2 K2M2 K2LM2 TC2
    [K2M2,K2LM2]=demo2d(Xn,ygridM2); 
    if(isnan(K2M2)==1)
        TC2=Xn(1);
    else
        fm2=Xn(K2LM2==max(K2LM2)); fm_min2=Xn(K2LM2<=0.2*max(K2LM2));
        if (Magco>=5.5)
            fm2=Xn(K2LM2==max(K2LM2)); fm_min2=Xn(K2LM2<=0.1*max(K2LM2));
        end
        fm_nextMax2=find(fm_min2>=fm2(1)); fm_sol2=find(Xn>=fm_min2(fm_nextMax2(1)));
        TC2=Xn(fm_sol2(1));    
        if (Magco<5)
            clear ftemp; ftemp=find(Xn>=fm2(1));
            TC2=Xn(ftemp(1));    %max curvature            
        end
    end
    if (strcmp(chekFit{nuM,3},'HB'))==1
        clear K2M2 K2LM2 Xnc ygridM2c
        if(tau*(tD0M2(1))<Xn(1,end))
            Xnc=Xn(Xn>=tau*(tD0M2(1)-1));ygridM2c=ygridM2(Xn>=tau*(tD0M2(1)-1));
        else
            Xnc=Xn;ygridM2c=ygridM2;
        end
              
        [K2M2,K2LM2]=demo2d(Xnc,ygridM2c);
        if(isnan(K2M2)==1)
            TC2=Xn(1);
        else
            fm2=Xnc(K2LM2==max(K2LM2)); fm_min2=Xnc(K2LM2<=0.2*max(K2LM2));
            if (Magco>=5.5)
                fm2=Xnc(K2LM2==max(K2LM2)); fm_min2=Xnc(K2LM2<=0.1*max(K2LM2));
            end
            fm_nextMax2=find(fm_min2>=fm2(1)); fm_sol2=find(Xnc>=fm_min2(fm_nextMax2(1)));
            TC2=Xnc(fm_sol2(1));    
            if (Magco<5)
                clear ftemp; ftemp=find(Xnc>=fm2(1));
                TC2=Xnc(ftemp(1));    %max curvature
            end
        end
         clear fm ym K2L K2 
        if(tD0M2(1)~=1)            
            [K2,K2L]=demo2d(Xn(1,1:tD0M2(1)),smooth(ygridM2(1:tD0M2(1),1)));
            fm=Xn(K2L==max(K2L));
            ymM2=ygridM2(Xn==fm(end));
        else
%             tD0M2(1)=tD0M2(1)-1; 
            fm=tau*tD0M2(1);
            ymM2=ygridM2(1);
        end        
    end
    clear ftemp; ftemp=find(timeMax>=TC2);
    PL.UpTime(nuM,:)=[fm(end) tDM2 TC2 tPlM2];
    PL.UpAmp(nuM,:)=[ymM2 yDM2 ygridM2(ftemp(1)) yPlM2];

    scatter(PL.UpTime(nuM,4),PL.UpAmp(nuM,4),60,'c','fill','MarkerEdgeColor','k')
    scatter(PL.UpTime(nuM,3),PL.UpAmp(nuM,3),60,'dm','fill','MarkerEdgeColor','k')
%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear ftc fmin fmax Pl_Level LPDT_M1 LPDT_M2 Pl_Level_min Pl_Level_max; 
    Pl_Level=[max(ygrid);max(ygridM1);max(ygridM2)];
    coefE(9,nuM)=std(Pl_Level);
    .......................................................................
    coefE(10,nuM)=std([coefE(6,nuM);TC1;TC2]);    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear fDi; fDi=find( Xn>=coefE(6,nuM));SlogR(1,nuM)=Dis_mean(fDi(1),1);
    .......................................................................   
    ylabel('LPDT','FontSize',14); xlabel('Time (s)','FontSize',14); 
    xlim([0 inf]);
   %======================================================================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tSTF=0:tau:((10/tau)-1)*tau;tSTF=tSTF';
    tSTF_Q=0:tau:((20/tau)-1)*tau;tSTF_Q=tSTF_Q';
    STF=zeros(length(tSTF),1);STF_P=10.^(PL.MainAmp(nuM,4));
    clear fm; fm=find(tSTF<=PL.MainTime(nuM,4)); STF(fm,1)=(STF_P/PL.MainTime(nuM,4)).*tSTF(fm,1);
    clear fm; fm=find(tSTF>=PL.MainTime(nuM,4) & tSTF<=(2*PL.MainTime(nuM,4))); 
    STF(fm,1)=2.*STF_P-(STF_P/PL.MainTime(nuM,4)).*tSTF(fm,1);
    .......................................................................
    p = get(gca, 'Position');
    h = axes('Parent', gcf, 'Position', [p(1)+(1.7*p(3)/4) p(2)+(0.25*p(4)/4) 0.8*p(3)/4 0.8*p(4)/4]);
    box on;
    plot(Xn,nuMStaLPDT,'color','b','LineWidth',2);
    title('Number of Stations','FontSize',8)
    grid on; 
    
    h = axes('Parent', gcf, 'Position', [p(1)+(2.7*p(3)/4) p(2)+(0.25*p(4)/4) 1.2*p(3)/4 1.2*p(4)/4]);
    clear fm; fm=find(tSTF<=(2.*PL.MainTime(nuM,4)+(tau/2))); STF=STF(fm,1);tSTF=tSTF(fm,1);
    box on;grid on; hold on;
    plot(tSTF,STF,'color','k','LineWidth',2);
    title('Q Filter is Off')
    %===========================================================anelastic
    %attenuation using crows code
    global Vr Vp Vs Ap
    Vr=0.9*Vs;
    Ap=(2*RadP)./(4*pi()*rho*(Vp.^3));
    clear logM0t_Mode logM0t_Curve d_logM0t logM0tP_Mode logM0tM_Mode ...
        logM0tP_Curve logM0tM_Curve MagTe MagSe DeltaSc0
  
    logM0t_Mode=(PL.MainAmp(nuM,4)-log10(Ap./(PL.MainTime(nuM,4)-PL.MainTime(nuM,1))));
    logM0t_Curve=(PL.MainAmp(nuM,3)-log10(Ap./(PL.MainTime(nuM,3)-PL.MainTime(nuM,1))));
    
    d_logM0t=coefE(9,nuM)-((coefE(10,nuM).*((PL.MainTime(nuM,4)-PL.MainTime(nuM,1))))./log(10));
    logM0tP_Mode=logM0t_Mode+d_logM0t;
    logM0tM_Mode=logM0t_Mode-d_logM0t;
    logM0tP_Curve=logM0t_Curve+d_logM0t;
    logM0tM_Curve=logM0t_Curve-d_logM0t;
    
    if Magco<3
        logM0tP_Mode=logM0t_Mode+2.*d_logM0t;
        logM0tM_Mode=logM0t_Mode-2.*d_logM0t;
        logM0tP_Curve=logM0t_Curve+2.*d_logM0t;
        logM0tM_Curve=logM0t_Curve-2.*d_logM0t;
    end
    if (d_logM0t<0.75)
        MagTe=([min(logM0t_Mode-0.75, logM0t_Curve-0.75),min(logM0t_Mode+0.75, logM0t_Curve+0.75)]-9.1)./1.5;
    else
        MagTe=([min(logM0tM_Mode, logM0tM_Curve),min(logM0tP_Mode, logM0tP_Curve)]-9.1)./1.5;    
    end
    MagTe=sort(MagTe);
    if (MagTe(1)<0)
        MagSe=0.1:0.1:MagTe(2);
    else
        MagSe=MagTe(1):0.1:MagTe(2);
    end
    if(Magco>5)
        DeltaSc0=10.^(-2:0.1:2); 
    else
        DeltaSc0=10.^(-3:0.1:1); 
    end
    
    Initial=[];Final=[];numQ=1;
    if (strcmp(QFilter,'on')==1)
    for l1=1:length(MagSe)
        for l2=1:length(DeltaSc0)
            Mag=MagSe(1,l1);
            DeltaSc=DeltaSc0(1,l2);
            Values=SaharQEvents(Q,Vp,Vs,Vr,Ap,Mag,DeltaSc,DisR);
            Initial(numQ,:)=Values(1,:);
            Final(numQ,:)=Values(2,:);
            numQ=numQ+1;
        end
    end
    
    CostF1_M=((Final(:,4)-(10^(PL.MainAmp(nuM,4)))).^2)./Final(:,4)...
     +((Final(:,3)-PL.MainTime(nuM,3)).^2)./Final(:,3);
    indCostF1_M=find(CostF1_M==min(CostF1_M));
    QCorI_M(nuM,:)=Initial(indCostF1_M(1),:);
    QCorF_M(nuM,:)=Final(indCostF1_M(1),:);
    
     CostF1_C=((Final(:,4)-(10^(PL.MainAmp(nuM,3)))).^2)./Final(:,4)...
     +((Final(:,3)-PL.MainTime(nuM,3)).^2)./Final(:,3);
    indCostF1_C=find(CostF1_C==min(CostF1_C));
    QCorI_C(nuM,:)=Initial(indCostF1_C(1),:);
    QCorF_C(nuM,:)=Final(indCostF1_C(1),:);
    ..................................................................STF_Q
    STF_Q=zeros(length(tSTF_Q),1);
    clear fm; fm=find(tSTF_Q<=QCorI_M(nuM,3)); STF_Q(fm,1)=(QCorI_M(nuM,4)/QCorI_M(nuM,3)).*tSTF_Q(fm,1);
    clear fm; fm=find(tSTF_Q>=QCorI_M(nuM,3) & tSTF_Q<=(2*QCorI_M(nuM,3))); 
    STF_Q(fm,1)=2.*QCorI_M(nuM,4)-(QCorI_M(nuM,4)/QCorI_M(nuM,3)).*tSTF_Q(fm,1);
    .......................................................................
    clear fm; fm=find(tSTF_Q<=2.*QCorI_M(nuM,3)); STF_Q=STF_Q(fm,1);tSTF_Q=tSTF_Q(fm,1);
    hold on;
    plot(tSTF_Q,STF_Q,':k','LineWidth',2);   
    title('Q Filter is On')
    end
    
    saveas(figure (2),fullfile(path_Fig,[num2str(g(l).name)  '-2.png']));
    close figure 2 
    %===========================================================anelastic
    figure (3)
    
    ColorList=jet(7);
    subplot 121;box on;hold on; grid on
    title('Origin Curves');
    if(Magco<1)
        plot(Xn,Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn,ygrid,'color',ColorList(1,:),'LineWidth',2); 
    elseif(Magco>=1 && Magco<2)
        plot(Xn,Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn,ygrid,'color',ColorList(2,:),'LineWidth',2);
    elseif(Magco>=2 && Magco<3)
        plot(Xn,Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn,ygrid,'color',ColorList(3,:),'LineWidth',2);
    elseif(Magco>=3 && Magco<4)
        plot(Xn,Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn,ygrid,'color',ColorList(4,:),'LineWidth',2);
    elseif (Magco>=4 && Magco<5)
        plot(Xn,Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn,ygrid,'color',ColorList(5,:),'LineWidth',2);
    elseif (Magco>=5 && Magco<6)
        plot(Xn,Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn,ygrid,'color',ColorList(6,:),'LineWidth',2);
    else
        plot(Xn,Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn,ygrid,'color',ColorList(7,:),'LineWidth',2);
    end 
    scatter(PL.MainTime(nuM,4),PL.MainAmp(nuM,4),60,'c','fill','MarkerEdgeColor','k')
    scatter(PL.MainTime(nuM,3),PL.MainAmp(nuM,3),60,'dm','fill','MarkerEdgeColor','k')
    xlim([0 inf])
    subplot 122; box on;hold on; grid on
    title('Corrected curves');
    if(Magco<1)
        plot(Xn-PL.MainTime(nuM,1),Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn-PL.MainTime(nuM,1),ygrid,'color',ColorList(1,:),'LineWidth',2); 
    elseif(Magco>=1 && Magco<2)
        plot(Xn-PL.MainTime(nuM,1),Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn-PL.MainTime(nuM,1),ygrid,'color',ColorList(2,:),'LineWidth',2);
    elseif(Magco>=2 && Magco<3)
        plot(Xn-PL.MainTime(nuM,1),Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn-PL.MainTime(nuM,1),ygrid,'color',ColorList(3,:),'LineWidth',2);
    elseif(Magco>=3 && Magco<4)
        plot(Xn-PL.MainTime(nuM,1),Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn-PL.MainTime(nuM,1),ygrid,'color',ColorList(4,:),'LineWidth',2);
    elseif (Magco>=4 && Magco<5)
        plot(Xn-PL.MainTime(nuM,1),Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn-PL.MainTime(nuM,1),ygrid,'color',ColorList(5,:),'LineWidth',2);
    elseif (Magco>=5 && Magco<6)
        plot(Xn-PL.MainTime(nuM,1),Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn-PL.MainTime(nuM,1),ygrid,'color',ColorList(6,:),'LineWidth',2);
    else
        plot(Xn-PL.MainTime(nuM,1),Yn,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on;
        plot(Xn-PL.MainTime(nuM,1),ygrid,'color',ColorList(7,:),'LineWidth',2);
    end
    scatter(PL.MainTime(nuM,4)-PL.MainTime(nuM,1),PL.MainAmp(nuM,4),60,'c','fill','MarkerEdgeColor','k')
    scatter(PL.MainTime(nuM,3)-PL.MainTime(nuM,1),PL.MainAmp(nuM,3),60,'dm','fill','MarkerEdgeColor','k')
    xlim([0 inf])
    Mweve(nuM,1)=Magco;EveL{nuM,1}=g(l).name;
    Depth(nuM,1)=depCo;
    nuM=nuM+1;
    cd ..
%     pause
end
disp('**********************************************************')
disp([num2str(nuM-1),' events are evaluated among ',num2str(length(g)-2),' events'])

% title(Region);
subplot 121; xlabel('Time (s)');ylabel('LPDT (SI unit)'); box on
subplot 122; xlabel('Time (s)');ylabel('LPDT (SI unit)'); box on
saveas(gca,fullfile(path_Fig,'All.png'));
cd (path_Fig)
close
