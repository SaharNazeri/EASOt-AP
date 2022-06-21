    function [pks,locs,Yn_temp]=SMLPXT(Magco,times,Ys)
    clear pks locs
    [pks,locs]=findpeaks(smooth(Ys),times,'MinPeakDistance',times(end)/3);
    if(isempty(pks)==1);locs(1,1)=times(end);end
    if (locs(1,1)>times(end)/3)
       [pks,locs]=findpeaks(smooth(Ys),times,'MinPeakDistance',times(end)/4);     
    end    
    if (Magco>=5.5)
        clear pks locs
        [pks,locs]=findpeaks(smooth(Ys),times,'MinPeakDistance',times(end)/2);
        if(isempty(pks)==1);locs(1,1)=times(end);end
        if (locs(1,1)>times(end)/2)
        [pks,locs]=findpeaks(smooth(Ys),times,'MinPeakDistance',times(end)/3);     
        end        
    end
    if(times(end)<=1)
        clear pks locs
        [pks,locs]=findpeaks(smooth(Ys),times,'MinPeakDistance',times(end)/2);
        if(isempty(pks)==1);locs(1,1)=times(end);end
        if (locs(1,1)>times(end)/2)
            [pks,locs]=findpeaks(smooth(Ys),times,'MinPeakDistance',times(end)/3);         
        end 
        if(isempty(pks)==1);locs(1,1)=times(end);end
        if (locs(1,1)<times(end)/3)
            clear YsT timesT
            YsT=Ys(times>locs(1,1));
            timesT=times(times>locs(1,1));
            [pks,locs]=findpeaks(smooth(YsT),timesT,'MinPeakDistance',timesT(end)/2);
        end
    end
    clear Xn_temp Yn_temp
    if isempty(locs)==1 || locs(1,1)==times(end) 
        Xn_temp=times;
        Yn_temp=Ys;
    else
        Xn_temp=times(1,times<=(2*locs(1,1)));
        Yn_temp=Ys(times<=(2*locs(1,1)),1);
        Yn_temp(end+1:length(Ys),1)=Yn_temp(end,1);
        Yn_temp=Yn_temp';
    end    
    
    
    end