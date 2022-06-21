function Ymax=FMax(Yin)
if size(Yin,1)==1;Yin=Yin';end 

clear Ymax0 window Ymax dd; 
Ymax0=Yin; window=1; Ymax(1)=Ymax0(1,1);
    for m=2:floor(length(Ymax0)/window)
        dd=Ymax0((m-1)*window+1:m*window,1);
        Ymax(m)=max(dd); 
        if ( Ymax(m-1) >= Ymax(m) )
            Ymax(m)=Ymax(m-1);
        end
    end
    
end