function [popt]=App1Val(Mw,LPXTin,SD1,SD2, SDN)
global path_Reg
global Vr Vp Ap
if size(Mw,1)==1;Mag=Mw';else; Mag=Mw; end
if size(LPXTin,1)==1;LPXT=LPXTin';else; LPXT=LPXTin; end  

stress_values = linspace(SD1,SD2, SDN);
for j=1:length(stress_values)
StressDrop=stress_values(j)*1e6;  %Pa
K=((7/16)^(1/3)) .* ((1./StressDrop).^(1/3)) .* ((1./Vr)-(2./(pi().*Vp)));
logM0=(3/2).*(LPXT-log10(Ap./K));

dobs=1.5*Mag+9.1;
MagND=logM0; 
     [p,S] = polyfit(dobs,MagND,1);
     popt(j,1)= p(1);
     popt(j,2) = p(2);
     popt(j,3) = StressDrop/1e6; %MPa
     [y_fit,delta] = polyval(p,dobs,S);
     y_diff = y_fit-dobs;
     RMS = rms(y_diff);
     popt(j,4) = RMS;

end
end