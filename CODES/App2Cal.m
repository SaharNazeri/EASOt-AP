function [M0p6,DeltaS6,a6,TX]=App2Cal(x0,Pl,TC)
Vp=5000; Vs=Vp/sqrt(3);Vr=0.9*Vs;
TX=2.*TC; PX=10.^(Pl);

a6=log10(TX./((1/Vr)+(2/(pi()*Vp))));
if (strcmp(x0,'D')==1)
    CX6=((2*0.52)*2)./(4*pi()*2700*(Vp^3)*(TX));
else
    CX6=((2*0.52)*4)./(4*pi()*2700*(Vp^3)*(TX.^2));
end
M0p6=PX./CX6;
DeltaS6=(7.*M0p6)./(16.*((10.^a6).^3));
DeltaS6=DeltaS6./(1e6);
M0p6=(log10(M0p6)-9.1)/1.5; %Mag
end
