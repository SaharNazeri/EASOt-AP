function [K2,K2L]=demo2d(t,y)
...........................................................................
clear PLdem YA XA L2 R2 K2 K2L;
Xn=t;
PLdem=[]; PLdem(:,2)=Xn;PLdem(:,5)=y;
...........................................................................
[B,I]=sort(PLdem(:,2)); 
PLdem=PLdem(I,:);
%  clear PL;PL=[]; PL(:,2)=Xn;PL(:,5)=ygrid;
...........................................................................
 window=1; YA(1)=PLdem(1,5);
 for m=2:floor(length(PLdem(:,5))/window)
     dd=PLdem((m-1)*window+1:m*window,5);
     YA(m)=max(dd); 
     if ( YA(m-1) >= YA(m) )
             YA(m)=YA(m-1);
     end
end
XA=PLdem(:,2)';
..........................................................................
x2 = XA;
% y2 = YA;
y2 = PLdem(:,5)';
X = [x2',y2'];
[L2,R2,K2] = curvature(X);
K2L(1,1:length(K2))=0.5.*sqrt(K2(:,1).^2+K2(:,2).^2);
fn=find(isnan(K2L)==1);K2L(1,fn)=0;
% figure;
% plot(L2,R2)
% title('Curvature radius vs. cumulative curve length')
% xlabel L
% ylabel R
% figure;
% h = scatter(PLdem(:,2),PLdem(:,5),10,K2L,'fill'); grid on; %axis equal
% set(h,'marker','.');
% xlabel x
% ylabel y
% title('2D curve with curvature vectors')
% hold on
% plot(x2,y2,'k','Linewidth',1.5)
% quiver(x2',y2',K2(:,1),K2(:,2));
% fm=find(K2L==max(K2L));
% scatter(PLdem(fm,2),PLdem(fm,5),100,'rd','fill','MarkerEdgeColor','k')
