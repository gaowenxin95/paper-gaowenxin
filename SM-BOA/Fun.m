
function z=Fun(u)

% F1 Sphere [-100,100]
%  z=sum(u.^2);

%F2 Schwel [-10,10]
%  z=sum(abs(u))+prod(abs(u));

%F3 [-100,100]
% dim=size(u,2);
% z=0;
% for i=1:dim
%   z=z+sum(u(1:i)-1)^2;
 % end


% F4 [-100,100]
% z=max(abs(u));

%F5 Rosenbrock  [-30,30]
%dim = 30;
%z=sum(100*(u(2:dim)-(u(1:dim-1).^2)).^2+(u(1:dim-1)-1).^2);


%F6 [-100,100]
%z=sum(abs((u+.5)).^2);

%F7 [-1.28,1.28]
%dim=size(u,2);
%z=sum([1:dim].*(u.^4))+rand;

%F8 [-500,500]
%z=sum(-u.*sin(sqrt(abs(u))));

%F9 Rastrign [-5.12,5.12]
% dim=size(u,2);
% z=sum(u.^2-10*cos(2*pi.*u))+10*dim;

%F10 Ackley [-32,32]
% dim=size(u,2);
% z=-20*exp(-.2*sqrt(sum(u.^2)/dim))-exp(sum(cos(2*pi.*u))/dim)+20+exp(1);

%F11 Griewank [-600,600]
% dim=size(u,2);
% z=sum(u.^2)/4000-prod(cos(u./sqrt([1:dim])))+1;

%F12  %[-50,50]
% dim=size(u,2);
% z=(pi/dim)*(10*((sin(pi*(1+(u(1)+1)/4)))^2)+sum((((u(1:dim-1)+1)./4).^2).*...
% (1+10.*((sin(pi.*(1+(u(2:dim)+1)./4)))).^2))+((u(dim)+1)/4)^2)+sum(Ufun(u,10,100,4));
% function o=Ufun(u,a,k,m)
% o=k.*((u-a).^m).*(u>a)+k.*((-u-a).^m).*(u<(-a));


%F13  %[-50,50]
%dim=size(u,2);
%z=.1*((sin(3*pi*u(1)))^2+sum((u(1:dim-1)-1).^2.*(1+(sin(3.*pi.*u(2:dim))).^2))+...
%((u(dim)-1)^2)*(1+(sin(2*pi*u(dim)))^2))+sum(Ufun(u,5,100,4));
%function o=Ufun(u,a,k,m)
%o=k.*((u-a).^m).*(u>a)+k.*((-u-a).^m).*(u<(-a));

%F14  [-65.536,-65.536]
%aS=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;,...
   % -32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];
%for j=1:25
  % bS(j)=sum((u'-aS(:,j)).^6);
%end
%z=(1/500+sum(1./([1:25]+bS))).^(-1);

% F15   %[-5,5]
%aK=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
%bK=[.25 .5 1 2 4 6 8 10 12 14 16];bK=1./bK;
%z=sum((aK-((u(1).*(bK.^2+u(2).*bK))./(bK.^2+u(3).*bK+u(4)))).^2);
%end

%F16    %[-5,5]
  %z=4*(u(1)^2)-2.1*(u(1)^4)+(u(1)^6)/3+u(1)*u(2)-4*(u(2)^2)+4*(u(2)^4);

  % F17   %[-5,5]
%z=(u(2)-(u(1)^2)*5.1/(4*(pi^2))+5/pi*u(1)-6)^2+10*(1-1/(8*pi))*cos(u(1))+10;



% F18    %[-5,5]
%z=(1+(u(1)+u(2)+1)^2*(19-14*u(1)+3*(u(1)^2)-14*u(2)+6*u(1)*u(2)+3*u(2)^2))*...
 % (30+(2*u(1)-3*u(2))^2*(18-32*u(1)+12*(u(1)^2)+48*u(2)-36*u(1)*u(2)+27*(u(2)^2)));


% F19     %[1,3]
%aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
%pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
%z=0;
%for i=1:4
    %z=z-cH(i)*exp(-(sum(aH(i,:).*((u-pH(i,:)).^2))));
%end


% F20       %[0,1]
%aH=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
%cH=[1 1.2 3 3.2];
%pH=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
%.2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
%z=0;
%for i=1:4
 % z=z-cH(i)*exp(-(sum(aH(i,:).*((u-pH(i,:)).^2))));
%end

% F21      %[0,10]
%aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
%cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
%z=0;
%for i=1:5
  %z=z-((u-aSH(i,:))*(u-aSH(i,:))'+cSH(i))^(-1);
%end


%F22      %[0,10]
% aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
% cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
% z=0;
% for i=1:7
%     z=z-((u-aSH(i,:))*(u-aSH(i,:))'+cSH(i))^(-1);
% end


% F23      %[0,10]
%aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
%cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
%z=0;
%for i=1:10
    %z=z-((u-aSH(i,:))*(u-aSH(i,:))'+cSH(i))^(-1);
%end




%F����[-10��10]

%x=u(1,1);  
%y=u(1,2);  
%temp=x^2+y^2;  
%z=0.5+(sin(sqrt(temp))^2-0.5)/(1+0.001*temp)^2;
%for i=1:30
  % dim=size(u,2);
   %temp= i-1/dim-1;
   %  z= sum((10.^(6*temp)).* u.^2);      
    
%end   
   
 
%F25 Apline
% function o = F25(x)
% dim = size(u,2);
% s = 0;
% for i = 1:dim
%     s = s + abs(u(i)*sin(u(i)) + 0.1*u(i));
% end
% z = s;
% end


%F26 Schaffer
% function z= F26(u)
% z = 0.5 + (sin(sqrt(u(:,1).^2 + u(:,2).^2)).^2 - 0.5) ./ (1 + 0.001.*(u(:,1).^2 + u(:,2).^2)).^2;
% end


%F27 powell
% function z = F27(x)
% dim = size(u,2);
% s= 0;
% for j = 1:dim/4
%     s = s + (u(:,4*j-3) + 10.*u(:,4*j-2)).^2 + 5.*(u(:,4*j-1) - u(:,4*j)).^2 + (u(:,4*j-2) - u(:,4*j-1)).^4 + 10.*(u(:,4*j-3) - u(:,4*j)).^4;
% end
% z = s;
% end

%F28 power sum
%     function z = F28(x)
% n = 4;
% b = [8,18,44,114];
% s = 0;
% for k = 1:n
%     f1 = 0;
%     for j = 1:n
%         f1 = f1 + u(:,j).^k;
%     end
%     s = s + (f1 - b(k)).^2;
% end
% z = s;
% end

%F29 Quartic
%     function z = F29(u)
% dim = size(u,2);
% s = 0;
% for j = 1:dim
%     s = s + j.*u(:,j).^4;
% end
% z = s + rand(1);
% end

%F30 Shubert
% function z = F30(x)
% f1 = 0;
% f2 = 0;
% for j = 1:5
%     f1 = f1 + j.* cos((j + 1).*u(:,1)+j);
%     f2 = f2 + j.* cos((j + 1).*u(:,2)+j);
% end
% z = f1.*f2;
% end

%F31 Sum Squares
%     function z= F31(x)
% dim = size(u,2);
% s = 0;
% for j = 1:dim
%     s = s + j*u(:,j).^2;
% end
% z = s;
% end

%F32 Beale
%     function z = F32(x)
%  %%Beale - [-4.5,4.5] - Dim2;
% z = (1.5 - u(:,1) + u(:,1).*u(:,2)).^2 + (2.25 - u(:,1) +u(:,1).*u(:,2).^2).^2 + (2.625 - u(:,1) + u(:,1).*u(:,2).^3).^2;
% end


%F33 Bohachevsky
%     function z= F33(x)
%  %Bohachevsky1 - [-100,100] - Dim2
% z= u(:,1).^2 + 2.*u(:,2).^2 - 0.3.*cos(3*pi.*u(:,1)) - 0.4.*cos(4*pi.*u(:,2)) + 0.7;
% end

%F34 Michalewicz
%     function z =F34(x)
% %- Michalewicz - [0,pi] - Dim10
% dim = size(u,2);
% s = 0;
% for j = 1:dim
%     s = s + sin(u(:,j)).*sin(j.*u(:,j).^2 / pi).^20;
% end
% z= -s;
% end

%F35 Matyas
   z = 0.26.*(u(:,1).^2 + u(:,2).^2) - 0.48.*u(:,1).*u(:,2);   


%F36 Booth
%     function z =F36(x)
% % - Booth - [-10,10] - Dim2
% z = (u(:,1) + 2.*u(:,2)- 7).^2 + (2.*u(:,1) + u(:,2) - 5).^2;
% end

%F37  Shekel
%     function z = F37(x)
% % - Shekel7 - [0,10] - Dim4
% dim = size(u,2);

% m = 7;
% a = ones(10,4);
% a(1,:) = 4.0*a(1,:);
% a(2,:) = 1.0*a(2,:);
% a(3,:) = 8.0*a(3,:);
% a(4,:) = 6.0*a(4,:);
% for j = 1:2
%     a(5,2*j-1) = 3.0; a(5,2*j) = 7.0;
%     a(6,2*j-1) = 2.0; a(6,2*j) = 9.0;
%     a(7,j)     = 5.0; a(7,j+2) = 3.0;
%     a(8,2*j-1) = 8.0; a(8,2*j) = 1.0;
%     a(9,2*j-1) = 6.0; a(9,2*j) = 2.0;
%     a(10,2*j-1)= 7.0; a(10,2*j)= 3.6;
% end
% c(1) = 0.1; c(2) = 0.2; c(3) = 0.2; c(4) = 0.4; c(5) = 0.4;
% c(6) = 0.6; c(7) = 0.3; c(8) = 0.7; c(9) = 0.5; c(10)= 0.5;
% s = 0;
% for j = 1:m
%     f1 = 0;
%     for i = 1:dim
%         f1 = f1 + (u(:,i) - a(j,i)).^2;
%     end
%     s = s + 1 ./ (f1 + c(j));
% end
% z = -s;
% end

%F38 Trid6
%     function  z =F38(x)
% % Trid6 - [-36,36] - Dim6
% dim = size(x,2);
% f1 = 0;
% f2 = 0;
% for j = 1:dim
%     f1 = f1 + (u(:,j) - 1).^2;
% end
% for j = 2:dim
%     f2 = f2 + u(:,j).*u(:,j-1);
% end
% z = f1 - f2;
% end
% 
% 
% %F39 Perm
% function  o= F39(x)
% % Perm - [-D,D] - Dim4
% n = 4;
% b = 0.5;
% s = 0;
% for k = 1:n
%     f1 = 0;
%     for j = 1:n
%         f1 = f1 + (j^k + b).*((u(:,j)./ j).^k - 1);
%     end
%     s = s + f1.^2;
% end
% o = s;
% end

% F40 Zakharov
%     function z = F40(x)
% % - Zakharov - [-5,10] - Dim10
%  dim = size(u,2);
% f1 = 0;
% f2 = 0;
% f3 = 0;
% for j = 1:dim
%     f1 = f1 + u(:,j).^2;
%     f2 = f2 + 0.5.*j.*u(:,j);
%     f3 = f3 + 0.5.*j.*u(:,j);
% end
% z = f1 + f2.^2 + f3.^4;

