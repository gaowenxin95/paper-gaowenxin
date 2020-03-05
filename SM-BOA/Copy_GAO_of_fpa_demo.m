% --------------------------------------------------------------------% 
% Flower pollenation algorithm (FPA), or flower algorithm             %
% Programmed by Xin-She Yang @ May 2012                               % 
% --------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: This demo program contains the very basic components of      %
% the flower pollination algorithm (FPA), or flower algorithm (FA),   % 
% for single objective optimization.    It usually works well for     %
% unconstrained functions only. For functions/problems with           % 
% limits/bounds and constraints, constraint-handling techniques       %
% should be implemented to deal with constrained problems properly.   %
%                                                                     %   
% Citation details:                                                   % 
%1)Xin-She Yang, Flower pollination algorithm for global optimization,%
% Unconventional Computation and Natural Computation,                 %
% Lecture Notes in Computer Science, Vol. 7445, pp. 240-249 (2012).   %
%2)X. S. Yang, M. Karamanoglu, X. S. He, Multi-objective flower       %
% algorithm for optimization, Procedia in Computer Science,           %
% vol. 18, pp. 861-868 (2013).                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [best,fmin,N_iter]=Copy_GAO_of_fpa_demo(para)
% Default parameters
if nargin<1,
   para=[30 0.8];
end
power_exponent=0.1;
sensory_modality=0.01;
n=para(1);           % Population size, typically 10 to 25
p=para(2);           % probabibility switch

% Iteration parameters
N_iter=500;            % Total number of iterations

% Dimension of the search variables
d=500;
Lb=-2*ones(1,d);
Ub=2*ones(1,d);
% Cr=0.15; %初始化交叉概率 

% Initialize the population/solutions
for i=1:n,
  Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
  Fitness(i)=Fun(Sol(i,:));
end

% Find the current best
[fmin,I]=min(Fitness);
best=Sol(I,:);
S=Sol; 

% Start the iterations -- Flower Algorithm 
for t=1:N_iter,
        % Loop over all bats/solutions
        for i=1:n,
          % Pollens are carried by insects and thus can move in
          % large scale, large distance.
          % This L should replace by Levy flights  
          % Formula: x_i^{t+1}=x_i^t+ L (x_i^t-gbest)
          Fnew=Fun(S(i,:));
          FP=(sensory_modality*(Fnew^power_exponent));   
          if rand<p,
          %% L=rand;
%           L=Levy(d);
%           dS=L.*(Sol(i,:)-best);
%           S(i,:)=Sol(i,:)+dS;
          dS = rand * rand * best - Sol(i,:);        %Eq. (2) in paper
          S(i,:)=Sol(i,:)+dS*FP;
          [sort_arr,ind_sm]=sort(Fitness); 
          xs=S(i,:); 
          xg=Sol(ind_sm(1),:);%  最优蝴蝶的位置 
          xb=Sol(ind_sm(2),:);%  次优蝴蝶的位置 
          xc=(xg+xb)./2;%  最优蝴蝶和次优蝴蝶位置的中心点    
                 xr=xc+(xc-xs); 
                 f_xg=Fun(xg); 
                 f_xs=Fun(xs); 
                 f_xr=Fun(xr); 
                 if f_xr<f_xg 
                     xe=xc+5*(xr-xc); 
                     f_xe=Fun(xe); 
                     if f_xe<f_xg 
                        S(i,:)=xe; 
                     else 
                        S(i,:)=xr; 
                     end 
                 end 
                 if f_xr>f_xs 
                     xt=xc+0.2*(xs-xc); 
                    f_xt=Fun(xt); 
                     if f_xt<f_xs 
                         S(i,:)=xt; 
                     end 
                 end 
                 if f_xs>f_xr>f_xg 
                     xw=xc-0.2*(xs-xc); 
                     f_xw=Fun(xw); 
                     if f_xw<f_xs 
                         S(i,:)=xw; 
                     else 
                         S(i,:)=xr; 
                     end 
                 end                 
          
          % Check if the simple limits/bounds are OK
          S(i,:)=simplebounds(S(i,:),Lb,Ub);
          
          % If not, then local pollenation of neighbor flowers 
          else
              epsilon=rand;
              % Find random flowers in the neighbourhood
%               JK=randperm(n);
%               % As they are random, the first two entries also random
%               % If the flower are the same or similar species, then
%               % they can be pollenated, otherwise, no action.
%               % Formula: x_i^{t+1}+epsilon*(x_j^t-x_k^t)
%               S(i,:)=S(i,:)+epsilon*(Sol(JK(1),:)-Sol(JK(2),:));
                 JK=randperm(n); 
                 %精英算子
%                  S(i,:)=Sol(i,:)+rand(1,d).*(Sol(JK(1),:)-Sol(JK(2),:))+rand(1,d).*(best-Sol(i,:));    
%               %%  交叉操作 
%                 r=rand(1,d); 
%                 Cr_I=find(r<Cr); 
%                 S(i,Cr_I)=Sol(JK(3),Cr_I); 
               dS=epsilon*epsilon*Sol(JK(1),:)-Sol(JK(2),:);
               S(i,:)=Sol(i,:)+dS*FP;                         %Eq. (3) in paper
        
          end
           % SCA optimize BOA
          for i = 1:n
            a = 2;
            r1=a-t*((a)/N_iter);
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();
            if r4<0.5
                
                S(i,:)= S(i,:)+(r1*sin(r2)*abs(r3*best-S(i,:)));
            else
                
                S(i,:)= S(i,:)+(r1*cos(r2)*abs(r3*best-S(i,:)));
            end
            
          
              % Check if the simple limits/bounds are OK
              S(i,:)=simplebounds(S(i,:),Lb,Ub);
          end
          
          % Evaluate new solutions
           Fnew=Fun(S(i,:));
          % If fitness improves (better solutions found), update then
            if (Fnew<=Fitness(i)),
                Sol(i,:)=S(i,:);
                Fitness(i)=Fnew;
           end
           
          % Update the current global best
          if Fnew<=fmin,
                best=S(i,:)   ;
                fmin=Fnew   ;
          end
        end
        % Display results every 100 iterations
        if round(t/100)==t/100,
        best
        fmin
        end
        cg_curve(t)=fmin
end

% Output/display
disp(['Total number of evaluations: ',num2str(N_iter*n)]);
disp(['Best solution=',num2str(best),'   fmin=',num2str(fmin)]);

semilogy(cg_curve,'r','LineWidth',2)
title('Convergence curve')
xlabel('Iteration');
ylabel('Best score obtained so far');
hold on
axis tight
grid off
box on
legend('PSO','WOA','SCA','BOA','SMSCABOA')

% disp(['Total number of evaluations: ',num2str(N_iter*n)]);
disp(['Best solution=',num2str(best),'   fmin=',num2str(fmin)]);



% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;

function y=sensory_modality_NEW(x,Ngen)
y=x+(0.025/(x*Ngen));

% Objective function and here we used Rosenbrock's 3D function
function z=Fun(u)



% F1 Sphere [-100,100]
%  z=sum(u.^2);

%F2 Schwel [-10,10]
%  z=sum(abs(u))+prod(abs(u));

%F3 [-100,100]
%dim=size(u,2);
%z=0;
%for i=1:dim
  %z=z+sum(u(1:i)-1)^2;
 % end


% F4 [-100,100]
% z=max(abs(u));

%F5 Rosenbrock  [-30,30]
%dim = 30;
%z=sum(100*(u(2:dim)-(u(1:dim-1).^2)).^2+(u(1:dim-1)-1).^2);


%F6 [-100,100]
%z=sum(abs((u+.5)).^2);

%F7 [-1.28,1.28]
% dim=size(u,2);
% z=sum([1:dim].*(u.^4))+rand;

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
dim=size(u,2);
z=(pi/dim)*(10*((sin(pi*(1+(u(1)+1)/4)))^2)+sum((((u(1:dim-1)+1)./4).^2).*...
(1+10.*((sin(pi.*(1+(u(2:dim)+1)./4)))).^2))+((u(dim)+1)/4)^2)+sum(Ufun(u,10,100,4));
function o=Ufun(u,a,k,m)
o=k.*((u-a).^m).*(u>a)+k.*((-u-a).^m).*(u<(-a));


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




%F其他[-10，10]

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
%   z = 0.26.*(u(:,1).^2 + u(:,2).^2) - 0.48.*u(:,1).*u(:,2);  


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
% dim = size(u,2);
% f1 = 0;
% f2 = 0;
% f3 = 0;
% for j = 1:dim
%     f1 = f1 + u(:,j).^2;
%     f2 = f2 + 0.5.*j.*u(:,j);
%     f3 = f3 + 0.5.*j.*u(:,j);
% end
% z = f1 + f2.^2 + f3.^4;
% end
