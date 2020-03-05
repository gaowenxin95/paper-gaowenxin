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

function [best,fmin,N_iter]=fpa_demo(para)
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
d=30;
Lb=-2*ones(1,d);
Ub=2*ones(1,d);
Cr=0.15; %初始化交叉概率 

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
          if rand>p,
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
                 S(i,:)=Sol(i,:)+rand(1,d).*(Sol(JK(1),:)-Sol(JK(2),:))+rand(1,d).*(best-Sol(i,:));    
              %%  交叉操作 
                r=rand(1,d); 
                Cr_I=find(r<Cr); 
                S(i,Cr_I)=Sol(JK(3),Cr_I); 

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

 semilogy(cg_curve,'Color','r')
title('Convergence curve')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid off
box on
legend('BOA')

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
% z=(1-u(1))^2+100*(u(2)-u(1)^2)^2+100*(u(3)-u(2)^2)^2;


%Ackely  
%  dim=30;
%  z=0;
%  for i=1:dim
%      z1=-20*exp(-0.2*sqrt(sum(u(i).^2)/dim))-exp(sum(cos(2*pi*u(i)))/dim)+20+exp(1);
%     z=z+z1;
%  end
%Apline 
%  dim=30;
%  z=0;
%  for i=1:dim
%     z1=abs(u(i)*sin(u(i)) + 0.1*u(i));
%     z=z+z1;
%  end

% % Bohachevsky 
%  dim=30;
%  z=0;
%  for i=1:dim
%     z= u(:,1).^2 + 2.*u(:,2).^2 - 0.3.*cos(3*pi.*u(:,1)) - 0.4.*cos(4*pi.*u(:,2)) + 0.7;
%     
%  end

%zakjarov 
%  dim=30;
%  z=0;
%  for i=1:dim
%     z1=sum(u(i).^2)+(1/2*sum(i*u(i))).^2+(1/2*sum(i*u(i))).^4;
%     z=z+z1;
% end

%Griewank 
% z=0;
% dim=30;
% for i=1:dim
%     z1=1+1/4000*sum(u(i).^2)-prod(cos(u(i)/sqrt(i)));
%     z=z1+z;
% end

% %matyas 
%  dim=30;
%  z=0;
%  for i=1:dim
%     z= 0.26.*(u(:,1).^2 + u(:,2).^2) - 0.48.*u(:,1).*u(:,2);
%     
%  end


% F7 Quartic function with noise
% dim=30;
% for i=1:dim
%      z=sum((1:dim).*(u(i).^4))+rand;
% end

% % schwefel 2.21
%  dim=30;
%  z=0;
%  for i=1:dim
%     z1=max(abs(u(i)));
%     z=z1;
%  end

%sphere 
 z=0;
 dim=30;
  for i=1:dim
     z1=sum(u(i).^2);
     z=z+z1;
 end

%sum square
% z=0;
% dim=30;
% for i=1:dim
%       z1=sum(i*u(i).^2);
%       z=z+z1;
% end

% f14 schaffer 
% j1=(sin(sqrt(u(:,1).^2+u(:,2).^2))).^2-0.5;
% j2=(1+0.001.*(u(:,1).^2+u(:,2).^2)).^2;
% z=j1./j2+0.5;

% powell
% dim=30;
% z1=0;
% for j=1:dim/4
%     z1=z1+ (u(:,4*j-3) + 10.*u(:,4*j-2)).^2 + 5.*(u(:,4*j-1) - u(:,4*j)).^2 + (u(:,4*j-2) - u(:,4*j-1)).^4 + 10.*(u(:,4*j-3) - u(:,4*j)).^4;
%     z=z1;
% end

% schwefel2.22
% dim=30;
% z=0;
% for i=1:dim
%     z1=sum(abs(u(i)))+prod(abs(u(i)));
%     z=z+z1;
% end

%rastraign 
%  dim=30;
%  z=0;
%  for i=1:dim
%      z1=sum(u(i).^2-10*cos(2*pi*u(i)))+10;
%      z=z+z1;
% end




 
%  F6 step
%  z=0;
%  dim=30;
% for i=1:dim
%     z1=sum((u(i)+0.5).^2);
%     z=z+z1
% end


% F18 Goldstein price
% z=(1+(u(1)+u(2)+1)^2*(19-14*u(1)+3*(u(1)^2)-14*u(2)+6*u(1)*u(2)+3*u(2)^2))*...
%     (30+(2*u(1)-3*u(2))^2*(18-32*u(1)+12*(u(1)^2)+48*u(2)-36*u(1)*u(2)+27*(u(2)^2)));


% %F28 power sum
% n = 4;
% b = [8,18,44,114];
% s = 0;
% for k = 1:n
%     z1 = 0;
%     for j = 1:n
%         z1 = z1 + u(:,j).^k;
%     end
%     s = s + (z1 - b(k)).^2;
% end
% z = s;


%sphere BAF9
% z=0;
 %dim=50;
 % for i=1:dim
      %for j=1:i
      %z2=0;
      %z1=sum(u(j));
      %z2=z2+z1;
      %end
      %z3=sum(z2.^2);
     % z=z+z3;      
 % end
 

% F13 Levy
% dim=30;
%  for i=1:dim
%     z1=.1*((sin(3*pi*u(1)))^2+sum((u(1:dim-1)-1).^2.*(1+(sin(3.*pi.*u(2:dim))).^2))+...
%         ((u(dim)-1)^2)*(1+(sin(2*pi*u(dim)))^2))+sum(Ufun(u,5,100,4));
%   end
% z=z1