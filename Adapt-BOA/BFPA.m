 %%---------Flower Pollination Algorithm for continuous optimization---%%
function [best,fmin,N_iter]=fpa_demo(para)
% ��ʼ������
if nargin<1
    para=[30 0.8]; %��Ⱥ��n,ת������p
end
n=para(1); % ��Ⱥ��С��10-25��ȡ����ĵ�һ����?
p=para(2); % ���ָ��ʣ�ȡ����ĵڶ�����?
N_iter=500; % ����������
d=30; % �����ռ���߱�����ά�������Եĺ�����ά��?
Lb=-600*ones(1,d); % �����½�
Ub=600*ones(1,d); % �����Ͻ�
for i=1:n
    Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);%�������һ��ȫ�����Ž⣬Sol(i,:)�洢����ֵ
    Fitness(i)=Fun(Sol(i,:)); %���㺯����Ӧ��ֵ
end
[fmin,I]=min(Fitness);
best=Sol(I,:);
S=Sol;
for t=1:N_iter,
    for i=1:n
        if rand<p,
            %��ʽx_i(t+1)=x_i(t)+L*(x_i(t)-gbest)�컨�ڷ�%
            L=Levy(d);
            dS=L.*(best-Sol(i,:));
            S(i,:)=Sol(i,:)+dS;
           % S(i,:)=S(i,:).*(1+(1-t/799)*trnd(800));
            S(i,:)=simplebounds(S(i,:),Lb,Ub);
        else
            %��ʽx_i(t+1)+epsilon*(x_j(t)-x_k(t))ͬ���ڷ�%
            epsilon=rand;
            JK=randperm(n);
            S(i,:)=Sol(i,:)+epsilon*(Sol(JK(1),:)-Sol(JK(2),:));
            S(i,:)=simplebounds(S(i,:),Lb,Ub);
        end
        Fnew=Fun(S(i,:));
        if(Fnew<=Fitness(i)),
            Sol(i,:)=S(i,:);
            Fitness(i)=Fnew;
        end
        if Fnew<=fmin,
            best=S(i,:);
            fmin=Fnew;
        end
    end
    if round(t)==t,%ÿ��100�ε�����ӡ���������?
        best
        fmin
    end
    result(t)=fmin;
end

disp(['Total number of evaluations:',num2str(N_iter)]);
disp(['Best solution=',num2str(best)]);
disp(['fmin=',num2str(fmin)]);

hold on
semilogy(1:N_iter,result,'-.m')
%plot(1:N_iter,result,'--r')
%hold on
%set(gcf,'position',[5,5,13,11])%%����ͼƬ��С
%set(gca,'FontSize',6); % �������ִ�С��ͬʱӰ���������ע��ͼ���������
% xlabel('iteration');
% ylabel('lg(Ŀ�꺯��ֵ)');
title('Convergence curve')
xlabel('Iteration');
ylabel('Best score obtained so far');
%ylabel('Ŀ�꺯��ֵ');
legend('CWBOA','BOA','WOA','FPA');

function s=simplebounds(s,Lb,Ub)
ns_tmp=s;
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);
s=ns_tmp;
%%---levy fights---%%
function L=Levy(d)
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
L=0.01*step;

function z=Fun(u)





%Ackely  
%  dim=30;
%  z=0;
%  for i=1:dim
%      z1=-20*exp(-0.2*sqrt(sum(u(i).^2)/dim))-exp(sum(cos(2*pi*u(i)))/dim)+20+exp(1);
%     z=z+z1;
%  end

%Apline 
% %  dim=30;
% %  z=0;
% %  for i=1:dim
% %     z1=abs(u(i)*sin(u(i)) + 0.1*u(i));
% %     z=z+z1;
% %  end

% % Bohachevsky 
%  dim=2;
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
%  z=0;
%  dim=30;
%   for i=1:dim
%      z1=sum(u(i).^2);
%      z=z+z1;
%  end

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
dim=30;
z=0;
for i=1:dim
    z1=sum(abs(u(i)))+prod(abs(u(i)));
    z=z+z1;
end

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























%z=abs(u(1).^u(2)+u(2).^u(1)-5*u(1)*u(2)*u(3)-85)+abs(u(1).^3-u(2).^u(3)-u(3).^u(2)-60)+abs(u(1).^u(3)+u(3).^u(1)-u(2)-2)
%z=abs(u(1).^2+u(2).^2+u(3).^2-3)+abs(u(1).^2+u(2).^2+u(1)*u(2)+u(1)+u(2)-5)+abs(u(1)+u(2)+u(3)-3)
% schwel F2 
% dim=30;
% z=0;
% for i=1:dim
%     z1=sum(abs(u(i)))+prod(abs(u(i)));
%     z=z+z1;
% end

%rastraign BAF7
%  dim=30;
%  z=0;
%  for i=1:dim
%      z=sum(u(i).^2-10*cos(2*pi*u(i)))+10; 
%  end

%apline 
%  dim=30;
%  z=0;
%  for i=1:dim
%     z1=abs(u(i)*sin(u(i)) + 0.1*u(i));
%     z=z+z1;
%  end

%Bohachevsky 
%  dim=2;
%  z=0;
%  for i=1:dim
%     z= u(:,1).^2 + 2.*u(:,2).^2 - 0.3.*cos(3*pi.*u(:,1)) - 0.4.*cos(4*pi.*u(:,2)) + 0.7;
%     
%  end

 %F4
%  z=0;
%  dim=30;
%  for i=1:dim
%      z1=max(abs(u(i)))
%      z=z+z1
%  end





%F7 Quartic function with noise
% z=0;
% dim=30;
% for i=1:dim
%     z1=sum(i*(u(i).^4))+rand;
% end
% z=z+z1

%Griewank BAF9
% z=0;
% dim=30;
% for i=1:dim
%     z1=1+1/4000*sum(u(i).^2)-prod(cos(u(i)/sqrt(i)));
%     z=z1+z;
% end

%rosenbrock BAF5
%  dim=30;
%  z=0;
%  for i=1:(dim-1)
%     z1=sum(100*(u(i+1)-u(i).^2).^2+(u(i)-1).^2);
%     z=z+z1;
%  end

%sum square
%  dim=30;
%  z=0;
%  for i=1:dim
%     z1=sum((u(i).^2*i));
%     z=z+z1;
%  end

%schwefel 2.26
%  dim=30;
%  z=0;
%  for i=1:dim
%     z1=max(abs(u(i)));
%     z=z1;
%  end

% %Bohachevsky 
%  dim=30;
%  z=0;
%  for i=1:dim
%     z= 0.26.*(u(:,1).^2 + u(:,2).^2) - 0.48.*u(:,1).*u(:,2);
%     
%  end




%ackely  BAF8
%  dim=30;
%  z=0;
%  for i=1:dim
%      z1=-20*exp(-0.2*sqrt(sum(u(i).^2)/dim))-exp(sum(cos(2*pi*u(i)))/dim)+20+exp(1);
%     z=z+z1;
%  end

% sphere BAF3
%  z=0;
%  dim=30;
%   for i=1:dim
%      z1=sum(u(i).^2);
%      z=z+z1;
%  end

%zakjarov BAF7
%  dim=30;
%  z=0;
%  for i=1:dim
%     z1=sum(u(i).^2)+(1/2*sum(i*u(i))).^2+(1/2*sum(i*u(i))).^4;
%     z=z+z1;
% end

% F18 Goldstein price

% z=(1+(u(1)+u(2)+1)^2*(19-14*u(1)+3*(u(1)^2)-14*u(2)+6*u(1)*u(2)+3*u(2)^2))*...
%     (30+(2*u(1)-3*u(2))^2*(18-32*u(1)+12*(u(1)^2)+48*u(2)-36*u(1)*u(2)+27*(u(2)^2)));
% 


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


%sphere BAF8
% z=0;
% dim=30;
% for i=1:dim
%       z1=sum(i*u(i).^2);
%       z=z+z1;
% end

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
 
% % F7 Quartic function with noise
% dim=30;
% for i=1:dim
%      z=sum((1:dim).*(u(i).^4))+rand;
% end

% % F13 Levy
% dim=2;
%  for i=1:dim
%      z=.1*((sin(3*pi*u(1)))^2+sum((u(1:dim-1)-1).^2.*(1+(sin(3.*pi.*u(2:dim))).^2))+...
%         ((u(dim)-1)^2)*(1+(sin(2*pi*u(dim)))^2))+sum(Ufun(u,5,100,4));
%   end


%F29 Quartic

% dim = 30;
% z1 = 0;
% for j = 1:dim
%     z1= z1 + j.*u(:,j).^4;
% end
% z= z1 + rand(1);




%F1 sphere D=30 need
%j1=u(1).^2+u(2).^2+u(3).^2+u(4).^2+u(5).^2+u(6).^2+u(7).^2+u(8).^2+u(9).^2+u(10).^2;
%j2=u(11).^2+u(12).^2+u(13).^2+u(14).^2+u(15).^2+u(16).^2+u(17).^2+u(18).^2+u(19).^2+u(20).^2;
%j3=u(21).^2+u(22).^2+u(23).^2+u(24).^2+u(25).^2+u(26).^2+u(27).^2+u(28).^2+u(29).^2+u(30).^2;
%z=j1+j2+j3;

%f9 rastraign d=30
% j1=(u(1).^2-10.*cos(2*pi.*u(1)))+(u(2).^2-10.*cos(2*pi.*u(2)))+(u(3).^2-10.*cos(2*pi.*u(3)))+(u(4).^2-10.*cos(2*pi.*u(4)))+(u(5).^2-10.*cos(2*pi.*u(5)))+(u(6).^2-10.*cos(2*pi.*u(6)))+(u(7).^2-10.*cos(2*pi.*u(7)))+(u(8).^2-10.*cos(2*pi.*u(8)))+(u(9).^2-10.*cos(2*pi.*u(9)))+(u(10).^2-10.*cos(2*pi.*u(10)));
% j2=(u(11).^2-10.*cos(2*pi.*u(11)))+(u(12).^2-10.*cos(2*pi.*u(12)))+(u(13).^2-10.*cos(2*pi.*u(13)))+(u(14).^2-10.*cos(2*pi.*u(14)))+(u(15).^2-10.*cos(2*pi.*u(15)))+(u(16).^2-10.*cos(2*pi.*u(16)))+(u(17).^2-10.*cos(2*pi.*u(17)))+(u(18).^2-10.*cos(2*pi.*u(18)))+(u(19).^2-10.*cos(2*pi.*u(19)))+(u(20).^2-10.*cos(2*pi.*u(20)));
% j3=(u(21).^2-10.*cos(2*pi.*u(21)))+(u(22).^2-10.*cos(2*pi.*u(22)))+(u(23).^2-10.*cos(2*pi.*u(23)))+(u(24).^2-10.*cos(2*pi.*u(24)))+(u(25).^2-10.*cos(2*pi.*u(25)))+(u(26).^2-10.*cos(2*pi.*u(26)))+(u(27).^2-10.*cos(2*pi.*u(27)))+(u(28).^2-10.*cos(2*pi.*u(28)))+(u(29).^2-10.*cos(2*pi.*u(29)))+(u(30).^2-10.*cos(2*pi.*u(30)));
% z=300+j1+j2+j3;

% F10 ackely d=30
% j1=-20*exp(-0.2*sqrt(1/30*(u(1).^2+u(2).^2+u(3).^2+u(4).^2+u(5).^2+u(6).^2+u(7).^2+u(8).^2+u(9).^2+u(10).^2+u(11).^2+u(12).^2+u(13).^2+u(14).^2+u(15).^2+u(16).^2+u(17).^2+u(18).^2+u(19).^2+u(20).^2+u(21).^2+u(22).^2+u(23).^2+u(24).^2+u(25).^2+u(26).^2+u(27).^2+u(28).^2+u(29).^2+u(30).^2)));
% j2=-exp(1/30*(cos(2*pi*u(1))+cos(2*pi*u(2))+cos(2*pi*u(3))+cos(2*pi*u(4))+cos(2*pi*u(5))+cos(2*pi*u(6))+cos(2*pi*u(7))+cos(2*pi*u(8))+cos(2*pi*u(9))+cos(2*pi*u(10))+cos(2*pi*u(11))+cos(2*pi*u(12))+cos(2*pi*u(13))+cos(2*pi*u(14))+cos(2*pi*u(15))+cos(2*pi*u(16))+cos(2*pi*u(17))+cos(2*pi*u(18))+cos(2*pi*u(19))+cos(2*pi*u(20))+cos(2*pi*u(21))+cos(2*pi*u(22))+cos(2*pi*u(23))+cos(2*pi*u(24))+cos(2*pi*u(25))+cos(2*pi*u(26))+cos(2*pi*u(27))+cos(2*pi*u(28))+cos(2*pi*u(29))+cos(2*pi*u(30))));
% z=j1+j2+20+exp(1);

%f4  Griewak d=30
% z=1+1/4000*(u(1).^2+u(2).^2+u(3).^2+u(4).^2+u(5).^2+u(6).^2+u(7).^2+u(8).^2+u(9).^2+u(10).^2+u(11).^2+u(12).^2+u(13).^2+u(14).^2+u(15).^2+u(16).^2+u(17).^2+u(18).^2+u(19).^2+u(20).^2+u(21).^2+u(22).^2+u(23).^2+u(24).^2+u(25).^2+u(26).^2+u(27).^2+u(28).^2+u(29).^2+u(30).^2)-cos(u(1)/sqrt(1))*cos(u(2)/sqrt(2))*cos(u(3)/sqrt(3))*cos(u(4)/sqrt(4))*cos(u(5)/sqrt(5))*cos(u(6)/sqrt(6))*cos(u(7)/sqrt(7))*cos(u(8)/sqrt(8))*cos(u(9)/sqrt(9))*cos(u(10)/sqrt(10))*cos(u(11)/sqrt(11))*cos(u(12)/sqrt(12))*cos(u(13)/sqrt(13))*cos(u(14)/sqrt(14))*cos(u(15)/sqrt(15))*cos(u(16)/sqrt(16))*cos(u(17)/sqrt(17))*cos(u(18)/sqrt(18))*cos(u(29)/sqrt(29))*cos(u(20)/sqrt(20))*cos(u(21)/sqrt(21))*cos(u(22)/sqrt(22))*cos(u(23)/sqrt(23))*cos(u(24)/sqrt(24))*cos(u(25)/sqrt(25))*cos(u(26)/sqrt(26))*cos(u(27)/sqrt(27))*cos(u(28)/sqrt(28))*cos(u(29)/sqrt(29))*cos(u(30)/sqrt(30));

%F5 roenbrock 30D
% j1=100*(u(2)-u(1).^2).^2+(u(1)-1).^2+100*(u(3)-u(2).^2).^2+(u(2)-1).^2+100*(u(4)-u(3).^2).^2+(u(3)-1).^2+100*(u(5)-u(4).^2).^2+(u(4)-1).^2;
% j2=100*(u(6)-u(5).^2).^2+(u(5)-1).^2+100*(u(7)-u(6).^2).^2+(u(6)-1).^2+100*(u(8)-u(7).^2).^2+(u(7)-1).^2+100*(u(9)-u(8).^2).^2+(u(8)-1).^2;
% j3=100*(u(10)-u(9).^2).^2+(u(9)-1).^2+100*(u(11)-u(10).^2).^2+(u(10)-1).^2+100*(u(12)-u(11).^2).^2+(u(11)-1).^2+100*(u(13)-u(12).^2).^2+(u(12)-1).^2;
% j4=100*(u(14)-u(13).^2).^2+(u(13)-1).^2+100*(u(15)-u(14).^2).^2+(u(14)-1).^2+100*(u(16)-u(15).^2).^2+(u(15)-1).^2+100*(u(17)-u(16).^2).^2+(u(16)-1).^2;
% j5=100*(u(18)-u(17).^2).^2+(u(17)-1).^2+100*(u(19)-u(18).^2).^2+(u(18)-1).^2+100*(u(20)-u(19).^2).^2+(u(19)-1).^2+100*(u(21)-u(20).^2).^2+(u(20)-1).^2;
% j6=100*(u(22)-u(21).^2).^2+(u(21)-1).^2+100*(u(23)-u(22).^2).^2+(u(22)-1).^2+100*(u(24)-u(23).^2).^2+(u(23)-1).^2+100*(u(25)-u(24).^2).^2+(u(24)-1).^2;
% j7=100*(u(26)-u(25).^2).^2+(u(25)-1).^2+100*(u(27)-u(26).^2).^2+(u(26)-1).^2+100*(u(28)-u(27).^2).^2+(u(27)-1).^2+100*(u(29)-u(28).^2).^2+(u(28)-1).^2+100*(u(30)-u(29).^2).^2+(u(29)-1).^2;
% z=j1+j2+j3+j4+j5+j6+j7;

% f6 schaffer 2
%j1=(sin(sqrt(u(:,1).^2+u(:,2).^2))).^2-0.5;
%j2=(1+0.001.*(u(:,1).^2+u(:,2).^2)).^2;
%z=j1./j2-0.5;

% f7 schwel d=30 need
% j1=abs(u(1))+abs(u(2))+abs(u(3))+abs(u(4))+abs(u(5))+abs(u(6))+abs(u(7))+abs(u(8))+abs(u(9))+abs(u(10))+abs(u(11))+abs(u(12))+abs(u(13))+abs(u(14))+abs(u(15))+abs(u(16))+abs(u(17))+abs(u(18))+abs(u(19))+abs(u(20))+abs(u(21))+abs(u(22))+abs(u(23))+abs(u(24))+abs(u(25))+abs(u(26))+abs(u(27))+abs(u(28))+abs(u(29))+abs(u(30));
% j2=abs(u(1))*abs(u(2))*abs(u(3))*abs(u(4))*abs(u(5))*abs(u(6))*abs(u(7))*abs(u(8))*abs(u(9))*abs(u(10))*abs(u(11))*abs(u(12))*abs(u(13))*abs(u(14))*abs(u(15))*abs(u(16))*abs(u(17))*abs(u(18))*abs(u(19))*abs(u(20))*abs(u(21))*abs(u(22))*abs(u(23))*abs(u(24))*abs(u(25))*abs(u(26))*abs(u(27))*abs(u(28))*abs(u(29))*abs(u(30));
% z=j1+j2;

% F8 zakjarov  30D
%j1=u(1).^2+u(2).^2+u(3).^2+u(4).^2+u(5).^2+u(6).^2+u(7).^2+u(8).^2+u(9).^2+u(10).^2+u(11).^2+u(12).^2+u(13).^2+u(14).^2+u(15).^2+u(16).^2+u(17).^2+u(18).^2+u(19).^2+u(20).^2+u(21).^2+u(22).^2+u(23).^2+u(24).^2+u(25).^2+u(26).^2+u(27).^2+u(28).^2+u(29).^2+u(30).^2;
%j2=(1/2*(1*u(1)+2*u(2)+3*u(3)+4*u(4)+5*u(5)+6*u(6)+7*u(7)+8*u(8)+9*u(9)+10*u(10)+11*u(11)+12*u(12)+13*u(13)+14*u(14)+15*u(15)+16*u(16)+17*u(17)+18*u(18)+19*u(19)+20*u(20)+21*u(21)+22*u(22)+23*u(23)+24*u(24)+25*u(25)+26*u(26)+27*u(27)+28*u(28)+29*u(29)+30*u(30))).^2;
%j3=(1/2*(1*u(1)+2*u(2)+3*u(3)+4*u(4)+5*u(5)+6*u(6)+7*u(7)+8*u(8)+9*u(9)+10*u(10)+11*u(11)+12*u(12)+13*u(13)+14*u(14)+15*u(15)+16*u(16)+17*u(17)+18*u(18)+19*u(19)+20*u(20)+21*u(21)+22*u(22)+23*u(23)+24*u(24)+25*u(25)+26*u(26)+27*u(27)+28*u(28)+29*u(29)+30*u(30))).^4;
%z=j1+j2+j3;


%%%%% ============ end ====================================
%%%%% ============ end ====================================
%d=2
% z=u(1).^2+u(2).^2-cos(18*u(1))-cos(18*u(2));
function z=Ufun(u,a,k,m)
z=k.*((u-a).^m).*(u>a)+k.*((-u-a).^m).*(u<(-a));


    

