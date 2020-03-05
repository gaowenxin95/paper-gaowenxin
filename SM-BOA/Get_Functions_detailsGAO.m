

function [lb,ub,dim,fobj] = Get_Functions_details(F)


switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=500;
        
    case 'F2'
        fobj = @F2;
        lb=-10;
        ub=10;
        dim=500;
        
    case 'F3'
        fobj = @F3;
        lb=-100;
        ub=100;
        dim=500;
        
    case 'F4'
        fobj = @F4;
        lb=-100;
        ub=100;
        dim=500;
        
    case 'F5'
        fobj = @F5;
        lb=-30;
        ub=30;
        dim=30;
        
    case 'F6'
        fobj = @F6;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F7'
        fobj = @F7;
        lb=-1.28;
        ub=1.28;
        dim=500;
        
    case 'F8'
        fobj = @F8;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F9'
        fobj = @F9;
        lb=-5.12;
        ub=5.12;
        dim=500;
        
    case 'F10'
        fobj = @F10;
        lb=-32;
        ub=32;
        dim=500;
        
    case 'F11'
        fobj = @F11;
        lb=-600;
        ub=600;
        dim=500;
        
    case 'F12'
        fobj = @F12;
        lb=-50;
        ub=50;
        dim=500;
        
    case 'F13'
        fobj = @F13;
        lb=-10;
        ub=10;
        dim=2;
        
    case 'F14'
        fobj = @F14;
        lb=-65.536;
        ub=65.536;
        dim=2;
        
    case 'F15'
        fobj = @F15;
        lb=-5;
        ub=5;
        dim=4;
        
    case 'F16'
        fobj = @F16;
        lb=-2;
        ub=2;
        dim=2;
        
    case 'F17'
     fobj = @F17;
     lb = -5.12;
     ub =  5.12;
     dim = 2;
        
    case 'F18'
        fobj = @F18;
        lb=-2;
        ub=2;
        dim=2;
        
%     case 'F19'
%         fobj = @F19;
%         lb=1;
%         ub=3;
%         dim=3;
     
    case 'F19'
        fobj = @F19;
        lb=0;
        ub=1;
        dim=6;     
        
    case 'F20'
        fobj = @F20;
        lb=0;
        ub=10;
        dim=4;    
        
    case 'F21'
        fobj = @F21;
        lb=0;
        ub=10;
        dim=4;    
        
    case 'F22'
        fobj = @F22;
        lb=0;
        ub=10;
        dim=4;         
%-----新添加的Easom Funtion-----------
    case 'F23'
        fobj = @F23;
        lb = -100;
        ub = 100;
        dim = 2;
 
%-------------------------------------
        
%-----新添加的Schaffer 多峰-----------
    case 'F24'
        fobj = @F24;
        lb = -100;
        ub = 100;
        dim = 2;
 
 %-------------------------------------
 %------Alpine函数 最小值0-------------   
    case 'F25'
        fobj = @F25;
        lb = -10;
        ub = 10;
        dim = 500;
        
    case 'F26'
        fobj = @F26;
        lb = -100;
        ub = 100;
        dim = 2;
        
     case 'F27'
        fobj = @F27;
        lb = -4;
        ub = 5;
        dim = 500;
     
     case 'F28'
        fobj = @F28;
        lb = 0;
        ub = 4;
        dim = 4;
        
     case 'F29'
        fobj = @F29;
        lb = -1.28;
        ub = 1.28;
        dim =32;
        
     case 'F30'
        fobj = @F30;
        lb=-10;
        ub=10;
        dim=2;
        
     case 'F31'
        fobj = @F31;
        lb=-10;
        ub=10;
        dim=500;
        
     case 'F32'
        fobj = @F32;
        lb=-4.5;
        ub=4.5;
        dim=2;
        
     case 'F33'
        fobj = @F33;
        lb = -100;
        ub = 100;
        dim = 2;
        
     case 'F34'
        fobj = @F34;
        lb =0;
        ub = pi;
        dim = 10;  
      
      case 'F35'
        fobj = @F35;
        lb=-10;
        ub=10;
        dim=30;
        
      case 'F36'
        fobj = @F36;
        lb=-10;
        ub=10;
        dim=2;
        
      case 'F37'
        fobj = @F37;
        lb=-10;
        ub=10;
        dim=4;
        
      case 'F38'
        fobj = @F38;
        lb=-100;
        ub=100;
        dim=30;
        
      case 'F39'
        fobj = @F39;
        lb=-5;
        ub=5;
        dim=3;
       
      case 'F40'
        fobj = @F40;
        lb=-5;
        ub=10;
        dim=30;
        
      case 'F41'
        fobj = @F41;
        lb=-10;
        ub=10;
        dim=500;
        
      case 'F42'
        fobj = @F42;
        lb=-65536;
        ub=65536;
        dim=2;
       
      case 'F43'
        fobj = @F43;
        lb=0;
        ub=1;
        dim=3;
        
        
      case 'F44'
        fobj = @F44;
        lb=-5;
        ub=5;
        dim=4;
       
      case 'F45'
        fobj = @F45;
        lb=0;
        ub=10;
        dim=2;
       
      case 'F46'
        fobj = @F46;
        lb=-50;
        ub=50;
        dim=30;
        
        
        
      case 'F47'
         fobj = @F47;
         lb = -5.12;
         ub =  5.12;
         dim = 5;
     
     
      case 'F48'
         fobj = @F48;
         lb = -pi;
         ub =  pi;
         dim = 2;
         
      case 'F49'
         fobj = @F49;
         lb = -100;
         ub =  100;
         dim = 30;
   
end

end

% F1 Sphere

function o = F1(x)
o=sum(x.^2);
end

% F2 Schwefel2.22

function o = F2(x)
o=sum(abs(x))+prod(abs(x));
end

% F3 Schwefel 1.2

function o = F3(x)
dim=size(x,2);
o=0;
for i=1:dim
    o=o+sum(x(1:i))^2;
end
end

% F4 Schwefel 2.21

function o = F4(x)
o=max(abs(x));
end

% F5 Rosenbrock

function o = F5(x)
dim=size(x,2);
o=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F6 Step

function o = F6(x)
o=sum(abs((x+.5)).^2);
end

% F7 Quartic function with noise

function o = F7(x)
dim=size(x,2);
o=sum((1:dim).*(x.^4))+rand;
end

% F8 Schwefel 2.26

function o = F8(x)
o=sum(-x.*sin(sqrt(abs(x))));
end

% F9 Rastrigin

function o = F9(x)
dim=size(x,2);
o=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F10 Ackley

function o = F10(x)
dim=size(x,2);
o=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F11 Griewank

function o = F11(x)
dim=size(x,2);
o=sum(x.^2)/4000-prod(cos(x./sqrt((1:dim))))+1;
end

% F12 

function o = F12(x)
dim=size(x,2);
o=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
(1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F13 Levy

function o = F13(x)
dim=size(x,2);
o=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F14 DE JONG

function o = F14(x)
aS=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;...
-32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];

for j=1:25
    bS(j)=sum((x'-aS(:,j)).^6);
end
o=(1/500+sum(1./([1:25]+bS))).^(-1);
end

% F15

function o = F15(x)
aK=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
bK=[.25 .5 1 2 4 6 8 10 12 14 16];bK=1./bK;
o=sum((aK-((x(1).*(bK.^2+x(2).*bK))./(bK.^2+x(3).*bK+x(4)))).^2);
end

% F16 Six-hump camel back 

function o = F16(x)
o=4*(x(1)^2)-2.1*(x(1)^4)+(x(1)^6)/3+x(1)*x(2)-4*(x(2)^2)+4*(x(2)^4);
end

% F17 Branins

function o = F17(x)
o=(x(2)-(x(1)^2)*5.1/(4*(pi^2))+5/pi*x(1)-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
end


% F18 Goldstein price

function o = F18(x)
o=(1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*(x(1)^2)-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*...
    (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*(x(1)^2)+48*x(2)-36*x(1)*x(2)+27*(x(2)^2)));
end

% % F19
% 
% function o = F19(x)
% aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
% pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
% o=0;
% for i=1:4
%     o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
% end
% end

%F20

function o = F19(x)
aH=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
cH=[1 1.2 3 3.2];
pH=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
.2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F21

function o = F20(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:5
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F22

function o = F21(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:7
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F23

function o = F22(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:10
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end


%F23 Easom
function o = F23(x)
o = -cos(x(1,1))*cos(x(1,2))*exp(-(x(1,1)-pi).^2-(x(1,2)-pi).^2);
end



% F24 多峰函数
% F24
% function o = F24(x)
% o = (sin(sqrt(sum(x.^2)))^2-0.5)/(1+0.001*(sum(x.^2)))^2-0.5;
% end
%------------------------------

%F25 Apline
function o = F25(x)
dim = size(x,2);
s = 0;
for i = 1:dim
    s = s + abs(x(i)*sin(x(i)) + 0.1*x(i));
end
o = s;
end


%F26 Schaffer
function o= F26(x)
o = 0.5 + (sin(sqrt(x(:,1).^2 + x(:,2).^2)).^2 - 0.5) ./ (1 + 0.001.*(x(:,1).^2 + x(:,2).^2)).^2;
end


%F27 powell
function o = F27(x)
dim = size(x,2);
s= 0;
for j = 1:dim/4
    s = s + (x(:,4*j-3) + 10.*x(:,4*j-2)).^2 + 5.*(x(:,4*j-1) - x(:,4*j)).^2 + (x(:,4*j-2) - x(:,4*j-1)).^4 + 10.*(x(:,4*j-3) - x(:,4*j)).^4;
end
o = s;
end

%F28 power sum
function o = F28(x)
n = 4;
b = [8,18,44,114];
s = 0;
for k = 1:n
    f1 = 0;
    for j = 1:n
        f1 = f1 + x(:,j).^k;
    end
    s = s + (f1 - b(k)).^2;
end
o = s;
end

%F29 Quartic
function o = F29(x)
dim = size(x,2);
s = 0;
for j = 1:dim
    s = s + j.*x(:,j).^4;
end
o = s + rand(1);
end

%F30 Shubert
function o = F30(x)
f1 = 0;
f2 = 0;
for j = 1:5
    f1 = f1 + j.* cos((j + 1).*x(:,1)+j);
    f2 = f2 + j.* cos((j + 1).*x(:,2)+j);
end
o = f1.*f2;
end

%F31 Sum Squares
function o= F31(x)
dim = size(x,2);
s = 0;
for j = 1:dim
    s = s + j*x(:,j).^2;
end
o = s;
end

%F32 Beale
function o = F32(x)
 %%Beale - [-4.5,4.5] - Dim2;
o = (1.5 - x(:,1) + x(:,1).*x(:,2)).^2 + (2.25 - x(:,1) +x(:,1).*x(:,2).^2).^2 + (2.625 - x(:,1) + x(:,1).*x(:,2).^3).^2;
end


%F33 Bohachevsky
function o= F33(x)
 %Bohachevsky1 - [-100,100] - Dim2
o= x(:,1).^2 + 2.*x(:,2).^2 - 0.3.*cos(3*pi.*x(:,1)) - 0.4.*cos(4*pi.*x(:,2)) + 0.7;
end

%F34 Michalewicz
function o =F34(x)
%- Michalewicz - [0,pi] - Dim10
dim = size(x,2);
s = 0;
for j = 1:dim
    s = s + sin(x(:,j)).*sin(j.*x(:,j).^2 / pi).^20;
end
o= -s;
end

%F35 Matyas
function o = F35(x)
% - Matyas - [-10,10] - Dim2
o = 0.26.*(x(:,1).^2 + x(:,2).^2) - 0.48.*x(:,1).*x(:,2);
end

%F36 Booth
function o =F36(x)
% - Booth - [-10,10] - Dim2
o = (x(:,1) + 2.*x(:,2)- 7).^2 + (2.*x(:,1) + x(:,2) - 5).^2;
end

%F37  Shekel
function o = F37(x)
% - Shekel7 - [0,10] - Dim4
dim = size(x,2);

m = 7;
a = ones(10,4);
a(1,:) = 4.0*a(1,:);
a(2,:) = 1.0*a(2,:);
a(3,:) = 8.0*a(3,:);
a(4,:) = 6.0*a(4,:);
for j = 1:2
    a(5,2*j-1) = 3.0; a(5,2*j) = 7.0;
    a(6,2*j-1) = 2.0; a(6,2*j) = 9.0;
    a(7,j)     = 5.0; a(7,j+2) = 3.0;
    a(8,2*j-1) = 8.0; a(8,2*j) = 1.0;
    a(9,2*j-1) = 6.0; a(9,2*j) = 2.0;
    a(10,2*j-1)= 7.0; a(10,2*j)= 3.6;
end
c(1) = 0.1; c(2) = 0.2; c(3) = 0.2; c(4) = 0.4; c(5) = 0.4;
c(6) = 0.6; c(7) = 0.3; c(8) = 0.7; c(9) = 0.5; c(10)= 0.5;
s = 0;
for j = 1:m
    f1 = 0;
    for i = 1:dim
        f1 = f1 + (x(:,i) - a(j,i)).^2;
    end
    s = s + 1 ./ (f1 + c(j));
end
o = -s;
end

%F38 Trid6
function  o =F38(x)
% Trid6 - [-36,36] - Dim6
dim = size(x,2);
f1 = 0;
f2 = 0;
for j = 1:dim
    f1 = f1 + (x(:,j) - 1).^2;
end
for j = 2:dim
    f2 = f2 + x(:,j).*x(:,j-1);
end
o = f1 - f2;
end


%F39 Perm
function  o= F39(x)
% Perm - [-D,D] - Dim4
n = 4;
b = 0.5;
s = 0;
for k = 1:n
    f1 = 0;
    for j = 1:n
        f1 = f1 + (j^k + b).*((x(:,j)./ j).^k - 1);
    end
    s = s + f1.^2;
end
o = s;
end

% F40 Zakharov
function o = F40(x)
% - Zakharov - [-5,10] - Dim10
dim = size(x,2);
f1 = 0;
f2 = 0;
f3 = 0;
for j = 1:dim
    f1 = f1 + x(:,j).^2;
    f2 = f2 + 0.5.*j.*x(:,j);
    f3 = f3 + 0.5.*j.*x(:,j);
end
o = f1 + f2.^2 + f3.^4;
end

%F41 dixonprice
function o = F41(x)
% dixonprice - [-10,10] - Dim30
dim = size(x,2);
s = 0;
for j = 2:dim
    s = s + j.*(2.*x(:,j).^2 - x(:,j-1)).^2;
end
o = s + (x(:,1) - 1).^2;
end

%F42 FoxHoles
function o= F42(x)
%- FoxHoles - [-65536,65536] - Dim2

ai0 = [-32, -16, 0, 16, 32];
a = [repmat(ai0, 1, 5);reshape(repmat(ai0, 5 , 1), 1, 25);];

f2 = 0;
for j = 1:25
    f1 = 0;
    for k = 1:2
        f1 = f1 + (x(:,k) - a(k,j)).^6;
    end
    f2 = f2 + 1 ./ (j + f1);
end

o = 1 ./ ((1/500) + f2);


end

%F43 - Hartman3
function o = F43(x)
%F39 - Hartman3 - [0,1] - Dim3
dim = size(x,2);
a(:,2)=10.0*ones(4,1);
for j=1:2
    a(2*j-1,1)=3.0; a(2*j,1)=0.1;
    a(2*j-1,3)=30.0; a(2*j,3)=35.0;
end
c(1)=1.0;c(2)=1.2;c(3)=3.0;c(4)=3.2;
p(1,1)=0.36890;p(1,2)=0.11700;p(1,3)=0.26730;
p(2,1)=0.46990;p(2,2)=0.43870;p(2,3)=0.74700;
p(3,1)=0.10910;p(3,2)=0.87320;p(3,3)=0.55470;
p(4,1)=0.03815;p(4,2)=0.57430;p(4,3)=0.88280;
s = 0;
for i = 1:4
    f1 = 0;
    for j = 1:dim
        f1 = f1 + a(i,j).*(x(:,j) - p(i,j)).^2;
    end
    s = s + c(i)*exp(-f1);
end
o = -s;
end



%F44 Kowalik
function o = F44(x)
% - Kowalik - [-5,5] - Dim4
a = [0.1957 0.1947 0.1735 0.1600 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246];
b = 1./ [0.25 0.5 1 2 4 6 8 10 12 14 16];
s = 0;
for j = 1:11
    f1 = x(:,1).*(b(j)^2 + b(j).*x(:,2));
    f2 = b(j)^2 + b(j).*x(:,3) + x(:,4);
    s = s + (a(j) - f1 ./ f2).^2;
end
o = s;
end


function [ a, c ] = langermanAC(  )
%LANGERMANA Summary of this function goes here
%   Detailed explanation goes here

a = [9.681 0.667 4.783 9.095 3.517 9.325 6.544 0.211 5.122 2.020 0.806
9.400 2.041 3.788 7.931 2.882 2.672 3.568 1.284 7.033 7.374 0.517
8.025 9.152 5.114 7.621 4.564 4.711 2.996 6.126 0.734 4.982 1.5
2.196 0.415 5.649 6.979 9.510 9.166 6.304 6.054 9.377 1.426 0.908
8.074 8.777 3.467 1.863 6.708 6.349 4.534 0.276 7.633 1.567 0.965
7.650 5.658 0.720 2.764 3.278 5.283 7.474 6.274 1.409 8.208 0.669
1.256 3.605 8.623 6.905 0.584 8.133 6.071 6.888 4.187 5.448 0.524
8.314 2.261 4.224 1.781 4.124 0.932 8.129 8.658 1.208 5.762 0.902
0.226 8.858 1.420 0.945 1.622 4.698 6.228 9.096 0.972 7.637 0.531
7.305 2.228 1.242 5.928 9.133 1.826 4.060 5.204 8.713 8.247 0.876
0.652 7.027 0.508 4.876 8.807 4.632 5.808 6.937 3.291 7.016 0.462
2.699 3.516 5.874 4.119 4.461 7.496 8.817 0.690 6.593 9.789 0.491
8.327 3.897 2.017 9.570 9.825 1.150 1.395 3.885 6.354 0.109 0.463
2.132 7.006 7.136 2.641 1.882 5.943 7.273 7.691 2.880 0.564 0.714
4.707 5.579 4.080 0.581 9.698 8.542 8.077 8.515 9.231 4.670 0.352
8.304 7.559 8.567 0.322 7.128 8.392 1.472 8.524 2.277 7.826 0.869
8.632 4.409 4.832 5.768 7.050 6.715 1.711 4.323 4.405 4.591 0.813
4.887 9.112 0.170 8.967 9.693 9.867 7.508 7.770 8.382 6.740 0.811
2.440 6.686 4.299 1.007 7.008 1.427 9.398 8.480 9.950 1.675 0.828
6.306 8.583 6.084 1.138 4.350 3.134 7.853 6.061 7.457 2.258 0.964
0.652 2.343 1.370 0.821 1.310 1.063 0.689 8.819 8.833 9.070 0.789
5.558 1.272 5.756 9.857 2.279 2.764 1.284 1.677 1.244 1.234 0.360
3.352 7.549 9.817 9.437 8.687 4.167 2.570 6.540 0.228 0.027 0.369
8.798 0.880 2.370 0.168 1.701 3.680 1.231 2.390 2.499 0.064 0.992
1.460 8.057 1.336 7.217 7.914 3.615 9.981 9.198 5.292 1.224 0.332
0.432 8.645 8.774 0.249 8.081 7.461 4.416 0.652 4.002 4.644 0.817
0.679 2.800 5.523 3.049 2.968 7.225 6.730 4.199 9.614 9.229 0.632
4.263 1.074 7.286 5.599 8.291 5.200 9.214 8.272 4.398 4.506 0.883
9.496 4.830 3.150 8.270 5.079 1.231 5.731 9.494 1.883 9.732 0.608
4.138 2.562 2.532 9.661 5.611 5.500 6.886 2.341 9.699 6.500 0.326];
    
c = a(:,end);
a(:,end) = [];
end


% F45 - Langerman2
function o = F45(x)
%F45 - Langerman2 - [0,10] - Dim2
dim = size(x,2);
s = 0;
[a,c] = langermanAC();
for i = 1:dim
    f1 = 0;
    for j = 1:dim
        f1 = f1 + (x(:,j) - a(i,j)).^2;
    end
    s = s + c(i) .* exp((-1/pi)*f1).*cos(pi*f1);
end
o = -s;
end








% F46 - Penalized 
function o = F46(x)
%F46 - Penalized - [-50,50] - Dim30

dim = size(x,2);
for j = 1:dim
    y(:,j) = 1 + 0.25.*(x(:,j) + 1);
end
f1 = 10.*sin(pi*y(:,1)).^2;
f2 = 0;
for j = 1:dim-1
    f2 = f2 + (y(:,j) - 1).^2.*(1 + 10.*sin(pi.*y(:,j+1)).^2);
end
o = (pi/dim).*(f1 + f2 + (y(:,dim) - 1).^2) + findU(x, 10, 100, 4,size(x,1));
end


%F47 - Stepint 
function o = F47(x)
%F47 - Stepint - [-5.12 5.12] - Dim5
dim = size(x,2);
s = 0;
for j = 1:dim
    s = s + ceil(x(:,j));
end
o = s + 25;
end

% F48 - Fletcher& Powell2
function o =F48(x)
%F48 - Fletcher & Powell2 - [-pi,pi] - Dim2

s = 0;
[ a,b,alph] = FletcherABalpha();
for i = 1:dim
    f1 = 0;
    f2 = 0;
    for j = 1:dim
        f1 = f1 + a(i,j)*sin(alph(j)) + b(i,j)*cos(alph(j));
        f2 = f2 + a(i,j)*sin(x(:,j)) + b(i,j)*cos(x(:,j));
    end
    s = s + (f1 - f2).^2;
end
o = s;
end




function o=Ufun(x,a,k,m)
o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end
