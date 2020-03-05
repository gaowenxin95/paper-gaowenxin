%_____________________________________________________________________________________________ %
%  Butterfly Optimization Algorithm (BOA) source codes demo V1.0                               %
%                                                                                              %
%  Author and programmer: Sankalap Arora                                                       %
%                                                                                              %
%         e-Mail: sankalap.arora@gmail.com                                                     %
%                                                                                              %
%  Main paper: Sankalap Arora, Satvir Singh                                                    %
%              Butterfly optimization algorithm: a novel approach for global optimization	   %
%              Soft Computing, in press,                                                       %
%              DOI: https://doi.org/10.1007/s00500-018-3102-4                                  %
%___________________________________________________________________________________________   %
%
function [fmin,best_pos,Convergence_curve]=BOA(n,N_iter,Lb,Ub,dim,fobj)

% n is the population size
% N_iter represnets total number of iterations
p=0.8;                       % probabibility switch
power_exponent=0.1;%a
sensory_modality=0.01;%c
pause=zeros(n,1);
max_pause=60;
%Initialize the positions of search agents
Sol=initialization(n,dim,Ub,Lb);

for i=1:n,
    Fitness(i)=fobj(Sol(i,:));
end

% Find the current best_pos
[fmin,I]=min(Fitness);
best_pos=Sol(I,:);
S=Sol; 

% Start the iterations -- Butterfly Optimization Algorithm 
for t=1:N_iter
  
        for i=1:n, % Loop over all butterflies/solutions
         
          %Calculate fragrance of each butterfly which is correlated with
          %objective function  Eq1
          Fnew=fobj(S(i,:));
          FP=(sensory_modality*(Fnew^power_exponent));   
    
          %Global or local search
%           L = Levy(dim);
          if rand<p,   
%               dis = rand * rand * best_pos - Sol(i,:);        %Eq. (2) in paper
             dis = rand * rand * best_pos - Sol(i,:);
             S(i,:)=Sol(i,:)+dis*FP;
             if pause(i,1)==max_pause, 
                 pause(i,:)=0; 
             end   
           else
              % Find random butterflies in the neighbourhood
              epsilon=rand;
              JK=randperm(n);
              dis=epsilon*epsilon*Sol(JK(1),:)-Sol(JK(2),:);
              S(i,:)=Sol(i,:)+dis*FP;                         %Eq. (3) in paper
        
          end
       
           % SCA optimize BOA
          for i = 1:n
            a = 2;
     
            r1=a-t*((a)/N_iter);
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();
            if r4<0.5
                
                nS(i,:)= S(i,:)+(r1*sin(r2)*abs(r3*best_pos-S(i,:)));
            else
                
                nS(i,:)= S(i,:)+(r1*cos(r2)*abs(r3*best_pos-S(i,:)));
            end
            
          end
          
          
            % Check if the simple limits/bounds are OK
            S(i,:)=simplebounds(nS(i,:),Lb,Ub);
          
            % Evaluate new solutions
            Fnew=fobj(nS(i,:));  %Fnew represents new fitness values
            
            % If fitness improves (better solutions found), update then
            if (Fnew<=Fitness(i)),
                Sol(i,:)=nS(i,:);
                Fitness(i)=Fnew;
                pause(i,1)=0;
            end
           
           % Update the current global best_pos
           if Fnew<=fmin,
                best_pos=nS(i,:);
                fmin=Fnew;
                
           else
                pause(i,1)=pause(i,1)+1; 
           end
         end
            
         Convergence_curve(t,1)=fmin;
         
         %Update sensory_modality
          sensory_modality=sensory_modality_NEW(sensory_modality, N_iter);
end

% Boundary constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb;
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub;
  % Update this new move 
  s=ns_tmp;

  
function y=sensory_modality_NEW(x,Ngen)
y=x+(0.025/(x*Ngen));



