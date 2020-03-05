clear all; 
clc; 
save_independent_run_minvalue=zeros(30,11); 
worst_sol=zeros(1,11); 
best_sol=zeros(1,11); 
mean_sol=zeros(1,11); 
std_sol=zeros(1,11); 
for i=1:11 
     [pg1 pg2 save_independent_run_minvalue(:,i)   Convergence_curve_avg(i,:)]=BOAO(i); 
end 
index=1:11; 
 
for i=1:11 
worst_sol(1,i)=max(save_independent_run_minvalue(:,i)); 
best_sol(1,i)=min(save_independent_run_minvalue(:,i)); 
mean_sol(1,i)=mean(save_independent_run_minvalue(:,i)); 
std_sol(1,i)=std(save_independent_run_minvalue(:,i)); 
end 
data=[index;best_sol;worst_sol;mean_sol;std_sol]; 
 
xlswrite('fpa.xls',data) 
save fpa_all_data; 
