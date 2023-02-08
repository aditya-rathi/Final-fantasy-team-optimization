clc
clear all
close all
syms  x1 x2 [20 4]
syms x1_flat x2_flat c_flat w_flat X_min  [1 80]
X = [x1;x2];

%% Setup LP
W = readmatrix('Exp_Score.csv');
C = readmatrix('Cost.csv');

W_flat = reshape(W,1,[]);
C_flat = reshape(C,1,[]);



% f = -sum((0.9*W_flat.*x1_flat+0.1*W_flat.*x2_flat),'all');
budget = 100;
numgoalkeepers = 2;
numdefenders = 5;
nummidfielders = 5;
numforwards = 3;
numplayers = 15;
numstarters = 11;

minstartingdefenders = 3;
minstartingmidfielders = 2;
minstartingforwards = 1;
numstartinggoalkeepers = 1;

g1 = sum(x1_flat.*C_flat + x2_flat.*C_flat)-budget;
h1 = sum(x1_flat)+sum(x2_flat)-numplayers;
h2 = sum(x1_flat, 'all')-numstarters;
h3 = sum(x1_flat(20*0+1:20*0+20))+sum(x2_flat(20*0+1:20*0+20))-numgoalkeepers;
h4 = sum(x1_flat(20*1+1:20*1+20))+sum(x2_flat(20*1+1:20*1+20))-numdefenders;
h5 = sum(x1_flat(20*2+1:20*2+20))+sum(x2_flat(20*2+1:20*2+20))-nummidfielders;
h6 = sum(x1_flat(20*3+1:20*3+20))+sum(x2_flat(20*3+1:20*3+20))-numforwards;
h7 = sum(x1_flat(20*0+1:20*0+20))-numstartinggoalkeepers;
g2 = sum(x1_flat(20*1+1:20*1+20))-numdefenders;
g3 = -sum(x1_flat(20*1+1:20*1+20))+minstartingdefenders;
g4 = sum(x1_flat(20*2+1:20*2+20))-nummidfielders;
g5 = -sum(x1_flat(20*2+1:20*2+20))+minstartingmidfielders;
g6 = sum(x1_flat(20*3+1:20*3+20))-numforwards;
g7 = -sum(x1_flat(20*3+1:20*3+20))+minstartingforwards;




for i = 1:20
gplayerlims(i) = sum(x1_flat([i,20+i,40+i,60+i])) + sum(x2_flat([i,20+i,40+i,60+i])) - 3;
end


% g= [g1;g2;g3;g4;g5;g6;g7;gplayerlims'];
% h = [h1;h2;h3;h4;h5;h6;h7];
% 
% 
%  F =  matlabFunction(f,'vars',{[x1_flat,x2_flat]});
%  G = matlabFunction(g,'vars',{[x1_flat,x2_flat]});
%  H = matlabFunction(h,'vars',{[x1_flat,x2_flat]});

%  nonlinfc= @(in1)deal(G(in1),H(in1));
 
 
 
 %% Run Different Scenarios
%  x0 = [zeros([20,4]); zeros( [20,4])];
%  opts = optimset('Algorithm','sqp','Display','off');
%  
%  runfmincon(F,x0, nonlinfc, opts);
%  
%   x0 = [zeros([20,4]); zeros( [20,4])];
%  opts = optimset('Algorithm','active-set','Display','off');
% 
%  
%  runfmincon(F,x0, nonlinfc, opts);
%  
%  
%  
%   x0 = [zeros([20,4]); zeros( [20,4])];
%  opts = optimset('Algorithm','interior-point','Display','off');
%  
%  runfmincon(F,x0, nonlinfc, opts);
%  
%   x0 = [ones([20,4]); ones( [20,4])];
%  opts = optimset('Algorithm','sqp','Display','off');
%  
%  runfmincon(F,x0, nonlinfc, opts);
%  
%   x0 = [ones([20,4]); ones( [20,4])];
%  opts = optimset('Algorithm','active-set','Display','off');
% 
%  
%  runfmincon(F,x0, nonlinfc, opts);
%  
%  
%  
%   x0 = [ones([20,4]); ones( [20,4])];
%  opts = optimset('Algorithm','interior-point','Display','off');
%  x0_flat = reshape(x0,1,[]);
%  
%  runfmincon(F,x0_flat, nonlinfc, opts);
%  
%  
%  
%  
%  
% 
%    x0 = [zeros([20,4]); zeros( [20,4])];
%  opts = optimset('Algorithm','sqp','Display','off','MaxFunEvals',inf,'MaxIter',inf);
%  runfmincon(F,x0, nonlinfc, opts);
%     x0 = [zeros([20,4]); zeros( [20,4])];
%  opts = optimset('Algorithm','active-set','Display','off','MaxFunEvals',inf,'MaxIter',inf);
%  runfmincon(F,x0, nonlinfc, opts);
%    x0 = [zeros([20,4]); zeros( [20,4])];
%  opts = optimset('Algorithm','interior-point','Display','off','MaxFunEvals',inf, 'MaxIter',inf);
%  runfmincon(F,x0, nonlinfc, opts);
 
 
 %% NLP
 s = (1-((sum(X_min.*x1_flat,'all'))/11))
f = -sum(X_min.*W_flat.*x1_flat,'all')-sum(s*X_min.*W_flat.*x2_flat,'all');


count = 1;
for i = 0.9:0.01:1
    
% h8(i) =X_min(i) -numstartinggoalkeepers;
% h9(i) = sum(X_min([i,20+i,40+i,60+i]));
% g8(i) = X_min(20*1+i)-numdefenders;
% g9(i) = -X_min(20*1+i)+minstartingdefenders;
% g10(i) = X_min(20*2+i)-nummidfielders;
% g11(i) = -X_min(20*2+i)+minstartingmidfielders;
% g12(i) = X_min(20*3+i)-numforwards;
% g13(i) = -X_min(20*3+i)+minstartingforwards;
% g8(i) = -X_min(i);
% g9(i) = X_min(i)-1;


 g10 = s-0.4;


for j=1:80
h8(j) = X_min(j)-i;
end

g= [g1;g2;g3;g4;g5;g6;g7;gplayerlims'];
h = [h1;h2;h3;h4;h5;h6;h7;h8'];


 F =  matlabFunction(f,'vars',{[x1_flat,x2_flat,X_min]});
 G = matlabFunction(g,'vars',{[x1_flat,x2_flat,X_min]});
 H = matlabFunction(h,'vars',{[x1_flat,x2_flat,X_min]});
 nonlinfc= @(in1)deal(G(in1),H(in1));





%% New scenarios
%  x0 = [zeros([20,4]);zeros([20,4]);  zeros( [20,4])];
%  opts = optimset('Algorithm','sqp','Display','off');
%  
%  runfmincon(F,x0, nonlinfc, opts);
%  
%   x0 = [zeros([20,4]); zeros([20,4]); zeros( [20,4])];
%  opts = optimset('Algorithm','active-set','Display','off');
% 
%  
%  runfmincon(F,x0, nonlinfc, opts);
%  
%  
%  
%   x0 = [zeros([20,4]); zeros([20,4]); zeros( [20,4])];
%  opts = optimset('Algorithm','interior-point','Display','off');
%  
%  runfmincon(F,x0, nonlinfc, opts);
%  
%   x0 = [ones([20,4]);ones([20,4]); ones( [20,4])];
%  opts = optimset('Algorithm','sqp','Display','off');
%  
%  runfmincon(F,x0, nonlinfc, opts);
%  
%   x0 = [ones([20,4]);ones([20,4]); ones( [20,4])];
%  opts = optimset('Algorithm','active-set','Display','off');
% 
%  
%  runfmincon(F,x0, nonlinfc, opts);
%  
%  
%  
%   x0 = [ones([20,4]);ones([20,4]); ones( [20,4])];
%  opts = optimset('Algorithm','interior-point','Display','off');
%  x0_flat = reshape(x0,1,[]);
%  
%  runfmincon(F,x0_flat, nonlinfc, opts);
 
 
 
 
 

  x0 = [ones([20,4]); ones([20,4]); ones( [20,4])];
 opts = optimset('Algorithm','sqp','Display','off','MaxFunEvals',1000,'MaxIter',1000);
  
 [x fval x1_out x2_out x_min_out lambda] = runfmincon(F,x0, nonlinfc, opts);
 iall(count) = i;
 xall{count} = x;
  fvalall(count) = fval;
  x1_outall{count} = x1_out;
  x2_outall{count} = x2_out;
  x_min_outall{count} = x_min_out;

  lambdaall{count} = [lambda.ineqnonlin; [];lambda.eqnonlin;[]; lambda.upper];
%  x0 = [x1_out; x2_out; x_min_out];
%  opts = optimset('Algorithm','sqp','Display','off','MaxFunEvals',inf,'MaxIter',inf);
%  [x fval x1_out x2_out x_min_out] = runfmincon(F,x0, nonlinfc, opts);
%   x0 = [x1_out; x2_out; x_min_out];
%  opts = optimset('Algorithm','sqp','Display','off','MaxFunEvals',inf,'MaxIter',inf);
%  [x fval x1_out x2_out x_min_out] = runfmincon(F,x0, nonlinfc, opts);
%   x0 = [x1_out; x2_out; x_min_out];
%  opts = optimset('Algorithm','sqp','Display','off','MaxFunEvals',inf,'MaxIter',inf);
%  [x fval x1_out x2_out x_min_out] = runfmincon(F,x0, nonlinfc, opts);
%   x0 = [x1_out; x2_out; x_min_out];
%  opts = optimset('Algorithm','sqp','Display','off','MaxFunEvals',inf,'MaxIter',inf);
%  [x fval x1_out x2_out x_min_out] = runfmincon(F,x0, nonlinfc, opts);
% writematrix(x1_out,'sqpx1.xls')
% writematrix(x2_out,'sqpx2.xls')
% writematrix(x_min_out,'sqpxmin.xls')
%   x0 = [rand([20,4]); rand([20,4]); rand( [20,4])];
%  opts = optimset('Algorithm','active-set','Display','off','MaxFunEvals',inf,'MaxIter',inf);
%  [x fval x1_out x2_out x_min_out] = runfmincon(F,x0, nonlinfc, opts);
%  writematrix(x1_out,'asx1.xls')
% writematrix(x2_out,'asx2.xls')
% writematrix(x_min_out,'asxmin.xls')
%    x0 = [zeros([20,4]); zeros([20,4]); zeros( [20,4])];
%  opts = optimset('Algorithm','interior-point','Display','off','MaxFunEvals',inf, 'MaxIter',inf);
%  runfmincon(F,x0, nonlinfc, opts);
 
count = count+1;



end
figure

plot(iall,-fvalall,'r.')
xlabel('x_{min} value');
ylabel('Total Expected Score');
title('Total Maximum Expected Score Possible as Expected Minutes Varies'); 
figure
hold on 
for i=1:51
plot (iall(i), sum(x1_outall{i}.*C,'all'),'b.');
% plot(iall(i), sum(x2_outall{i}.*C,'all'),'r.');
hold on 
end
xlabel('x_{min} value');
ylabel('Total spent on starting 11 (x1*C)');
title('Amount Spent on the Starting 11 Squad in Optimal Solution, as Expected Minutes Varies');


%  writematrix(x1_outall{1},'x1r5.xls')
%  writematrix(x2_outall{1},'x2r5.xls')
%  writematrix(x_min_outall{1},'xminr05.xls')
%   writematrix(lambdaall{1},'lambda05.xls')
%   writematrix(x1_outall{26},'x1r75.xls')
%  writematrix(x2_outall{26},'x2r75.xls')
%  writematrix(x_min_outall{26},'xminr75.xls')
%  writematrix(lambdaall{26},'lambda75.xls')
%    writematrix(x1_outall{41},'x1r90.xls')
%  writematrix(x2_outall{41},'x2r90.xls')
%  writematrix(x_min_outall{41},'xminr90.xls')
%   writematrix(lambdaall{41},'lambda90.xls')
%     writematrix(x1_outall{51},'x1r100.xls')
%  writematrix(x2_outall{51},'x2r100.xls')
%  writematrix(x_min_outall{51},'xminr100.xls')
%    writematrix(lambdaall{51},'lambda100.xls')
 
 %% Runfmincon function
 function [x fval x1_out x2_out x_min_out lambda] = runfmincon(F,x0, nonlinfc, opts)
  x0_flat = reshape(x0,1,[]);
[x, fval, exitflag, output, lambda, grad, hess]  = fmincon(F,x0_flat,[],[],[],[],...
    zeros(1,length(x0_flat)),[],...
    nonlinfc,opts);
disp('-----------------------------------------------------');
disp([ 'Algorithm:', opts.Algorithm]);
disp([ 'Exitflag:', num2str(exitflag)]);
disp('x0');

table(x0)
disp (['Fval = ' ,num2str(fval)]);


x1_out = x(1:80);
x2_out = x(81:160);
x1_out=reshape(x1_out,20,[]);
x1_out(x1_out<1E-5) = 0;
table(x1_out)

x2_out=reshape(x2_out,20,[]);
x2_out(x2_out<1E-5) = 0;
table(x2_out)


if length(x)>160
    x_min_out = x(161:240);
  x_min_out=reshape(x_min_out,20,[]);
x_min_out(x_min_out<1E-5) = 0;
s = 1-((sum(x_min_out.*x1_out,'all'))/11)
su = sum(x_min_out,'all')
table(x_min_out)
 end
  ineq = table(lambda.ineqnonlin,'VariableNames',{'Lambda inequality nonlinear constraints'})
eqn = table(lambda.eqnonlin,'VariableNames',{'Lambda equality nonlinear constraints'})

 end
% f_new = double(subs(f,{x1_flat,reshape(x1_out,1,[])))






