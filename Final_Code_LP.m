clc
clear all
close all
syms  x1 x2 [20 4]
syms x1_flat x2_flat c_flat w_flat [1 80]
X = [x1;x2];

%% Setup LP
W = readmatrix('Exp_Score.csv');
C = readmatrix('Cost.csv');

W_flat = reshape(W,1,[]);
C_flat = reshape(C,1,[]);



f = -sum((0.9*W_flat.*x1_flat+0.1*W_flat.*x2_flat),'all');
budget = 100;
numgoalkeepers = 2;
numdefenders = 5;
nummidfielders = 5;
numforwards = 3;
numplayers = 15;
numstarters = 11;

minstartingdefenders = 3;
minstartingmidfielders = 2;
minstartingforwards = 3;
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
gplayerlims(i) = sum(x1_flat([i,20+i,40+i,60+i]))+ sum(x2_flat([i,20+i,40+i,60+i]))-3;
end


g= [g1;g2;g3;g4;g5;g6;g7;gplayerlims'];
h = [h1;h2;h3;h4;h5;h6;h7];


 F =  matlabFunction(f,'vars',{[x1_flat,x2_flat]});
 G = matlabFunction(g,'vars',{[x1_flat,x2_flat]});
 H = matlabFunction(h,'vars',{[x1_flat,x2_flat]});

 nonlinfc= @(in1)deal(G(in1),H(in1));
 
 
 
 
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
 
 
 
 
 
 
 
 
 
 x0 = [ones([20,4]); ones( [20,4])];
 opts = optimset('Algorithm','sqp','Display','off','MaxFunEvals',inf);
 
 [x fval,exitflag x1_out x2_out] = runfmincon(F,x0, nonlinfc, opts);
 writematrix(x1_out,'sqplpx1.xls')
writematrix(x2_out,'sqplpx2.xls')
 
  x0 = [rand([20,4]); rand( [20,4])];
 opts = optimset('Algorithm','active-set','Display','off','MaxFunEvals',inf);

 
 [x fval,exitflag x1_out x2_out] = runfmincon(F,x0, nonlinfc, opts);
 writematrix(x1_out,'acx1.xls')
writematrix(x2_out,'acx2.xls')
%  
%  
%  
%   x0 = [rand([20,4]); rand( [20,4])];
%  opts = optimset('Algorithm','interior-point','Display','off','MaxFunEvals',1000000);
%  
%  runfmincon(F,x0, nonlinfc, opts);

 
 
 %% NLP
W = readmatrix('Exp_Score.csv');
 
 
 %% Runfmincon function
 function [x fval,exitflag x1_out x2_out] = runfmincon(F,x0, nonlinfc, opts)
  x0_flat = reshape(x0,1,[]);
[x fval,exitflag] = fmincon(F,x0_flat,[],[],[],[],...
    zeros(1,160),[],...
    nonlinfc,opts);
disp('-----------------------------------------------------');
disp([ 'Algorithm:', opts.Algorithm]);
disp('x0');
table(x0)
disp (['Fval = ' ,num2str(fval)]);
exitflag

x1_out = x(1:80);
x2_out = x(81:160);
x1_out=reshape(x1_out,20,[]);
x1_out(x1_out<1E-5) = 0;
table(x1_out)

x2_out=reshape(x2_out,20,[]);
x2_out(x2_out<1E-5) = 0;
table(x2_out)

 end







