clc
clear all
close all
syms  x1 x2 [20 4]
syms x1_flat x2_flat c_flat w_flat  [1 80]
syms r [1 1]
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
% s = sum((1-(X_min.*x1_flat/11)),'all')
f = -sum(W_flat.*x1_flat,'all')- sum(W_flat.*x2_flat,'all');






% for i = 1:80
% % h8(i) =X_min(i) -numstartinggoalkeepers;
% % h9(i) = sum(X_min([i,20+i,40+i,60+i]));
% % g8(i) = X_min(20*1+i)-numdefenders;
% % g9(i) = -X_min(20*1+i)+minstartingdefenders;
% % g10(i) = X_min(20*2+i)-nummidfielders;
% % g11(i) = -X_min(20*2+i)+minstartingmidfielders;
% % g12(i) = X_min(20*3+i)-numforwards;
% % g13(i) = -X_min(20*3+i)+minstartingforwards;
% g8(i) = -X_min(i);
% g9(i) = X_min(i)-1;
% end
% g8 = -r;
% g9 = r - 1;
% g10 = 1-((sum(X_min.*x1_flat,'all'))/11) - (4/11);





g= [g1;g2;g3;g4;g5;g6;g7;gplayerlims'];
h = [h1;h2;h3;h4;h5;h6;h7];


 F =  matlabFunction(f,'vars',{[x1_flat,x2_flat]});
 G = matlabFunction(g,'vars',{[x1_flat,x2_flat]});
 H = matlabFunction(h,'vars',{[x1_flat,x2_flat]});
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
 
 
 
 
 

  x0 = [0.5*ones([80,1]); 0.5*ones([80,1]) ];
 opts = optimset('Algorithm','sqp','Display','off','MaxFunEvals',inf,'MaxIter',inf,'Display','iter');
 [x fval x1_out x2_out] = runfmincon(F,x0, nonlinfc, opts);
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
writematrix(x1_out,'sqpx1.xls')
writematrix(x2_out,'sqpx2.xls')
writematrix(r_out,'sqpxmin.xls')
%   x0 = [rand([20,4]); rand([20,4]); rand( [20,4])];
%  opts = optimset('Algorithm','active-set','Display','off','MaxFunEvals',inf,'MaxIter',inf);
%  [x fval x1_out x2_out x_min_out] = runfmincon(F,x0, nonlinfc, opts);
%  writematrix(x1_out,'asx1.xls')
% writematrix(x2_out,'asx2.xls')
% writematrix(x_min_out,'asxmin.xls')
%    x0 = [zeros([20,4]); zeros([20,4]); zeros( [20,4])];
%  opts = optimset('Algorithm','interior-point','Display','off','MaxFunEvals',inf, 'MaxIter',inf);
%  runfmincon(F,x0, nonlinfc, opts);
 







 
 
 %% Runfmincon function
 function [x fval x1_out x2_out ] = runfmincon(F,x0, nonlinfc, opts)
  x0_flat = reshape(x0,1,[])
[x fval exitflag output] = fmincon(F,x0_flat,[],[],[],[],...
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

% r_out = x(161);
% table(r_out)
 
 end
% f_new = double(subs(f,{x1_flat,reshape(x1_out,1,[])))