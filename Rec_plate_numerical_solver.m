%%% Numerical solver for rectangular plate problem adressed by Jiong Wang (2018) doi: 10.1016/j.jmps.2017.10.017

%%% Output represents the variation of thickness $h1$ and growth factor $\lambda$ (L here) 

%%%% Note you must include fminsearchbnd.m in the directory otherwise code will not run 

clc;
clear all;
close all
format long
% load('more_data_L-h.mat')     % Output of the numerical solver you can plot the variation (plot(h1,L)) by loading this data  

tic
global h

up=1.01:0.001:1.16;
h1=0.015:0.01:0.2;    % thickness of plate
LL=zeros(length(up),length(h1));
fvalmat=zeros(length(up),length(h1));
ff=zeros(length(up), 2);
L=zeros(length(h1),1);
f=zeros(length(h1),1);

for  jj=1:length(h1)
    h=h1(jj);

for i=1:length(up)
    
    options = optimset('display', 'on','MaxIter',1000,'TolFun',1e-20,'TolX', 1e-20);

    [ll,fval]=fminsearchbnd(@(l)rec_plate_objective_func(l),1, 1, up(i),options);

    ff(i,1)=ll;                          % minimum value of growth factor for particular value of h1
    ff(i,2)=fval;                        % minimised objective function value.
end

    LL(:,jj)=ff(:,1);                    

    fvalmat(:,jj)=ff(:,2);               

    locsfval=find(fvalmat(:,jj)<1e-15);  % Points where we get first minimised objective fuction value for paticular h1

    L(jj)=LL(locsfval(1),jj);            % All first minimised values of $\lambda$ corresponds to h1 value until h1(end)

    f(jj)=wang_objective(L(jj));         % Objective function value corresponds to $\lambda$ (L here)
end

plot(h1,L,'*b','MarkerSize',10)
xlabel('h','FontSize',16,'FontWeight','bold');
ylabel('\lambda','FontSize',16,'FontWeight','bold');
set(gca,'FontSize',16)
ylim([1 1.09])
xlim([0 h1(end)])

% % save('more_data_L-h.mat','L','f')  %%% Output file of data for L corresponds to different h  
% % savefig('L vs h.fig')

toc