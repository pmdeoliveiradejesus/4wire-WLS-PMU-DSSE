% PMU-based system state estimation for multigrounded distribution systems
%  mcc -m NewDSSE.m -o DSSE4wireTime
clear all
close all
clc
disp('-------------------')
disp('Advanced 4-wire DSSE')
disp('Generalized m layer 2^n node 4wire distribution test system')
disp('Performance Test')
disp('Paulo De Oliveira, Nestor Rodriguez, David Celeita, Gustavo Ramos')
disp('pm.deoliveiradejes@uniandes.edu.co')
disp('V.1.0 July 3, 2019')
disp('-------------------')
prompt='What is the number of layers (1-7)[default=5]? ';
l=input(prompt);
    if isempty(l)
     l = 5;
    end
prompt2='What is the size of the sample (1-1000) [default=250]? ';
Ni=input(prompt2);
    if isempty(Ni)
     Ni = 250;
    end
tic    
%% General parameters
econv=10^-4;
%l=4;%number of layers
Nei=1;% minimum noise level
Nef=4;% maximum noise level
%Ni=25;% number of iterations
%% Runs the  4-wire Power Flow for each layer
for layer=1:l
fprintf('Running 4-wire l-layer Kersting NEV Power Flow for layer %d.\n',layer);      
[zm,LossesPm,Ym,Ynewm,Ynew2m,rm,injectm,vm,SLm]=KerstingGeneric_powerflow(layer);
zx{layer}=num2cell(zm);
LossesPmx{layer}=num2cell(LossesPm);
Ymx{layer}=num2cell(Ym);
Ynewx{layer}=num2cell(Ynewm);
Ynew2x{layer}=num2cell(Ynew2m);
rx{layer}=num2cell(rm);
injectmx{layer}=num2cell(injectm);
vmx{layer}=num2cell(vm);
SLmx{layer}=num2cell(SLm);
end
%% Begins the iterative process
for layer=1:l
for noisel=Nei:Nef  
for it=1:Ni
fprintf('layer/error_type/iter %d %d %d.\n',layer,noisel,it);   
nb=2^layer;
z0=cell2mat(zx{layer});
LossesP=cell2mat(LossesPmx{layer});
Y=cell2mat(Ymx{layer});
Ynew=cell2mat(Ynewx{layer});
Ynew2=cell2mat(Ynew2x{layer});
r=cell2mat(rx{layer});
inject=cell2mat(injectmx{layer});
v=cell2mat(vmx{layer});
SL=cell2mat(SLmx{layer});
%% 
ustat=length(Y)+nb-1;% State Vars
m=length(z0);% number of measurements
%% Noise generator, altering the OpenDSS solution zm
nl=3*noisel/100;%noise level from 3 to 12%
lowerbound=-1;
upperbound=1;
for j=1:m
        zalt(j,1)=z0(j,1)*(1+nl*unifrnd(lowerbound,upperbound)); 
       %zalt(j,1)=z0(j,1)*(1+nl*normrnd(lowerbound,upperbound));         
end
kk=0;
for k=4:4:8*nb
    kk=kk+1;
zalt(m+kk)=zalt(6*nb+k);
end
%% Meter data accuracy and weights calculation
sigma0=.03;%acuraccy of the meters
SIGMAv0=sigma0*10; %Accuracy (dev stad 1*sigma0% on a scale 10000V)
SIGMAi0=sigma0*.400; %Accuracy (dev stad 1*sigma0% on a scale 400A)
SIGMAin0=sigma0*.100; %Accuracy (dev stad 1*sigma0% on a scale 100A)
SIGMAz0=3*sigma0*5; %Accuracy (dev stad 3*sigma0% on a scale 5 ohm)
SigV=ones(1,6*(nb))*SIGMAv0^-2;
SigI=[];
for i=1:nb 
SigI=horzcat(ones(1,3)*SIGMAi0^-2,ones(1,1)*SIGMAin0^-2, ones(1,3)*SIGMAi0^-2,ones(1,1)*SIGMAin0^-2,SigI);
end
SigZ=ones(1,nb-1)*SIGMAz0^-2;
We=horzcat(SigV,SigI,SigZ);
for k=1:2*nb
We(m+k)=SIGMAin0^-2;
end
W=diag(We);
WL=horzcat(SigV,SigI,SigZ)';
%% Build h(x) and H(x) matrix (Outside de while loop, only constant parameters
Ox=zeros(3,4);
A=[1 0 0 -1;0 1 0 -1;0 0 1 -1];
for k1=1:2*nb
    for k2=1:2*nb
   M{k1,k2}=Ox; 
    end
end   
for k=1:2*nb
   M{k,k}=A; 
end
M2=cell2mat(M);
[m2c,m2w]=size(M2);
%Jacobian precalculation
H=zeros(m+2*nb,ustat);
Yone=zeros(nb-1,nb*8);
H=vertcat(M2,Ynew,Yone,Ynew2);
jj=0;
for k=m-nb+2:m
  jj=jj+1;
  H(k,length(v)+jj)=1;
 end
e=1;iteration=1;
%% 4-Wire Distrinution System State Estimation Procedure
while e>econv   
% calculate h(x) Non linear elements
h1=M2*v;
% kk1=0;
h2=Ynew*v;
kk=0;
for  k=8:4:4*nb
kk=kk+1;
h2(k)=h2(k)-inv(r(kk))*v(k,1);
h2(k+4*nb)=h2(k+4*nb)-inv(r(kk))*v(k+4*nb,1);
end
h3=eye(nb-1)*r;
h4=Ynew2*v;
h=[h1',h2',h3',h4']';
%% Update the Jacobian with non-linear elements
j1=0;
 for k=m2c+8:4:m2c+nb*4
    j1=j1+1;  
    H(k,4+4*j1)=H(k,4+4*j1)-inv(r(j1));
 end
 j2=0;
 for k=m2c+4*nb+8:4:m2c+nb*8
  j2=j2+1;    
   H(k,4+4*nb+4*j2)=H(k,4+4*nb+4*j2)-inv(r(j2));
 end
  j3=0;
 for k3=m2c+8:4:m2c+nb*4
  j3=j3+1;    
   H(k3,8*nb+j3)=v(k3-m2c,1)*(r(j3))^(-2);
  end 
  j4=0;
 for k4=m2c+4*nb+8:4:m2c+nb*8
  j4=j4+1;    
   H(k4,8*nb+j4)=v(k4-m2c,1)*(r(j4))^(-2);
  end
%% Two matrices!
 Bs=H'*diag(We);
 Gs=Bs*H;
 isposdef = all(eig(Gs)) > 0;%is Gs positive definite?
 %% Inverse as a function of Bs and Gs
  time0=cputime;
R = chol(Gs)';%Choelsky decomposition L'*L=Gs
t=Bs*(h-zalt);
u(1)=0;%Forward substitution
flag=0;
 for i=1:ustat  
 u(i)=inv(R(i,i))*(t(i)-flag);
 flag=0;
 if i<ustat
 for j=1:i   
 flag=R(i+1,j)*u(j)+flag; 
 end
 end
 end
 dx(ustat,1)=0;%backward substitution
 flag=0;
 for i=ustat:-1.0:1
 dx(i,1)=inv(R(i,i))*(u(i)-flag);
 flag=0;
 if i>1
  for j=ustat:-1:i 
  flag=R(j,i-1)*dx(j,1)+flag; 
  end
 end
 end
  time(iteration)=cputime-time0;
% time3=cputime;
% dx=-inv(H'*W*H)*H'*W*(zalt-h); Direct without Cholesky
% time4(iter)=cputime-time3
x=vertcat(v,r);
x=x-dx;
for k=1:8*nb
v(k,1)=x(k);
end
for k=1:nb-1
r(k,1)=x(k+8*nb);
end
iteration=iteration+1;
e=max(abs(dx));
%% End of the DSSE
end 

% Estimated active power losses
for k=1:4*nb;
vc(k,1)=complex(v(k),v(k+4*nb)) ;  
end
for k=1:8*nb
ij(k,1)=h(k+6*nb);  
end
for k=1:4*nb
ijc(k,1)=complex(ij(k),ij(k+4*nb)) ;  
end
P=real(vc).*real(ijc)+imag(vc).*imag(ijc);
LossesPe(noisel,layer,it)=sum(P);
InverseTime(noisel,layer,it)=mean(time);
IterTime(noisel,layer,it)=iteration;
J(noisel,layer,it)=(zalt-h)'*W*(zalt-h);
% if layer==1;
%     if noisel==1; 
% Jhistogram(it)=(zalt-h)'*W*(zalt-h);
%     end
% end

%Confidence(noisel,layer,it)=J(noisel,layer,it)-chi2inv(.01,m-ustat);
pValue(noisel,layer,it)=1-chi2cdf(J(noisel,layer,it),m-ustat);%Confidence level >.99% not suspicious bad data
nn(layer)=ustat;
%[layer,noisel,it]
end
end
LssP(layer)=LossesP*1000;%kW
end

 for noisel=Nei:Nef
     for layer=1:l
meanLossesPe(noisel,layer)=mean(LossesPe(noisel,layer,:)); 
%meanLosseseP(noisel,layer)=mean(LosseseP(noisel,layer,:));
meanInvT(noisel,layer)=mean(InverseTime(noisel,layer,:));
meanIterT(noisel,layer)=mean(IterTime(noisel,layer,:));
meanResidual(noisel,layer)=mean(J(noisel,layer,:));
%meanConf(noisel,layer)=mean(Confidence(noisel,layer,:));
meanpValue(noisel,layer)=mean(pValue(noisel,layer,:));
meanConvTime(noisel,layer)=meanInvT(noisel,layer)*meanIterT(noisel,layer);
     end
 end
save Statvars meanInvT meanIterT meanResidual meanLossesPe meanConvTime LssP
toc
%% Display Output 
figure('Color','w','units','normalized','outerposition',[0 0 1 1],'name','Inverse Time',...
    'numbertitle','off');
plot(nn,1000*meanInvT(1,:),'-r',nn,1000*meanInvT(2,:),'-b',nn,1000*meanInvT(3,:),'-k',nn,1000*meanInvT(4,:),'-g')
set(gca,'FontSize',24)
tit = title('4-wire DSSE Jacobian matrix factorization average time', 'FontSize', 24);
set(tit,'Interpreter','latex');
xaxis = xlabel({'Number of state variables'}, 'FontSize', 24);
set(xaxis,'Interpreter','latex');
yaxis = ylabel({'Average factorization time (ms)'}, 'FontSize', 24);
set(yaxis,'Interpreter','latex');
leg = legend({'mean time for error ${\epsilon}=3 \%$';'mean time for error ${\epsilon}=6 \%$';'mean time for error ${\epsilon}=9 \%$'...
    ;'mean time for error ${\epsilon}=12 \%$'}, 'FontSize', 24, 'Location','southeast');
set(leg,'Interpreter','latex');
legend boxoff 

figure('Color','w','units','normalized','outerposition',[0 0 1 1],'name','Number of iterations',...
    'numbertitle','off');
plot(nn,meanIterT(1,:),'-r',nn,meanIterT(2,:),'-b',nn,meanIterT(3,:),'-k',nn,meanIterT(4,:),'-g')
set(gca,'FontSize',24)
tit = title('4-wire DSSE Number of iterations', 'FontSize', 24);
set(tit,'Interpreter','latex');
xaxis = xlabel({'Number of state variables'}, 'FontSize', 24);
set(xaxis,'Interpreter','latex');
yaxis = ylabel({'Number of iterations to converge e=10-4'}, 'FontSize', 24);
set(yaxis,'Interpreter','latex');
leg = legend({'Number of iterations for error ${\epsilon}=3 \%$';'Number of iterations for error ${\epsilon}=6 \%$';'Number of iterations for error ${\epsilon}=9 \%$'...
    ;'Number of iterations for error ${\epsilon}=12 \%$'}, 'FontSize', 24, 'Location','southeast');
set(leg,'Interpreter','latex');

figure('Color','w','units','normalized','outerposition',[0 0 1 1],'name','meanConvTime',...
    'numbertitle','off');
plot(nn,meanConvTime(1,:),'-r',nn,meanConvTime(2,:),'-b',nn,meanConvTime(3,:),'-k',nn,meanConvTime(4,:),'-g')
set(gca,'FontSize',24,'YScale', 'log')
tit = title('meanConvTime', 'FontSize', 24);
set(tit,'Interpreter','latex');
xaxis = xlabel({'Number of state variables'}, 'FontSize', 24);
set(xaxis,'Interpreter','latex');
yaxis = ylabel({'meanConvTime (seconds)'}, 'FontSize', 24);
set(yaxis,'Interpreter','latex');
leg = legend({'mean time for error ${\epsilon}=3 \%$';'mean time for error ${\epsilon}=6 \%$';'mean time for error ${\epsilon}=9 \%$'...
    ;'mean time for error ${\epsilon}=12 \%$'}, 'FontSize', 24, 'Location','southeast');
set(leg,'Interpreter','latex');
legend boxoff 



figure('Color','w','units','normalized','outerposition',[0 0 1 1],'name','LossesP',...
    'numbertitle','off');
plot(nn,meanLossesPe(1,:),'-r',nn,meanLossesPe(2,:),'-b',nn,meanLossesPe(3,:),'-k',nn,meanLossesPe(4,:),'-g',nn,LssP(:),'.b')
set(gca,'FontSize',24)
tit = title('LossesP', 'FontSize', 24);
set(tit,'Interpreter','latex');
xaxis = xlabel({'Number of state variables'}, 'FontSize', 24);
set(xaxis,'Interpreter','latex');
yaxis = ylabel({'Avg. ctive Losses (kW)'}, 'FontSize', 24);
set(yaxis,'Interpreter','latex');
leg = legend({'mean time for error ${\epsilon}=3 \%$';'mean time for error ${\epsilon}=6 \%$';'mean time for error ${\epsilon}=9 \%$'...
    ;'mean time for error ${\epsilon}=12 \%$';'Actual Active Losses',}, 'FontSize', 24, 'Location','northeast');
set(leg,'Interpreter','latex');
legend boxoff 

% figure
% hist(Jhistogram,100);
% meanInvT
