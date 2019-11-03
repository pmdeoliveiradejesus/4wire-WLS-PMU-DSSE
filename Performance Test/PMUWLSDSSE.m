%To compile this prograM: mcc -m PMUWLSDSSE.m -o PMUWLSDSSE
clear all
close all
clc
disp('-------------------')
disp('Advanced 4-wire PMU-WLS-DSSE')
disp('Generalized m layer 2^n node 4wire distribution test system')
disp('Paulo De Oliveira, Nestor Rodriguez, David Celeita, Gustavo Ramos')
disp('pm.deoliveiradejes@uniandes.edu.co')
disp('V.1.0 July 3, 2019')
disp('-------------------')   
%% General parameters
layer=7;%number of layers
nl=.03;%noise level .03=3%
econv=10^-4; %convergence criteria
%% To get 2-bus IEEE Paper Resuls set layer=1
%    layer=1;
%% Runs the  4-wire Power Flow for each layer
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

%% Begins the iterative process
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
lowerbound=-1;
upperbound=1;
for j=1:m
    zalt(j,1)=z0(j,1)*(1+nl*(-lowerbound+(lowerbound+upperbound)*rand(1,1))); 
        %zalt(j,1)=z0(j,1)*(1+nl*unifrnd(lowerbound,upperbound)); 
        %zalt(j,1)=z0(j,1)*(1+nl*normrnd(lowerbound,upperbound));         
end
if layer==1
% Only for l=1 and nl=3% - Paper example
   zalt=[5.13685271701070;1.71619517122123;-6.73743427864737;4.99451346817954;1.24238044981855;-6.72906146331143;4.99988680620950;-6.6900866719488;1.86374451728886;4.67983140587769;-6.90885103136665;1.92974574908040;0.430833963684793;-0.0702232038049941;-0.233552074214137;-0.100686290465729;-0.439101497362902;0.0692593051264286;0.230127146593019;0.0979942874546346;0.139171102396880;-0.513222174262098;0.263484226235493;0.107544229788664;-0.133781221303544;0.509903503818708;-0.273085400259189;-0.108938112961564;4.94990000000000;-0.100686290465729;0.0979942874546346;0.107544229788664;-0.108938112961564];
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
% dx=-inv(H'*W*H)*H'*W*(zalt-h); %Direct without Cholesky
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
%% Key performance indexes
InverseTime=mean(time)% Gain matrix factorization time
ConvergenceTime=sum(time)% Gain matrix factorization time
IterNumber=iteration%number of iteration to reach convergence
J=(zalt-h)'*W*(zalt-h) %Residual
Confidence=J-chi2inv(.01,m-ustat);
pValue =1-chi2cdf(J,m-ustat)%Confidence level >.99% not suspicious bad data
% Estimated active power losses
for k=1:4*nb;
vc(k,1)=complex(v(k),v(k+4*nb)) ;  
 vcm(k,1)=complex(vm(k),vm(k+4*nb)) ;  
end
for k=1:8*nb
ij(k,1)=h(k+6*nb);  
end
for k=1:4*nb
ijc(k,1)=complex(ij(k),ij(k+4*nb)) ;  
end
P=real(vc).*real(ijc)+imag(vc).*imag(ijc);
LossesPe=sum(P);
%% Verify i=[Y].v and i=^[Y*].v
sum(ij-Ym*v);
kk=0;
for k=8:4:4*nb
    kk=kk+1;
Ynew(k,k)=Ynew(k,k)-inv(r(kk));
end
kk=0;
for k=4*nb+8:4:8*nb
        kk=kk+1;
Ynew(k,k)=Ynew(k,k)-inv(r(kk));
end
sum(ij-Ynew*v);


