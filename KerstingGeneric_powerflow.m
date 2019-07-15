function [z0,LossesP,Y,Ynew,Ynew2,r,inject,v,SL] = run_powerflow(l)
% Test case; Kersting NEV
% Kersting, W.H. A three-phase unbalanced line model with grounded neutrals
% through a resistance.  In Proceedings of the  2008 IEEE Power and Energy 
% Society General Meeting-PESGM, Pittsburgh, PA, USA, 20--24 July 2008; 
% pp.  12651-12652.
% Two nodes, distance 1.13miles
nb=2^l;
%% Kesting NEV Database - Generic n^l bus power flow
a=complex(cos(2*pi/3),sin(2*pi/3));
a2=a^2;
Vo=12.47/sqrt(3);%Nominal voltage kV
Sbase=3;%MVA Base at High Voltage
Vbase=12.47/sqrt(3);%kV base
Ibase=1000*Sbase/Vbase;%Amperes
Zbase=Vbase^2/Sbase;%Zbase high in ohms
f=60;%Frequency Hz
rvd=100; % Resitivity ohm-m
eta=1.6093;%Impedances are Given in ohm/mile
mu0=4*pi*eta/10000;%H/mile
w=2*pi*f;%angular frequency
De=2160*sqrt(rvd/f); %Carson's correction factor 
re=(pi/4)*4*eta*pi*f*0.0001; %Carson's correction ground loop resistance (ohm)
%% DL Data
Length=6000*0.000189393939/(l);%section length in miles
%Length=0.1*60000*0.000189393939/(nb-1);%section length in miles
resist1=0.5;%ohms
resist2=5.0;%ohms
GMRf=0.0244;%feet phase
rf=1*0.306;%ohm/mile
GMRn=0.00814;%feet
rn=1*0.5920;%ohm/mile neutral
%Spacings
Dab=2.5;%feet
Dbc=4.5;%feet
Dac=Dab+Dbc;%feet
Dcn=(4*4+3*3)^.5;%feet
Dbn=(4*4+1.5*1.5)^.5;%feet
Dan=(4*4+4*4)^.5;%feet
hqa=29;%feet
hqb=29;%feet
hqc=29;%feet
for i=1:4:4
    zp(i,i)=rf+w*mu0/8+sqrt(-1)*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
    zp(i+1,i+1)=rf+w*mu0/8+sqrt(-1)*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
    zp(i+2,i+2)=rf+w*mu0/8+sqrt(-1)*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
    zp(i+3,i+3)=(rn+re+sqrt(-1)*(mu0*f*(log(2160*sqrt(rvd/f)*inv(GMRn)))));%ohm/mile
    zp(i,i+1)=w*mu0/8+sqrt(-1)*((mu0*w/(2*pi))*(log(De/(Dab))));%ohm/mile
    zp(i,i+2)=w*mu0/8+sqrt(-1)*((mu0*w/(2*pi))*(log(De/(Dac))));%ohm/mile
    zp(i,i+3)=w*mu0/8+sqrt(-1)*((mu0*w/(2*pi))*(log(De/(Dan))));%ohm/mile
    zp(i+1,i+2)=w*mu0/8+sqrt(-1)*((mu0*w/(2*pi))*(log(De/(Dbc))));%ohm/mile
    zp(i+1,i+3)=w*mu0/8+sqrt(-1)*((mu0*w/(2*pi))*(log(De/(Dbn))));%ohm/mile
    zp(i+2,i+3)=w*mu0/8+sqrt(-1)*((mu0*w/(2*pi))*(log(De/(Dcn))));%ohm/mile
    zp(i+1,i)=zp(i,i+1);%ohm/mile
    zp(i+2,i)=zp(i,i+2);%ohm/mile
    zp(i+3,i)=zp(i,i+3);%ohm/mile
    zp(i+2,i+1)=zp(i+1,i+2);%ohm/mile
    zp(i+3,i+1)=zp(i+1,i+3);%ohm/mile
    zp(i+3,i+2)=zp(i+2,i+3);%ohm/mile
end
yshunt=zeros(l*4);
% Line susceptance in siemens/mile taken from Sunderman paper
% yshunt=0*eta*[3.52623E-6 -1.14173E-6 -4.37331E-7 -5.25407E-7 ;
% -1.14173E-6 3.71665E-6  -7.26846E-7  -6.68079E-7 ;
% -4.37331E-7 -7.26846E-7 3.35207E-6  -6.71857E-7 ;
% -5.25407E-7 -6.68079E-7 -6.71857E-7 3.32103E-6 ;]*sqrt(-1);
Dbc=4.5;%feet
Dac=Dab+Dbc;%feet
Dcn=(4*4+3*3)^.5;%feet
Dbn=(4*4+1.5*1.5)^.5;%feet
Dan=(4*4+4*4)^.5;%feet
hqa=29;%feet
hqb=29;%feet
hqc=29;%feet
%% Load Demands
fpa=0.9; Sa=3.0;%power factor, MVA
fpb=0.95;Sb=3.5;%power factor, MVA
fpc=0.85;Sc=2.5;%power factor, MVA
SL(1)=-((Sa*fpa+sqrt(-1)*Sa*(1-fpa^2)^.5)  )/(nb  -1);% MVA
SL(1+1)=-((Sb*fpb+sqrt(-1)*Sb*(1-fpb^2)^.5) ) /(nb  -1);% MVA
SL(1+2)=-((Sc*fpc+sqrt(-1)*Sc*(1-fpc^2)^.5)  )/(nb -1);% MVA
SL(1+3)=0;
MVAsc3=500; %MVA aprox Icc=21kA
MVAsc1=500;
alpha=3;
beta=3;
%%% 
%% OpenDSS Load Flow
Z1m=(Vbase*sqrt(3))^2/MVAsc3;
R1=Z1m/sqrt(1+alpha^2);
X1=alpha*R1;
Z1=complex(R1,X1);
gamma=3*(Vbase*sqrt(3))^2/MVAsc1;
a=beta^2+1;
b=4*R1+4*beta*X1;
c=4*R1^2+4*X1^2-gamma^2;
R0=(-b+sqrt(b^2-4*a*c))/(2*a);
R02=(-b-sqrt(b^2-4*a*c))/(2*a);
X0=beta*R0;
Z0=complex(R0,X0);
a=-0.5+j*sqrt(3)*.5;
As=[1 1 1;1 a^2 a; 1 a a^2];
Zg=(As)*diag([Z0;Z1;Z1])*inv(As);
anglev=45;
for i=1:4:4
Vof(i,1)=pol2rec(Vo,anglev);%(Vo);% kVolts
Vof(i+1,1)=pol2rec(Vo,-120+anglev);%(Vo*a2);% kVolts
Vof(i+2,1)=pol2rec(Vo,120+anglev);%(Vo*a);% kVolts
Vof(i+3,1)=0;
end
%% Ybus generation
Yg=inv(Zg);% siemens
Yg(4,1)=-Yg(1,1)-Yg(1,2)-Yg(1,3);% siemens
Yg(4,2)=-Yg(2,2)-Yg(2,1)-Yg(2,3);% siemens
Yg(4,3)=-Yg(3,3)-Yg(1,3)-Yg(2,3);% siemens
Yg(1,4)=-Yg(1,1)-Yg(1,2)-Yg(1,3);% siemens
Yg(2,4)=-Yg(2,2)-Yg(2,1)-Yg(2,3);% siemens
Yg(3,4)=-Yg(3,3)-Yg(3,1)-Yg(3,2);% siemens
Yg(4,4)=inv(resist1)-Yg(4,1)-Yg(4,2)-Yg(4,3);% siemens
Zp=zp*Length;
Yp=inv(Zp) ;% siemens
%% Calculating all in per unit
V=[];
for i=1:nb 
V=vertcat(Vof,V);
end
 %% Building the augmented Ybus
%Yload=conj(S2)./abs(1)^2/3;% Define an addmittance load, (arbitrary)
%Connetcted to ground though resist2
Yloadx=(-conj(diag(SL))./abs(Vbase)^2/3) ;% siemens 
Yloadx(4,1)=-Yloadx(1,1);% siemens
Yloadx(4,2)=-Yloadx(2,2);% siemens
Yloadx(4,3)=-Yloadx(3,3);% siemens
Yloadx(1,4)=-Yloadx(1,1);% siemens
Yloadx(2,4)=-Yloadx(2,2);% siemens
Yloadx(3,4)=-Yloadx(3,3);% siemens
Yloadx(4,4)=inv(resist2)-Yloadx(4,1)-Yloadx(4,2)-Yloadx(4,3);% siemens 
O=zeros(4,4);
for k1=1:nb-1
    for k2=1:nb-1
   WW{k1,k2}=O; 
    end
end
        
for k=1:nb-1
   WW{k,k}=+Yloadx; 
end
WW2=cell2mat(WW);
 
for k1=1:nb
    for k2=1:nb
   Yopendss{k1,k2}=O; 
    end
end
        
for k=1:nb
   Yopendss{k,k}=Yp; 
end
Yopendss{1,1}=Yopendss{1,1}+Yg;
for k=2:nb
  Yopendss{k,k}=Yopendss{k,k}+Yloadx; 
end        
 for k=2:nb-nb/2
   Yopendss{k,k}=Yopendss{k,k}+2*Yp; 
 end   
 
 Index1(1)=1;
 for jx=2:nb
 Index2(jx-1)=jx;
 end
 jyy=2;
 for jy=2:2:nb-1
 Index1(jy)=jyy;
 Index1(jy+1)=jyy;
 jyy=jyy+1;
 end   
 for k=1:nb-1
 Yopendss{Index1(k),Index2(k)}=-Yp; 
 Yopendss{Index2(k),Index1(k)}=-Yp; 
 end
 Yopendss2=cell2mat(Yopendss);     
for k1=1:nb
    for k2=1:nb
   Ynormal{k1,k2}=O; 
    end
end      
for k=1:nb
   Ynormal{k,k}=Yp; 
end
 for k=2:nb-nb/2
   Ynormal{k,k}=Ynormal{k,k}+2*Yp; 
 end   
 Index1(1)=1;
 for jx=2:nb
 Index2(jx-1)=jx;
 end
 jyy=2;
 for jy=2:2:nb-1
 Index1(jy)=jyy;
 Index1(jy+1)=jyy;
 jyy=jyy+1;
 end   
 for k=1:nb-1
 Ynormal{Index1(k),Index2(k)}=-Yp; 
 Ynormal{Index2(k),Index1(k)}=-Yp; 
 end
 Y=cell2mat(Ynormal);
econv=10^-12;e=1;iter=1;
Zimplicit=inv(Yopendss2);
%% Load Flow Process
while e>econv
    Ic1=Yg*Vof;
    for i=1:4:4*(nb-1)
        I2(i,1)=conj(SL(1)/(V(i+4)-V(i+7)));% compensation injected currents at bus loads
        I2(i+1,1)=conj(SL(2)/(V(i+5)-V(i+7)));% compensation injected currents at bus loads
        I2(i+2,1)=conj(SL(3)/(V(i+6)-V(i+7)));% compensation injected currents at bus loads
        I2(i+3,1)=-I2(i,1)-I2(i+1,1)-I2(i+2,1)-V(i+7)/(resist2);
    end 
    for k=1:4*(nb-1)
        V2(k,1)=V(k+4);
    end
    I2c= I2 + WW2*V2 ;
    Ic=vertcat(Ic1,I2c);
    Vaux=Zimplicit*Ic;
    e=max(abs(V-Vaux));
    ex=sum(abs((Ic-Yopendss2*V)));
    V=Vaux;
    iter=iter+1;
end
%OpenDSS Results
%% End of OpenDSS callulation 
Va=[V(1);V(2);V(3);V(4);];
Vb=[V(5);V(6);V(7);V(8);];
I1=Yp*Va-Yp*Vb;
I=vertcat(I1,I2);
%sum((I-Y*V)); %verification
P=real(V).*real(I)+imag(V).*imag(I);
LossesP=sum(P);
G=real(Y);B=imag(Y);
YY=vertcat(horzcat(G,-B),horzcat(B,G));%pu
II=vertcat(real(I),imag(I));%pu
VV=vertcat(real(V),imag(V));%pu
%sum((II-YY*VV));
m=0;
for i=1:4:8*(nb)
    zm(i-m,1)=VV(i)-VV(i+3);
    zm(i+1-m,1)=VV(i+1)-VV(i+3);
    zm(i+2-m,1)=VV(i+2)-VV(i+3);
    m=m+1;
end
for k=1:8*nb
    zm(end+1,1)=II(k);
end
for i=1:nb -1
    zm(end+1)=resist2;
end
n=length(zm);
for k=1:nb-1
r(k,1)=resist2;
end
Y=YY;
Ynew=YY; 
for k=4:4:8*nb
Ynew(k,:)=-Ynew(k-3,:)- Ynew(k-2,:)-Ynew(k-1,:);
end
Ynew(4,4)=Ynew(4,4)-inv(resist1);
Ynew(4+4*nb,4+4*nb)=Ynew(4+4*nb,4+4*nb)-inv(resist1);
z0=zm;
inject=II;
v=VV;
kk=0;
for k=4:4:8*nb
    kk=kk+1;
Ynew2(kk,:)=YY(k,:);
end
end %end OpenDSS powerflow routine


