mvc=abs(vc);
mvcm=abs(vcm);
kk=0;
for k=4:4:4*nb
  kk=kk+1;
 xmvc(kk)=mvc(k)*1000; 
  xmvcm(kk)=mvcm(k)*1000; 
end
plot(xmvc)
 hold
 plot(xmvcm)