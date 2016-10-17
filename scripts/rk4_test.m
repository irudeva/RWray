

umx=1.29777516e-06;
vmx=1.30e-06;

umy=2.44499843e-05;
vmy=2.20e-08;

qxy=1.30070776e-16;
qxx=8.22972235e-16;

qyy=1.44923526e-15;
qxy=1.30070776e-16;

k=1.62498225546e-07;
l=-1.55155582347e-07;

Ks2=k*k+l*l;

dkdt=-k*umx-l*vmx+(qxy*k-qxx*l)/Ks2;
dldt=-k*umy-l*vmy+(qyy*k-qxy*l)/Ks2;

[t,k1]=ode23(@(t,k1) dkdt,[0 600],k)         
[t,l1]=ode23(@(t,k1) dldt,[0 600],l)

um=14.87;
vm=-1.68;

qx=-1.17e-11;
qy=1.20e-11;

omega = um*k+vm*l+(qx*l-qy*k)/Ks2
          
          