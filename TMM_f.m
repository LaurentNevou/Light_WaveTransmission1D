%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Igor A. Sukhoivanov and Igor V. Guryev
% Photonic Cristals: Physical and Practical Modeling
% Chap3: Fundamentals of Computation of Photonic Crystal Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[A,B,psi]=TMM_f(zz,zv,nt,nL,nR,lambda)

AmplitudeInput=1;
k=2*pi/lambda;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left bondary condition
M=[];

M(1,1:3)=[1 -1 -1];
M(2,1:3) = [-1i*k*nL -1i*k*nt(1)  1i*k*nt(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filling the matrix

for j=1:length(nt)-1

  M(j*2+1,2*j:2*j+3) = [exp(1i*k*nt(j)*zz(j)) exp(-1i*k*nt(j)*zz(j)) -exp(1i*k*nt(j+1)*zz(j)) -exp(-1i*k*nt(j+1)*zz(j))];

  M(j*2+2,2*j:2*j+3) = ...
  [1i*k*nt(j)*exp(1i*k*nt(j)*zz(j)) -1i*k*nt(j)*exp(-1i*k*nt(j)*zz(j)) -1i*k*nt(j+1)*exp(1i*k*nt(j+1)*zz(j)) 1i*k*nt(j+1)*exp(-1i*k*nt(j+1)*zz(j))];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Right bondary condition

M(length(nt)*2+1,2*length(nt):2*length(nt)+2)= [exp(1i*k*nt(end)*zz(j+1)) exp(-1i*k*nt(end)*zz(j+1)) -exp(1i*k*nR*zz(j+1))];

M(length(nt)*2+2,2*length(nt):2*length(nt)+2)= [1i*k*nt(end)*exp(1i*k*nt(end)*zz(j+1)) -1i*k*nt(end)*exp(-1i*k*nt(end)*zz(j+1)) -1i*k*nR*exp(1i*k*nR*zz(j+1))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=zeros(length(M),1);
D(1)=-sqrt(AmplitudeInput);
D(2)=-sqrt(AmplitudeInput)*1i*k*nL;

AB=inv(M)*D;
A=[1 ; AB(2:2:end)];
B=[AB(1:2:end-1) ; 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psi=[];
for j=1:length(nt)
  
  psi= [ psi  A(j+1)*exp(1i*k*nt(j)*zv{j}) + B(j+1)*exp(-1i*k*nt(j)*zv{j}) ];
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%