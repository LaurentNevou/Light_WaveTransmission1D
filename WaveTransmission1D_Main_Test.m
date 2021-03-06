%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% last update 7May2019, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code "_Test" is comparing the full algorithm based on the tranfers matrix
% with a simplified version that I call "formula"

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lambda=(600:2:1000)*1e-9;
%lambda=(850:0.5:1050)*1e-9;
lambda=(600:2:1000)*1e-9;

dz=5e-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Choose your structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input_ARcoating
input_BraggMirror
%input_FabryPerotCavity
%input_VCSEL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotTransmission=0;
% OR
plotReflexion=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z and the optical index n

t  = layer(:,1);
nt = layer(:,2);

for j=1:length(t)
  
  if j==1
    zz(1) = t(1);
    zv{1} = 0:dz:t(1); 
    z     = zv{1};
    n     =(zv{j}*0+1) * nt(j);
  else
    zz(j) = zz(end)+t(j);
    zv{j} = (zz(end-1)+dz):dz:zz(end);
    z     = [ z  zv{j} ];
    n     = [ n   (zv{j}*0+1) * nt(j)  ];
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l=1:length(lambda)

  [AA,BB,psi] = TMM_f(zz,zv,nt,nL,nR,lambda(l));
  
  A(:,l)=AA;
  B(:,l)=BB;
  PSI(:,l)=psi.';
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%X0fig=-1800; Y0fig=100;
X0fig=100; Y0fig=100;
Wfig=1500;Hfig=1000;

figure('Name','Results','position',[X0fig Y0fig Wfig Hfig])

FS=15;
LW=2;
idx=find(abs(lambda-lambda0)==min(abs(lambda-lambda0)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1,'fontsize',FS)
hold on;grid on;

plot(z*1e6,n,'b','linewidth',LW)

xlim([0 z(end)]*1e6)
ylim([0 4])
xlabel('z (um)')
ylabel('optical index')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3,'fontsize',FS)
hold on;grid on;

plot(z*1e6,real(PSI(:,idx)),'b.-','linewidth',LW)
plot(z*1e6,imag(PSI(:,idx)),'g.-','linewidth',LW)

plot(z*1e6,(abs(PSI(:,idx))).^2,'r.-','linewidth',LW)

xlim([0 z(end)]*1e6)
title(strcat('@lambda=',num2str(lambda(idx)*1e9),'nm'))
xlabel('z (um)')
ylabel('Electrical field (a.u.)')
legend('real(E)','imag(E)','|E|^2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(lambda)>1

subplot(1,2,2,'fontsize',FS)
hold on;grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = abs(B(1,:)).^2;
T = (nR/nL) * abs(A(end,:)).^2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotTransmission==1 && plotReflexion==0
  
  plot(lambda*1e9,T,'g-','linewidth',LW)
  plot(lambda*1e9,Tf,'m--','linewidth',LW)
  legend('Transmission: AN+1','Transmission Formula')
  ylabel('Transmission')
  
elseif plotTransmission==0 && plotReflexion==1
  
  plot(lambda*1e9,R,'g-','linewidth',LW)
  plot(lambda*1e9,Rf,'m--','linewidth',LW)
  legend('Reflexion: B0','Reflexion Formula')
  ylabel('Reflexion')
  
elseif (plotTransmission==1 && plotReflexion==1) || (plotTransmission==0 && plotReflexion==0)

  display('ERROR: choose between Transmission OR Reflexion plot')  
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlabel('lambda (nm)')
xlim([lambda(1) lambda(end)]*1e9)
ylim([0 1.15])
title('Spectrum')

end
