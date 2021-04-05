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

% zz(1)=0 => exp(0)=1
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
% the amplitude is defined in intensity domain therefore, I used the sqrt() for the Electric field
% Nevertheless, it doesn t have any influence
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			A0						B0						A1						B1							A2						B2			.......			  	AN						BN						  AN+1					  BN+1				  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      exp(ik.n0.x0)	 	   exp(-ik.n0.x0)	      -exp(ik.n1.x0)		     -exp(-ik.n1.x0)                0						0								0						0							0						0				  %
%ik.n0.exp(ik.n0.x0)	-ik.n0.exp(-ik.n0.x0)	-ik.n1.exp(ik.n1.x0)		ik.n1.exp(-ik.n1.x0)                0						0								0						0							0						0				  %
%																																																																	  %
%			0						0				   exp(ik.n1.x1)		      exp(-ik.n1.x1)	      -exp(ik.n2.x1)	    	 -exp(-ik.n2.x1)					0						0							0						0				  %
%			0						0			 ik.n1.exp(ik.n1.x1)	   -ik.n1.exp(-ik.n1.x1)	-ik.n2.exp(ik.n2.x1)		ik.n2.exp(-ik.n2.x1)					0						0							0						0				  %
%																																																																	  %
%			.						.						.						.							.						.								.						.							.						.				  %
%			.						.						.						.							.						.								.						.							.						.				  %
%			.						.						.						.							.						.								.						.							.						.				  %
%			0						.				   exp(ik.nj.xj)		      exp(-ik.nj.xj)	     -exp(ik.nj+1.xj)           -exp(-ik.nj+1.xj)					.						.							0						0				  %
%			0						.		 	 ik.nj.exp(ik.nj.xj)	   -ik.nj.exp(-ik.nj.xj) -ik.nj+1.exp(ik.nj+1.xj)      ik.nj.exp(-ik.nj+1.xj)					.						.							0						0				  %
%																																																																	  %
%			.						.						.						.							.						.								.						.							.						.				  %
%			.						.						.						.							.						.								.						.							.						.				  %
%			.						.						.						.							.						.								.						.							.						.				  %
%																																																																	  %
%			0						0						.						.					  exp(ik.nN-1.xN-1)	    	 exp(-ik.nN-1.xN-1)			-exp(ik.nN.xN-1)		-exp(-ik.nN.xN-1)					0						0				  %
%			0						0						.						.			  ik.nN-1.exp(ik.nN-1.xN-1) -ik.nN-1.exp(-ik.nN-1.xN-1)	  -ik.nN.exp(ik.nN.xN-1)   ik.nN.exp(-ik.nN.xN-1)					0						0				  %
%																																																																	  %
%			0						0						.						.							0						0						 exp(ik.nN.xN)			 exp(-ik.nN.xN)	 			-exp(ik.nN+1.xN)	   	   -exp(-ik.nN+1.xN)	  %
%			0						0						.						.							0						0				   ik.nN.exp(ik.nN.xN)	  -ik.nN.exp(-ik.nN.xN)		-ik.nN+1.exp(ik.nN+1.xN)	ik.nN+1.exp(-ik.nN+1.xN)	  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
