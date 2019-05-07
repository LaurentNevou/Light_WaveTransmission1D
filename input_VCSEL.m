%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comment: Why the TMM does not match perfectly to the formula?
%
% -> the reflectivity of the DBR is directly linked to the amount of period. This
% amount is different in both cases. Fixe in the formula while convolute with 
% both miror side in the TMM.
% -> The distance between the cavity mode are differents in both cases because in 
% the formula, the cavity length l3 is fixe while in the TMM, the cavity has an effective 
% length that is larger than l3 since the wave penetrate inside the DBR miror
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VCSEL structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nL=3;
nR=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n1=3;
n2=3.6;
lambda0=940e-9;      % Central wavelength

l1=lambda0/(4*abs(n1));   % thickness at lambda/4
l2=lambda0/(4*abs(n2));   % thickness at lambda/4

n3=3+0.0i;
l3=5 * lambda0/(2*abs(n3));

alpha3=4*pi*imag(n3)./lambda;



layer=[

l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2

l3   n3

l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1
l2   n2
l1   n1

];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Formula computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r12=(n1-n2)/(n1+n2);
r21=(n2-n1)/(n1+n2);

t12=2*n1/(n1+n2);
t21=2*n2/(n1+n2);

D1(1,1,:)= exp(+1i*2*pi*n1*l1./lambda);       %% take care on the sign here
D1(2,2,:)= exp(-1i*2*pi*n1*l1./lambda);       %% take care on the sign here

D2(1,1,:)= exp(+1i*2*pi*n2*l2./lambda);       %% take care on the sign here
D2(2,2,:)= exp(-1i*2*pi*n2*l2./lambda);       %% take care on the sign here

P1=(1/t12)*[1 r12 ; r12 1];
P2=(1/t21)*[1 r21 ; r21 1];

for j=1:length(lambda)
  S(:,:,j)=D2(:,:,j)*P2*D1(:,:,j)*P1;
end

Nperiod = floor( length(layer(:,1))/4 );
for j=1:length(lambda)
  SN(:,:,j)=S(:,:,j)^Nperiod;
end

for j=1:length(lambda)
  Rb(j)=(abs(SN(1,2,j)/SN(2,2,j)))^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TT=1-Rb;
theta=pi;
phi=2*pi*n3*l3./lambda;
delta=2*(phi-theta);


% see the book of Vincenzo Savona for the formula

tfp = TT ./ ( 1 -  Rb .* exp(2i*(phi-theta)) ) ;
Tfp = (abs(tfp)).^2 .* exp(-alpha3*l3);

rfp = - ( (sqrt(Rb)-sqrt(Rb).*exp(2i*(phi-theta))) ./ (1-Rb.*exp(2i*(phi-theta))));
Rfp = (abs(rfp)).^2;

Tf=Tfp;
Rf=Rfp;



