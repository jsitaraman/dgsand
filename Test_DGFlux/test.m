close all; clc; clear all; 
% % p = 1
% A(1,1) = 0.004883;
% A(1,2) = 0.001790;
% A(1,3) = 0.002441;
% A(2,1) = 0.001790;
% A(2,2) = 0.001628;
% A(2,3) = 0.001790;
% A(3,1) = 0.002441;
% A(3,2) = 0.001790;
% A(3,3) = 0.004883;
% 
% R(1) = 0.000625;
% R(2) = 0.001250;
% R(3) = 0.000625;


% p = 2
nbasis = 6; 
A(1,1) = 0.004883;
A(1,2) = 0.001790;
A(1,3) = 0.002441;
A(1,4) = 0.000846;
A(1,5) = 0.000423;
A(1,6) = 0.001009;
A(2,1) = 0.001790;
A(2,2) = 0.001628;
A(2,3) = 0.001790;
A(2,4) = 0.000521;
A(2,5) = 0.000521;
A(2,6) = 0.000423;
A(3,1) = 0.002441;
A(3,2) = 0.001790;
A(3,3) = 0.004883;
A(3,4) = 0.000423;
A(3,5) = 0.000846;
A(3,6) = 0.001009;
A(4,1) = 0.000846;
A(4,2) = 0.000521;
A(4,3) = 0.000423;
A(4,4) = 0.000228;
A(4,5) = 0.000114;
A(4,6) = 0.000155;
A(5,1) = 0.000423;
A(5,2) = 0.000521;
A(5,3) = 0.000846;
A(5,4) = 0.000114;
A(5,5) = 0.000228;
A(5,6) = 0.000155;
A(6,1) = 0.001009;
A(6,2) = 0.000423;
A(6,3) = 0.001009;
A(6,4) = 0.000155;
A(6,5) = 0.000155;
A(6,6) = 0.000342;

R(1) = 0.000625;
R(2) = .00125;
R(3) = .000625;
R(4) = .000313;
R(5) = .000313;
R(6) = .000104;

% invert
qref = inv(A)*R'; 

% regularization
Areg = [A; eye(nbasis)];
AregA = Areg'*Areg; 
breg = [R'; zeros(nbasis,1)];
Aregb = Areg'*breg;
qreg = inv(AregA)*Aregb;

% truncated SVD: A = U*S*V'
[U,S,V] = svd(A); 
n = 5;
for i = n+1:length(S)
    S(i,i) = 0.0; 
end
Asvd = U*S*V';
qsvd = inv(Asvd)*R';


