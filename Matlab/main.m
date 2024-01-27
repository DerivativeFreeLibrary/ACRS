clear all;
close all;
clc;

tol = 1.d-4;
maxit = 100000000;
lb=[-10;-10;-10;-10];
ub=[ 10; 10; 10; 10];

% random sequence initialization
rng(137885)

options = struct('WKdim',25*10,'eps',1.0d-0,'maxiter',3000,'omega',100,'verbose',1);

[pout,fout,iter,nf,tcpu,diff]=acrs(@powell,lb,ub,[],options);

