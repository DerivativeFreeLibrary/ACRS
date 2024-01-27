function [x,fval,iter,nf,tcpu,diff]=acrs(f,L,U,xin,options)
%[x,fval,iter,nf,tcpu]=acrs_lsq(f,L,U,options,model,y,w,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acrs is a global optimization algorithm based on the paper  P. Brachetti,
% G. Di Pillo, M. De Felice Ciccoli, S. Lucidi �A New Version of the 
% Price�s Algorithm for Global Optimization� Journal of Global Optimization (1997) 10: 165-184
%INPUT 
% f = function to be minimized
% L lower bound, it is a vector 1xN (the algorithm needs a lower and upper bound, they cannot be infinity) 
% U upper bound, it is a vector 1xN
% xin can contain a good point that the user already knows. It can be the empty vector if no point is known 
% options is a struct to set:
%   M = dimension of the working set used by ACRS. The default is 25*N
%   eps = tolerance in the stopping criterion. The default is 1.0d-4
%   maxiter = maximum allowed number of iterations. The default is 10000.
%   verbose = it can be 0 or 1, the default is 1. If it is zero the algorithm does not print 
%   anything during the execution. If it is 1, then the algorithm prints at each iteration the iteration number, 
%   the current number of function evaluations, the minimum function value
%   contained in the working set, the maximum function value contained in
%   the working set, the difference between maximum and minimum value in
%   the working set.
% If it is empty, then the default value are used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% x best point found
% fval objective function value in x
% iter number of iterations required to find x
% nf number of function evaluations required to find x
% tcpu cpu time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(L);
tc=cputime;
nf =0;
if ((isempty(L))|(isempty(U))|(nargin<5)),
    disp('error in the input parameters');
    return;
end

if ((isempty(options))),
    
    M=25*N;
    eps=1.0d-4;
    K=10000;
    omega = N;
    verbose =1; 

else
    
    M=options.WKdim;
    eps=options.eps;
    K=options.maxiter;
    omega = options.omega;
    verbose =options.verbose;

end

r=rand(N,M);

nMinLoc=0;
for i=1:N
    for j=1:M
    S(i,j)=[L(i)+(U(i)-L(i))*r(i,j)];
    end
end

if (~isempty(xin)),
    S(:,1)=xin;
end
 for i=1:M
    fprintf('%4d ',i);
    if(verbose >= 2)
        for j = 1:N
            fprintf('%12.3e ',S(j,i:i));
        end
    end
    F(i)=f(S(:,i:i));

    if((F(i) > 9.e+3) || (F(i) < 0.1))
        fprintf('%15.6e\n',F(i));
    else
        fprintf('%15.6f\n',F(i));
    end
    %fprintf('%15.6e\n',F(i));
    nf=nf+1;
 end
[fmax0,Ikmax]=max(F);
[fmin0,Ikmin]=min(F);
fmin=fmin0;
fmax=fmax0;
xkmax=S(:,Ikmax:Ikmax);
xkmin=S(:,Ikmin:Ikmin);
k=0;
if (verbose > 0), 
    disp('iter       nf            fmin             fmax-fmin');
end
while (((fmax-fmin)>=eps)& (k<K))


    k=k+1;
    p=randperm(M);
    
    for i=1:N+1
       	if i<N+1
        	Sk(:,i:i)=[S(:,p(i):p(i))];
            Fv(i)=F(p(i));
            %Fv(i)=F(i);
      	else if ~(Ikmax == p(N+1))
         	Sk(:,N+1:N+1)=[S(:,Ikmax:Ikmax)];
            Fv(N+1)=F(Ikmax);
            else
                Sk(:,i:i)=[S(:,p(i):p(i))];
                Fv(i)=F(p(i));
            %   Fv(i)=F(i);
            end
		end

    end
    fhik=[omega*(((fmax-fmin)^2)/(fmax0-fmin0))];
    for j=1:N+1 
       n(j)=1/(Fv(j)-fmin+fhik);
    end
    for j=1:N+1
 
        w(j)=[n(j)/sum(n)];
         
     end
        
     ck=(w*Sk')'; 
     %[fmaxn1,Ikmaxn1]=max(Fv);
     alfak=1-((fmax-w*Fv')/(fmax-fmin+fhik));
     %ribaltamento va proiettato sui box
     xk=ck -(alfak*(S(:,Ikmax:Ikmax)-ck));
     for i=1:N
         if (xk(i)>U(i)),
             xk(i)=U(i);
         else if (xk(i)<L(i)),
                 xk(i)=L(i);
             end
         end
     end
     alk = rand(1);

     fxk=f(xk);
     nf=nf+1;

     if (fxk<=fmax)
         S(:,Ikmax:Ikmax)=xk;
         F(Ikmax)=fxk   ;
     end
    [fmax,Ikmax]=max(F);
    [fmin,Ikmin]=min(F);
    P(k)=fmin;
%     if (k>1),
%         if ~(P(k) == P(k-1)),
%             fid=fopen('bestpoint.txt','w');
%             
%             if((fmin > 9.e+3) || (fmin < 0.1))
%                 fprintf(fid,'best f: %15.6e\n',fmin);
%             else
%                 fprintf(fid,'best f: %15.6f\n',fmin);
%             end
%             fprintf(fid,'best x: ');
%             for j = 1:N
%                 fprintf(fid,'%15.6e ',xk(j));
%             end
%             fprintf(fid,'\n');
%     
%             %fprintf(fid,'best f %6.2f \n',fmin);
%             %fprintf(fid,'best point %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',xk);
%             fclose(fid);
%         end
%     end
    iter=k;
    if (verbose > 0), 

        fprintf(' %5d  %5d  ',iter,nf);
        if((fmin > 9.e+3) || (fmin < 0.1))
            fprintf('%15.6e  ',fmin);
        else
            fprintf('%15.6f  ',fmin);
        end
        if((fmax > 9.e+3) || (fmax < 0.1))
            fprintf('%15.6e  ',fmax);
        else
            fprintf('%15.6f  ',fmax);
        end
        if((fmax-fmin > 9.e+3) || (fmax-fmin < 0.1))
            fprintf('%15.6e\n',fmax-fmin);
        else
            fprintf('%15.6f\n',fmax-fmin);
        end
        % disp([num2str(iter),'  ', num2str(nf), '  ' , num2str(fmin), '  ' , num2str(fmax), '  ', num2str(fmax-fmin)]);
    
    end
end
if ((fmax-fmin)<eps),
    disp('Stopping criterion satisfied')
elseif (k>=K)
    disp('Maximum number of iterations reached')
end

 x=S(:,Ikmin:Ikmin);
 fval = F(Ikmin);
 numminloc=nMinLoc+M+iter;
 diff=fmax-fmin;
 tcpu = cputime-tc;
