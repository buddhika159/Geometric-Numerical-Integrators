function [q, p, NW_itrs, Max_Error, method_time] = IRK4(gradq, gradp, gradqq, gradqp, gradpp, q_int, p_int, dt, itr, tol, Max_NW_Iter)
% INPUTS ------------------------------------------------------------------
% gradq - Funciton handle for the gradiant vactor of the Hamiltonian w.r.t. q
% gradp - Funciton handle for the gradiant vactor of the Hamiltonian w.r.t. p
% gradqq - Funciton handle for the second derivative matrix of the Hamiltonian w.r.t q and q 
% gradqp - Funciton handle for the second derivative matrix of the Hamiltonian w.r.t q and p 
% gradpp - Funciton handle for the second derivative matrix of the Hamiltonian w.r.t p and p 
% q_int - Initial guess for q
% p_int - Initial guess for p
% dt    - step size
% itr   - Number of iterations
% tol   - tolerence for accuracy of the method
% Max_NW_Iter - Maximum number of Newton Iterations using

% OUTPUTS -----------------------------------------------------------------
% q - Approximates for q for the time frame 0:dt:dt*itr
% p - Approximates for p for the time frame 0:dt:dt*itr



%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

dim = length(q_int);


         
HH =@(z) [gradqq(z(0*dim+1:1*dim),z(1*dim+1:2*dim)) gradqp(z(0*dim+1:1*dim),z(1*dim+1:2*dim))
          gradqp(z(0*dim+1:1*dim),z(1*dim+1:2*dim)) gradpp(z(0*dim+1:1*dim),z(1*dim+1:2*dim))]; 


HHH =@(z) [HH(z(0*dim+1:2*dim))	zeros(2*dim)
            zeros(2*dim) HH(z(2*dim+1:4*dim))];


gradH =@(z) [gradq(z(0*dim+1:1*dim),z(1*dim+1:2*dim))
             gradp(z(0*dim+1:1*dim),z(1*dim+1:2*dim))];

gradHH =@(z) [gradH(z(0*dim+1:2*dim));
              gradH(z(2*dim+1:4*dim))];



z = zeros(2*dim,itr+1);  
Errors = zeros(itr,1);
NW_itrs = 0;

z(:,1) = [q_int;p_int];



J = [zeros(dim)   eye(dim)
      -eye(dim) zeros(dim)];        
 
a11 = 1/4;
a12 = 1/4 - sqrt(3)/6;
a21 = 1/4 + sqrt(3)/6;
a22 = 1/4;

b1 = 1/2;
b2 = 1/2;

JJ = [a11*J a12*J
      a21*J a22*J];

%--------------------------------------------------------------------------
% 4th order 2 stage Implicit Runge Kutta
%--------------------------------------------------------------------------  
tic;   
for i=1:itr 
    z_int = z(:,i);   
    z_sol = [z(:,i);z(:,i)];
    
    error = 1;
    NW_i = 0;
    while (error > tol && NW_i < Max_NW_Iter) 
       NW = - (eye(4*dim) - dt*JJ*HHH(z_sol))\(z_sol - [z_int;z_int] - dt*JJ*gradHH(z_sol));
       z_sol = z_sol + NW;

       error = norm(NW);
        NW_i =  NW_i + 1;
    end

    z(:,i+1) = z_int + dt*[b1*J b2*J]*gradHH(z_sol);
    Errors(i) = error;    
    NW_itrs = NW_itrs + NW_i;   

end
method_time = toc;

q = z(0*dim+1:1*dim,:);
p = z(1*dim+1:2*dim,:);                          
            
Max_Error = max(Errors);


disp('4th order IRK4 done')

end


