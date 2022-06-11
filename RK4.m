function [q, p, method_time] = RK4(gradq, gradp, q_int, p_int, dt, itr)
% INPUTS ------------------------------------------------------------------
% gradq - Funciton handle for the gradiant vactor of the Hamiltonian w.r.t. q
% gradp - Funciton handle for the gradiant vactor of the Hamiltonian w.r.t. p

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


gradH =@(z) [gradq(z(0*dim+1:1*dim),z(1*dim+1:2*dim))
             gradp(z(0*dim+1:1*dim),z(1*dim+1:2*dim))];

J = [zeros(dim)   eye(dim)
      -eye(dim) zeros(dim)];                                                                         
                                                                       

F = @(z) J * gradH(z);  % Hamiltonian Vector Field

q = zeros(dim,itr+1);  
p = zeros(dim,itr+1);  

q(:,1) = q_int;
p(:,1) = p_int;

%--------------------------------------------------------------------------
% 4th order Explicit Runge Kutta
%--------------------------------------------------------------------------  
tic;   
for i=1:itr 
    Z1 = [q(:,i); p(:,i)];
    
    Z2 = Z1 + dt/2 * F(Z1);
    Z3 = Z1 + dt/2 * F(Z2);
    Z4 = Z1 + dt * F(Z3);
    Z  = Z1 + dt/6 * (F(Z1) + 2*F(Z2) + 2*F(Z3) + F(Z4));  % main equation
    
    q(:,i+1) = Z(0*dim+1:1*dim);
    p(:,i+1) = Z(1*dim+1:2*dim);   

end
method_time = toc;



disp('4th order RK4 done')

end


