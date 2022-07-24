clear all
close all
clc

%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

omega = 100;
time_range = 1000
dt = 0.001
tol = 1e-10
Max_NW_Iter = 10;
dim = 5;

no_sim = 20;

itr = time_range/dt;


q_int = [3;0.01;0.01;0.01;0.01]; 
p_int = [1;   0;   0;   0;   0]; 


%--------------------------------------------------------------------------
% Hamiltonian and gradiants of the system
%--------------------------------------------------------------------------
syms q [1 dim] real
syms p [1 dim] real

Q = [q1 q2 q3 q4 q5];
P = [p1 p2 p3 p4 p5];

H(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5) = (-1).*p1.^2.*p2.^2+(-1).*p2.^2.*p3.^2+(-1).*p3.^2.*p4.^2+(-1).*p4.^2.*p5.^2+p2.^2.*q1.^2+(1/4).*(p1.^2+q1.^2).^2+(-4).*p1.*p2.*q1.*q2+p1.^2.*q2.^2+p3.^2.*q2.^2+(-1).*q1.^2.*q2.^2+(1/4).*(p2.^2+q2.^2).^2+(-4).*p2.*p3.*q2.*q3+p2.^2.*q3.^2+p4.^2.*q3.^2+(-1).*q2.^2.*q3.^2+(1/4).*(p3.^2+q3.^2).^2+(-4).*p3.*p4.*q3.*q4+p3.^2.*q4.^2+p5.^2.*q4.^2+(-1).*q3.^2.*q4.^2+(1/4).*(p4.^2+q4.^2).^2+(-4).*p4.*p5.*q4.*q5+p4.^2.*q5.^2+(-1).*q4.^2.*q5.^2+(1/4).*(p5.^2+q5.^2).^2;

gradq  = matlabFunction(gradient(H, Q));
gradp  = matlabFunction(gradient(H, P));
gradqq = matlabFunction(jacobian(gradient(H,Q),Q));
gradqp = matlabFunction(jacobian(gradient(H,Q),P));
gradpp = matlabFunction(jacobian(gradient(H,P),P));

gradq  = @(q,p) gradq (q(1),q(2),q(3),q(4),q(5),p(1),p(2),p(3),p(4),p(5));
gradp  = @(q,p) gradp (q(1),q(2),q(3),q(4),q(5),p(1),p(2),p(3),p(4),p(5));
gradqq = @(q,p) gradqq(q(1),q(2),q(3),q(4),q(5),p(1),p(2),p(3),p(4),p(5));
gradqp = @(q,p) gradqp(q(1),q(2),q(3),q(4),q(5),p(1),p(2),p(3),p(4),p(5));
gradpp = @(q,p) gradpp(q(1),q(2),q(3),q(4),q(5),p(1),p(2),p(3),p(4),p(5));


%--------------------------------------------------------------------------
% Geometric integrators
%--------------------------------------------------------------------------
method = [{'Tao 2'};{'Tao 4'};{'Tao 6'};{'Strang-proj 2'};{'Strang-proj 4'};{'Strang-proj 6'}];

Tao_args    = {gradq, gradp, q_int, p_int, dt, itr, omega};
Strang_args = {gradq, gradp, q_int, p_int, dt, itr};

time_cals   = zeros(6,no_sim);

warning('off')

for i = 1:no_sim
    [~, ~,                      Defect_1,   time_1] = Tao_method.base2 (Tao_args{:});   
    [~, ~,                      Defect_2,   time_2] = Tao_method.TJ4   (Tao_args{:});   
    [~, ~,                      Defect_3,   time_3] = Tao_method.TJ6   (Tao_args{:});

    [~, ~,                                  time_4] = Strang_proj.base2(Strang_args{:});   
    [~, ~,                                  time_5] = Strang_proj.TJ4  (Strang_args{:});   
    [~, ~,                                  time_6] = Strang_proj.TJ6  (Strang_args{:});

    time_cals(:,i)           = [  time_1;  time_2;  time_3;  time_4;  time_5;  time_6];
    
    table(method, time_cals(:,i))

end

time_cal   = mean(time_cals,2);

tol

table_1 = table(method, time_cal)


% table(method, time_cal, NW_itr, time_cal_Broyden, NW_itr_Broyden)





