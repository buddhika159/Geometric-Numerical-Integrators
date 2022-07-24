clear all
close all
clc

%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

omega = 7;
time_range = 100
dt = 0.01
tol = 1e-10
Max_NW_Iter = 100;
dim = 20;

no_sim = 10;


itr = time_range/dt;


% L = randi([-100 100],1,10)'/5;
% 
% for i=1:length(L)
%     if L(i) == 0
%         L(i) = randi([-20 20],1);
%     end
% end
% 
% 
% X = randi([-10 10],1,10)'/2;
% Y = randi([-10 10],1,10)'/2;
% [L X Y]'


% L = [ -0.5000    0.3000    0.6000    0.7000   -0.2000   -0.8000   -0.9000   -0.3000    0.7000   -0.6000]';
% X = [  3.0000  -10.0000    6.0000    9.0000         0    7.0000   -8.0000    5.0000    9.0000    7.0000]';
% Y = [ -5.0000   -6.0000         0   -2.0000         0   10.0000    2.0000    9.0000         0   -1.0000]';
% [L X Y]'

L = [-14.8000  -18.8000   17.6000   -8.0000   -8.2000   -6.8000   -1.4000    6.0000  -11.0000   13.8000]';
X = [  0.5000    3.5000   -1.5000   -0.5000   -4.5000   -3.5000    1.5000   -2.0000    4.0000   -4.0000]';
Y = [  5.0000    0.5000    2.0000    5.0000   -2.0000   -1.0000   -0.5000    3.0000    3.5000   -4.0000]';
[L X Y]'

sL = sqrt(abs(L));


q_int = sL.*X; 
p_int = sign(L).*sL.*Y;


%--------------------------------------------------------------------------
% Hamiltonian and gradiants of the system
%--------------------------------------------------------------------------
syms q [1 dim] real
syms p [1 dim] real

Q = [q1 q2 q3 q4 q5 q6 q7 q8 q9 q10];
P = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10];

H(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) = ...
    (-1/4).*pi.^(-1).*(L(1).*L(2).*log(((q1.*sL(1).^(-1)+(-1).*q2.*sL(2).^(-1)).^2+(p1.*sign(L(1)).*sL(1).^(-1)+(-1).*p2.*sign(L(2)).*sL(2).^(-1)).^2).^(1/2))+L(1).*L(3).*log(((q1.*sL(1).^(-1)+(-1).*q3.*sL(3).^(-1)).^2+(p1.*sign(L(1)).*sL(1).^(-1)+(-1).*p3.*sign(L(3)).*sL(3).^(-1)).^2).^(1/2))+L(2).*L(3).*log(((q2.*sL(2).^(-1)+(-1).*q3.*sL(3).^(-1)).^2+(p2.*sign(L(2)).*sL(2).^(-1)+(-1).*p3.*sign(L(3)).*sL(3).^(-1)).^2).^(1/2))+L(1).*L(4).*log(((q1.*sL(1).^(-1)+(-1).*q4.*sL(4).^(-1)).^2+(p1.*sign(L(1)).*sL(1).^(-1)+(-1).*p4.*sign(L(4)).*sL(4).^(-1)).^2).^(1/2))+L(2).*L(4).*log(((q2.*sL(2).^(-1)+(-1).*q4.*sL(4).^(-1)).^2+(p2.*sign(L(2)).*sL(2).^(-1)+(-1).*p4.*sign(L(4)).*sL(4).^(-1)).^2).^(1/2))+L(3).*L(4).*log(((q3.*sL(3).^(-1)+(-1).*q4.*sL(4).^(-1)).^2+(p3.*sign(L(3)).*sL(3).^(-1)+(-1).*p4.*sign(L(4)).*sL(4).^(-1)).^2).^(1/2))+L(1).*L(5).*log(((q1.*sL(1).^(-1)+(-1).*q5.*sL(5).^(-1)).^2+(p1.*sign(L(1)).*sL(1).^(-1)+(-1).*p5.*sign(L(5)).*sL(5).^(-1)).^2).^(1/2))+L(2).*L(5).*log(((q2.*sL(2).^(-1)+(-1).*q5.*sL(5).^(-1)).^2+(p2.*sign(L(2)).*sL(2).^(-1)+(-1).*p5.*sign(L(5)).*sL(5).^(-1)).^2).^(1/2))+L(3).*L(5).*log(((q3.*sL(3).^(-1)+(-1).*q5.*sL(5).^(-1)).^2+(p3.*sign(L(3)).*sL(3).^(-1)+(-1).*p5.*sign(L(5)).*sL(5).^(-1)).^2).^(1/2))+L(4).*L(5).*log(((q4.*sL(4).^(-1)+(-1).*q5.*sL(5).^(-1)).^2+(p4.*sign(L(4)).*sL(4).^(-1)+(-1).*p5.*sign(L(5)).*sL(5).^(-1)).^2).^(1/2))+L(1).*L(6).*log(((q1.*sL(1).^(-1)+(-1).*q6.*sL(6).^(-1)).^2+(p1.*sign(L(1)).*sL(1).^(-1)+(-1).*p6.*sign(L(6)).*sL(6).^(-1)).^2).^(1/2))+L(2).*L(6).*log(((q2.*sL(2).^(-1)+(-1).*q6.*sL(6).^(-1)).^2+(p2.*sign(L(2)).*sL(2).^(-1)+(-1).*p6.*sign(L(6)).*sL(6).^(-1)).^2).^(1/2))+L(3).*L(6).*log(((q3.*sL(3).^(-1)+(-1).*q6.*sL(6).^(-1)).^2+(p3.*sign(L(3)).*sL(3).^(-1)+(-1).*p6.*sign(L(6)).*sL(6).^(-1)).^2).^(1/2))+L(4).*L(6).*log(((q4.*sL(4).^(-1)+(-1).*q6.*sL(6).^(-1)).^2+(p4.*sign(L(4)).*sL(4).^(-1)+(-1).*p6.*sign(L(6)).*sL(6).^(-1)).^2).^(1/2))+L(5).*L(6).*log(((q5.*sL(5).^(-1)+(-1).*q6.*sL(6).^(-1)).^2+(p5.*sign(L(5)).*sL(5).^(-1)+(-1).*p6.*sign(L(6)).*sL(6).^(-1)).^2).^(1/2))+L(1).*L(7).*log(((q1.*sL(1).^(-1)+(-1).*q7.*sL(7).^(-1)).^2+(p1.*sign(L(1)).*sL(1).^(-1)+(-1).*p7.*sign(L(7)).*sL(7).^(-1)).^2).^(1/2))+L(2).*L(7).*log(((q2.*sL(2).^(-1)+(-1).*q7.*sL(7).^(-1)).^2+(p2.*sign(L(2)).*sL(2).^(-1)+(-1).*p7.*sign(L(7)).*sL(7).^(-1)).^2).^(1/2))+L(3).*L(7).*log(((q3.*sL(3).^(-1)+(-1).*q7.*sL(7).^(-1)).^2+(p3.*sign(L(3)).*sL(3).^(-1)+(-1).*p7.*sign(L(7)).*sL(7).^(-1)).^2).^(1/2))+L(4).*L(7).*log(((q4.*sL(4).^(-1)+(-1).*q7.*sL(7).^(-1)).^2+(p4.*sign(L(4)).*sL(4).^(-1)+(-1).*p7.*sign(L(7)).*sL(7).^(-1)).^2).^(1/2))+L(5).*L(7).*log(((q5.*sL(5).^(-1)+(-1).*q7.*sL(7).^(-1)).^2+(p5.*sign(L(5)).*sL(5).^(-1)+(-1).*p7.*sign(L(7)).*sL(7).^(-1)).^2).^(1/2))+L(6).*L(7).*log(((q6.*sL(6).^(-1)+(-1).*q7.*sL(7).^(-1)).^2+(p6.*sign(L(6)).*sL(6).^(-1)+(-1).*p7.*sign(L(7)).*sL(7).^(-1)).^2).^(1/2))+L(1).*L(8).*log(((q1.*sL(1).^(-1)+(-1).*q8.*sL(8).^(-1)).^2+(p1.*sign(L(1)).*sL(1).^(-1)+(-1).*p8.*sign(L(8)).*sL(8).^(-1)).^2).^(1/2))+L(2).*L(8).*log(((q2.*sL(2).^(-1)+(-1).*q8.*sL(8).^(-1)).^2+(p2.*sign(L(2)).*sL(2).^(-1)+(-1).*p8.*sign(L(8)).*sL(8).^(-1)).^2).^(1/2))+L(3).*L(8).*log(((q3.*sL(3).^(-1)+(-1).*q8.*sL(8).^(-1)).^2+(p3.*sign(L(3)).*sL(3).^(-1)+(-1).*p8.*sign(L(8)).*sL(8).^(-1)).^2).^(1/2))+L(4).*L(8).*log(((q4.*sL(4).^(-1)+(-1).*q8.*sL(8).^(-1)).^2+(p4.*sign(L(4)).*sL(4).^(-1)+(-1).*p8.*sign(L(8)).*sL(8).^(-1)).^2).^(1/2))+L(5).*L(8).*log(((q5.*sL(5).^(-1)+(-1).*q8.*sL(8).^(-1)).^2+(p5.*sign(L(5)).*sL(5).^(-1)+(-1).*p8.*sign(L(8)).*sL(8).^(-1)).^2).^(1/2))+L(6).*L(8).*log(((q6.*sL(6).^(-1)+(-1).*q8.*sL(8).^(-1)).^2+(p6.*sign(L(6)).*sL(6).^(-1)+(-1).*p8.*sign(L(8)).*sL(8).^(-1)).^2).^(1/2))+L(7).*L(8).*log(((q7.*sL(7).^(-1)+(-1).*q8.*sL(8).^(-1)).^2+(p7.*sign(L(7)).*sL(7).^(-1)+(-1).*p8.*sign(L(8)).*sL(8).^(-1)).^2).^(1/2))+L(1).*L(9).*log(((q1.*sL(1).^(-1)+(-1).*q9.*sL(9).^(-1)).^2+(p1.*sign(L(1)).*sL(1).^(-1)+(-1).*p9.*sign(L(9)).*sL(9).^(-1)).^2).^(1/2))+L(2).*L(9).*log(((q2.*sL(2).^(-1)+(-1).*q9.*sL(9).^(-1)).^2+(p2.*sign(L(2)).*sL(2).^(-1)+(-1).*p9.*sign(L(9)).*sL(9).^(-1)).^2).^(1/2))+L(3).*L(9).*log(((q3.*sL(3).^(-1)+(-1).*q9.*sL(9).^(-1)).^2+(p3.*sign(L(3)).*sL(3).^(-1)+(-1).*p9.*sign(L(9)).*sL(9).^(-1)).^2).^(1/2))+L(4).*L(9).*log(((q4.*sL(4).^(-1)+(-1).*q9.*sL(9).^(-1)).^2+(p4.*sign(L(4)).*sL(4).^(-1)+(-1).*p9.*sign(L(9)).*sL(9).^(-1)).^2).^(1/2))+L(5).*L(9).*log(((q5.*sL(5).^(-1)+(-1).*q9.*sL(9).^(-1)).^2+(p5.*sign(L(5)).*sL(5).^(-1)+(-1).*p9.*sign(L(9)).*sL(9).^(-1)).^2).^(1/2))+L(6).*L(9).*log(((q6.*sL(6).^(-1)+(-1).*q9.*sL(9).^(-1)).^2+(p6.*sign(L(6)).*sL(6).^(-1)+(-1).*p9.*sign(L(9)).*sL(9).^(-1)).^2).^(1/2))+L(7).*L(9).*log(((q7.*sL(7).^(-1)+(-1).*q9.*sL(9).^(-1)).^2+(p7.*sign(L(7)).*sL(7).^(-1)+(-1).*p9.*sign(L(9)).*sL(9).^(-1)).^2).^(1/2))+L(8).*L(9).*log(((q8.*sL(8).^(-1)+(-1).*q9.*sL(9).^(-1)).^2+(p8.*sign(L(8)).*sL(8).^(-1)+(-1).*p9.*sign(L(9)).*sL(9).^(-1)).^2).^(1/2))+L(1).*L(10).*log(((q1.*sL(1).^(-1)+(-1).*q10.*sL(10).^(-1)).^2+(p1.*sign(L(1)).*sL(1).^(-1)+(-1).*p10.*sign(L(10)).*sL(10).^(-1)).^2).^(1/2))+L(2).*L(10).*log(((q2.*sL(2).^(-1)+(-1).*q10.*sL(10).^(-1)).^2+(p2.*sign(L(2)).*sL(2).^(-1)+(-1).*p10.*sign(L(10)).*sL(10).^(-1)).^2).^(1/2))+L(3).*L(10).*log(((q3.*sL(3).^(-1)+(-1).*q10.*sL(10).^(-1)).^2+(p3.*sign(L(3)).*sL(3).^(-1)+(-1).*p10.*sign(L(10)).*sL(10).^(-1)).^2).^(1/2))+L(4).*L(10).*log(((q4.*sL(4).^(-1)+(-1).*q10.*sL(10).^(-1)).^2+(p4.*sign(L(4)).*sL(4).^(-1)+(-1).*p10.*sign(L(10)).*sL(10).^(-1)).^2).^(1/2))+L(5).*L(10).*log(((q5.*sL(5).^(-1)+(-1).*q10.*sL(10).^(-1)).^2+(p5.*sign(L(5)).*sL(5).^(-1)+(-1).*p10.*sign(L(10)).*sL(10).^(-1)).^2).^(1/2))+L(6).*L(10).*log(((q6.*sL(6).^(-1)+(-1).*q10.*sL(10).^(-1)).^2+(p6.*sign(L(6)).*sL(6).^(-1)+(-1).*p10.*sign(L(10)).*sL(10).^(-1)).^2).^(1/2))+L(7).*L(10).*log(((q7.*sL(7).^(-1)+(-1).*q10.*sL(10).^(-1)).^2+(p7.*sign(L(7)).*sL(7).^(-1)+(-1).*p10.*sign(L(10)).*sL(10).^(-1)).^2).^(1/2))+L(8).*L(10).*log(((q8.*sL(8).^(-1)+(-1).*q10.*sL(10).^(-1)).^2+(p8.*sign(L(8)).*sL(8).^(-1)+(-1).*p10.*sign(L(10)).*sL(10).^(-1)).^2).^(1/2))+L(9).*L(10).*log(((q9.*sL(9).^(-1)+(-1).*q10.*sL(10).^(-1)).^2+(p9.*sign(L(9)).*sL(9).^(-1)+(-1).*p10.*sign(L(10)).*sL(10).^(-1)).^2).^(1/2)));

gradq  = matlabFunction(gradient(H, Q));
gradp  = matlabFunction(gradient(H, P));
gradqq = matlabFunction(jacobian(gradient(H,Q),Q));
gradqp = matlabFunction(jacobian(gradient(H,Q),P));
gradpp = matlabFunction(jacobian(gradient(H,P),P));

gradq =@(q,p) gradq (q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),q(9),q(10),p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10));
gradp =@(q,p) gradp (q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),q(9),q(10),p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10));
gradqq=@(q,p) gradqq(q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),q(9),q(10),p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10));
gradqp=@(q,p) gradqp(q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),q(9),q(10),p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10));
gradpp=@(q,p) gradpp(q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),q(9),q(10),p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10));


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






