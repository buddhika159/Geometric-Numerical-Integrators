 clear all
close all
clc

%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

omega = 7;
time_range = 10000
dt = 0.1
tol = 1e-13
Max_NW_Iter = 100;
dim = 10;


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


L = [ -0.5000    0.3000    0.6000    0.7000   -0.2000   -0.8000   -0.9000   -0.3000    0.7000   -0.6000]';
X = [  3.0000  -10.0000    6.0000    9.0000         0    7.0000   -8.0000    5.0000    9.0000    7.0000]';
Y = [ -5.0000   -6.0000         0   -2.0000         0   10.0000    2.0000    9.0000         0   -1.0000]';
[L X Y]'

% L = [-14.8000  -18.8000   17.6000   -8.0000   -8.2000   -6.8000   -1.4000    6.0000  -11.0000   13.8000]';
% X = [  0.5000    3.5000   -1.5000   -0.5000   -4.5000   -3.5000    1.5000   -2.0000    4.0000   -4.0000]';
% Y = [  5.0000    0.5000    2.0000    5.0000   -2.0000   -1.0000   -0.5000    3.0000    3.5000   -4.0000]';
% [L X Y]'

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
% Invariants of the system
%--------------------------------------------------------------------------
H  = matlabFunction(H);
Q=@(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ...
    q1.*L(1).*sL(1).^(-1)+q2.*L(2).*sL(2).^(-1)+q3.*L(3).*sL(3).^(-1)+q4.*L(4).*sL(4).^(-1)+q5.*L(5).*sL(5).^(-1)+q6.*L(6).*sL(6).^(-1)+q7.*L(7).*sL(7).^(-1)+q8.*L(8).*sL(8).^(-1)+q9.*L(9).*sL(9).^(-1)+q10.*L(10).*sL(10).^(-1);
P=@(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ...
    p1.*L(1).*sign(L(1)).*sL(1).^(-1)+p2.*L(2).*sign(L(2)).*sL(2).^(-1)+p3.*L(3).*sign(L(3)).*sL(3).^(-1)+p4.*L(4).*sign(L(4)).*sL(4).^(-1)+p5.*L(5).*sign(L(5)).*sL(5).^(-1)+p6.*L(6).*sign(L(6)).*sL(6).^(-1)+p7.*L(7).*sign(L(7)).*sL(7).^(-1)+p8.*L(8).*sign(L(8)).*sL(8).^(-1)+p9.*L(9).*sign(L(9)).*sL(9).^(-1)+p10.*L(10).*sign(L(10)).*sL(10).^(-1);
I=@(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ...
    L(1).*(p1.^2.*sL(1).^(-2)+q1.^2.*sL(1).^(-2))+L(2).*(p2.^2.*sL(2).^(-2)+q2.^2.*sL(2).^(-2))+L(3).*(p3.^2.*sL(3).^(-2)+q3.^2.*sL(3).^(-2))+L(4).*(p4.^2.*sL(4).^(-2)+q4.^2.*sL(4).^(-2))+L(5).*(p5.^2.*sL(5).^(-2)+q5.^2.*sL(5).^(-2))+L(6).*(p6.^2.*sL(6).^(-2)+q6.^2.*sL(6).^(-2))+L(7).*(p7.^2.*sL(7).^(-2)+q7.^2.*sL(7).^(-2))+L(8).*(p8.^2.*sL(8).^(-2)+q8.^2.*sL(8).^(-2))+L(9).*(p9.^2.*sL(9).^(-2)+q9.^2.*sL(9).^(-2))+L(10).*(p10.^2.*sL(10).^(-2)+q10.^2.*sL(10).^(-2));



%--------------------------------------------------------------------------
% Geometric integrators
%--------------------------------------------------------------------------

Tao_args    = {gradq, gradp, q_int, p_int, dt, itr, omega};
Strang_args = {gradq, gradp, q_int, p_int, dt, itr};
semi_args   = {gradq, gradp, q_int, p_int, dt, itr, tol, Max_NW_Iter};
Gauss_args  = {gradq, gradp, gradqq, gradqp, gradpp, q_int, p_int, dt, itr, tol, Max_NW_Iter};

semi = semiexplicit;

% [ q_1,  p_1,                   Defect_1,  time_1] = Tao_method.base2  (Tao_args{:});
[ q_2,  p_2,                              time_2] = Strang_proj.base2 (Strang_args{:});  
[ q_3,  p_3,  NW_3,  Error_3,  Defect_3,  time_3] = semi.base2        (semi_args{:});  
[ q_4,  p_4,  NW_4,  Error_4,             time_4] = IMid2             (Gauss_args{:});

% [ q_5,  p_5,                   Defect_5,  time_5] = Tao_method.TJ4    (Tao_args{:});
[ q_6,  p_6,                              time_6] = Strang_proj.TJ4   (Strang_args{:}); 
[ q_7,  p_7,  NW_7,  Error_7,  Defect_7,  time_7] = semi.TJ4          (semi_args{:});   
[ q_8,  p_8,  NW_8,  Error_8,  Defect_8,  time_8] = semi.Su4          (semi_args{:});   
[ q_9,  p_9,  NW_9,  Error_9,             time_9] = IRK4              (Gauss_args{:});  
% [ q_9,  p_9,                              time_9] = RK4               (Strang_args{:}); 

% [q_10, p_10,                  Defect_10, time_10] = Tao_method.TJ6    (Tao_args{:});
[q_11, p_11,                             time_11] = Strang_proj.TJ6   (Strang_args{:});  
[q_12, p_12, NW_12, Error_12, Defect_12, time_12] = semi.TJ6          (semi_args{:});   
[q_13, p_13, NW_13, Error_13, Defect_13, time_13] = semi.Su6          (semi_args{:});  
[q_14, p_14, NW_14, Error_14, Defect_14, time_14] = semi.Yo6          (semi_args{:});  


% 
% method = [{'Tao 2'};{'Strang_proj 2'};{'semiexplicit 2'};{'Midpoint'};
%           {'Tao 4'};{'Strang_proj 4'};{'semiexplicit 4'};{'semiexplicit-S 4'};{'IRK4'};
%           {'Tao 6'};{'Strang_proj 6'};{'semiexplicit 6'};{'semiexplicit-S 6'};{'semiexplicit-Y 6'}];
% 
% NW_itrs     = [     nan;     nan;    NW_3;    NW_4;     nan;     nan;    NW_7;    NW_8;    NW_9;      nan;     nan;    NW_12;    NW_13;    NW_14]/itr;
% Max_Errors  = [     nan;     nan; Error_3; Error_4;     nan;     nan; Error_7; Error_8; Error_9;      nan;     nan; Error_12; Error_13; Error_14];
% Max_Defects = [Defect_1;     nan;Defect_3;     nan;Defect_5;     nan;Defect_7;Defect_8;     nan;Defect_10;     nan;Defect_12;Defect_13;Defect_14];
% time_cal    = [  time_1;  time_2;  time_3;  time_4;  time_5;  time_6;  time_7;  time_8;  time_9;  time_10; time_11;  time_12;  time_13;  time_14];
% 
% 
% table(method, time_cal)
% 
% 


%--------------------------------------------------------------------------
% Computation of errors of the system invariants
%--------------------------------------------------------------------------

z_int = num2cell([q_int;p_int]);

H_int = H(z_int{:});
Q_int = Q(z_int{:});
P_int = P(z_int{:});
I_int = I(z_int{:});


% z_1  = [arrayfun(@(x) {q_1(x,:)},1:dim), arrayfun(@(x) {p_1(x,:)},1:dim)];
z_2  = [arrayfun(@(x) {q_2(x,:)},1:dim), arrayfun(@(x) {p_2(x,:)},1:dim)];
z_3  = [arrayfun(@(x) {q_3(x,:)},1:dim), arrayfun(@(x) {p_3(x,:)},1:dim)];
z_4  = [arrayfun(@(x) {q_4(x,:)},1:dim), arrayfun(@(x) {p_4(x,:)},1:dim)];
% z_5  = [arrayfun(@(x) {q_5(x,:)},1:dim), arrayfun(@(x) {p_5(x,:)},1:dim)];
z_6  = [arrayfun(@(x) {q_6(x,:)},1:dim), arrayfun(@(x) {p_6(x,:)},1:dim)];
z_7  = [arrayfun(@(x) {q_7(x,:)},1:dim), arrayfun(@(x) {p_7(x,:)},1:dim)];
z_8  = [arrayfun(@(x) {q_8(x,:)},1:dim), arrayfun(@(x) {p_8(x,:)},1:dim)];
z_9  = [arrayfun(@(x) {q_9(x,:)},1:dim), arrayfun(@(x) {p_9(x,:)},1:dim)];
% z_10 = [arrayfun(@(x) {q_10(x,:)},1:dim), arrayfun(@(x) {p_10(x,:)},1:dim)];
z_11 = [arrayfun(@(x) {q_11(x,:)},1:dim), arrayfun(@(x) {p_11(x,:)},1:dim)];
z_12 = [arrayfun(@(x) {q_12(x,:)},1:dim), arrayfun(@(x) {p_12(x,:)},1:dim)];
z_13 = [arrayfun(@(x) {q_13(x,:)},1:dim), arrayfun(@(x) {p_13(x,:)},1:dim)];
z_14 = [arrayfun(@(x) {q_14(x,:)},1:dim), arrayfun(@(x) {p_14(x,:)},1:dim)];

% H_Tao2_error         = (H(z_1{:}) - H_int)/H_int;
H_Strang_proj2_error = (H(z_2{:}) - H_int)/H_int;
H_semi2_error        = (H(z_3{:}) - H_int)/H_int;
H_mid_error          = (H(z_4{:}) - H_int)/H_int;
% H_Tao4_error         = (H(z_5{:}) - H_int)/H_int;
H_Strang_proj4_error = (H(z_6{:}) - H_int)/H_int;
H_Triple4_error      = (H(z_7{:}) - H_int)/H_int;
H_Suzuki4_error      = (H(z_8{:}) - H_int)/H_int;
H_IRK4_error         = (H(z_9{:}) - H_int)/H_int;
% H_Tao6_error         = (H(z_10{:}) - H_int)/H_int;
H_Strang_proj6_error = (H(z_11{:}) - H_int)/H_int;
H_Triple6_error      = (H(z_12{:}) - H_int)/H_int;
H_Suzuki6_error      = (H(z_13{:}) - H_int)/H_int;
H_Yoshida6_error     = (H(z_14{:}) - H_int)/H_int;


% Q_Tao2_error         = (Q(z_1{:}) - Q_int)/Q_int;
Q_Strang_proj2_error = (Q(z_2{:}) - Q_int)/Q_int;
Q_semi2_error        = (Q(z_3{:}) - Q_int)/Q_int;
Q_mid_error          = (Q(z_4{:}) - Q_int)/Q_int;
% Q_Tao4_error         = (Q(z_5{:}) - Q_int)/Q_int;
Q_Strang_proj4_error = (Q(z_6{:}) - Q_int)/Q_int;
Q_Triple4_error      = (Q(z_7{:}) - Q_int)/Q_int;
Q_Suzuki4_error      = (Q(z_8{:}) - Q_int)/Q_int;
Q_IRK4_error         = (Q(z_9{:}) - Q_int)/Q_int;
% Q_Tao6_error         = (Q(z_10{:}) - Q_int)/Q_int;
Q_Strang_proj6_error = (Q(z_11{:}) - Q_int)/Q_int;
Q_Triple6_error      = (Q(z_12{:}) - Q_int)/Q_int;
Q_Suzuki6_error      = (Q(z_13{:}) - Q_int)/Q_int;
Q_Yoshida6_error     = (Q(z_14{:}) - Q_int)/Q_int;


% P_Tao2_error         = (P(z_1{:}) - P_int)/P_int;
P_Strang_proj2_error = (P(z_2{:}) - P_int)/P_int;
P_semi2_error        = (P(z_3{:}) - P_int)/P_int;
P_mid_error          = (P(z_4{:}) - P_int)/P_int;
% P_Tao4_error         = (P(z_5{:}) - P_int)/P_int;
P_Strang_proj4_error = (P(z_6{:}) - P_int)/P_int;
P_Triple4_error      = (P(z_7{:}) - P_int)/P_int;
P_Suzuki4_error      = (P(z_8{:}) - P_int)/P_int;
P_IRK4_error         = (P(z_9{:}) - P_int)/P_int;
% P_Tao6_error         = (P(z_10{:}) - P_int)/P_int;
P_Strang_proj6_error = (P(z_11{:}) - P_int)/P_int;
P_Triple6_error      = (P(z_12{:}) - P_int)/P_int;
P_Suzuki6_error      = (P(z_13{:}) - P_int)/P_int;
P_Yoshida6_error     = (P(z_14{:}) - P_int)/P_int;


% I_Tao2_error         = (I(z_1{:}) - I_int)/I_int;
I_Strang_proj2_error = (I(z_2{:}) - I_int)/I_int;
I_semi2_error        = (I(z_3{:}) - I_int)/I_int;
I_mid_error          = (I(z_4{:}) - I_int)/I_int;
% I_Tao4_error         = (I(z_5{:}) - I_int)/I_int;
I_Strang_proj4_error = (I(z_6{:}) - I_int)/I_int;
I_Triple4_error      = (I(z_7{:}) - I_int)/I_int;
I_Suzuki4_error      = (I(z_8{:}) - I_int)/I_int;
I_IRK4_error         = (I(z_9{:}) - I_int)/I_int;
% I_Tao6_error         = (I(z_10{:}) - I_int)/I_int;
I_Strang_proj6_error = (I(z_11{:}) - I_int)/I_int;
I_Triple6_error      = (I(z_12{:}) - I_int)/I_int;
I_Suzuki6_error      = (I(z_13{:}) - I_int)/I_int;
I_Yoshida6_error     = (I(z_14{:}) - I_int)/I_int;




%--------------------------------------------------------------------------
% Figures
%--------------------------------------------------------------------------
time = 0:dt:dt*itr;



figure

set(0,'DefaultAxesFontSize', 12)

subplot(3,4,1)
plot(time,H_Strang_proj2_error); hold on
plot(time,H_semi2_error)
xlabel('time','interpreter','latex')
legend('Strang-proj 2','semiexplicit 2')

title('$\frac{H-H_{0}}{H_{0}}$','interpreter','latex','FontSize', 16)

subplot(3,4,2)
plot(time,Q_Strang_proj2_error); hold on
plot(time,Q_semi2_error)
xlabel('time','interpreter','latex')


title('$\frac{Q-Q_{0}}{Q_{0}}$','interpreter','latex','FontSize', 16)

subplot(3,4,3)
plot(time,P_Strang_proj2_error); hold on
plot(time,P_semi2_error)
xlabel('time','interpreter','latex')

 
title('$\frac{P-P_{0}}{P_{0}}$','interpreter','latex','FontSize', 16)

subplot(3,4,4)
plot(time,I_Strang_proj2_error); hold on
plot(time,I_semi2_error)
xlabel('time','Interpreter','latex')


title('$\frac{I-I_{0}}{I_{0}}$','interpreter','latex','FontSize', 16)


subplot(3,4,5)
plot(time,H_Strang_proj4_error); hold on
plot(time,H_Triple4_error) 
xlabel('time','Interpreter','latex')
legend('Strang-proj 4','semiexplicit 4');


subplot(3,4,6)
plot(time,Q_Strang_proj4_error); hold on
plot(time,Q_Triple4_error) 
xlabel('time','Interpreter','latex')


subplot(3,4,7)
plot(time,P_Strang_proj4_error); hold on
plot(time,P_Triple4_error) 
xlabel('time','Interpreter','latex')


subplot(3,4,8)
plot(time,I_Strang_proj4_error); hold on
plot(time,I_Triple4_error) 
xlabel('time','Interpreter','latex')


subplot(3,4,9)
plot(time,H_Strang_proj6_error); hold on
plot(time,H_Triple6_error)
xlabel('time','Interpreter','latex')
legend('Strang-proj 6','semiexplicit 6')


subplot(3,4,10)
plot(time,Q_Strang_proj6_error); hold on
plot(time,Q_Triple6_error)
xlabel('time','Interpreter','latex')


subplot(3,4,11)
plot(time,P_Strang_proj6_error); hold on
plot(time,P_Triple6_error)
xlabel('time','Interpreter','latex')


subplot(3,4,12)
plot(time,I_Strang_proj6_error); hold on
plot(time,I_Triple6_error)
xlabel('time','Interpreter','latex')


set(findobj(gcf,'type','axes'), ...
    'FontName'    , 'latex'   , ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'     , ...
    'YMinorGrid'  , 'off'     , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         )

set(findall(gcf,'type','legend'), ...
    'FontSize'      , 11        , ...
    'Interpreter'   , 'latex'   , ...
    'ItemTokenSize' , [10,5]    , ...
    'Location'      , 'best'    )


set(gcf,'Position',[150 150 1150 550])






figure

set(0,'DefaultAxesFontSize', 12)

subplot(3,4,1)
plot(time,H_semi2_error); hold on
plot(time,H_mid_error)
xlabel('time','interpreter','latex')
legend('semiexplicit 2','midpoint')

title('$\frac{H-H_{0}}{H_{0}}$','interpreter','latex','FontSize', 16)

subplot(3,4,2)
plot(time,Q_semi2_error); hold on
plot(time,Q_mid_error)
xlabel('time','interpreter','latex')


title('$\frac{Q-Q_{0}}{Q_{0}}$','interpreter','latex','FontSize', 16)

subplot(3,4,3)
plot(time,P_semi2_error); hold on
plot(time,P_mid_error)
xlabel('time','interpreter','latex')

 
title('$\frac{P-P_{0}}{P_{0}}$','interpreter','latex','FontSize', 16)

subplot(3,4,4)
plot(time,I_semi2_error); hold on
plot(time,I_mid_error)
xlabel('time','Interpreter','latex')


title('$\frac{I-I_{0}}{I_{0}}$','interpreter','latex','FontSize', 16)


subplot(3,4,5)
plot(time,H_Triple4_error) ; hold on
plot(time,H_Suzuki4_error)
plot(time,H_IRK4_error)
xlabel('time','Interpreter','latex')
legend('semiexplicit 4','semiexplicit-S 4','IRK4')


subplot(3,4,6)
plot(time,Q_Triple4_error) ; hold on
plot(time,Q_Suzuki4_error)
plot(time,Q_IRK4_error)
xlabel('time','Interpreter','latex')


subplot(3,4,7)
plot(time,P_Triple4_error) ; hold on
plot(time,P_Suzuki4_error)
plot(time,P_IRK4_error)
xlabel('time','Interpreter','latex')


subplot(3,4,8)
plot(time,I_Triple4_error); hold on 
plot(time,I_Suzuki4_error)
plot(time,I_IRK4_error)
xlabel('time','Interpreter','latex')


subplot(3,4,9)
plot(time,H_Triple6_error); hold on
plot(time,H_Suzuki6_error)
plot(time,H_Yoshida6_error)
xlabel('time','Interpreter','latex')
legend('semiexplicit 6','semiexplicit-S 6','semiexplicit-Y 6')


subplot(3,4,10)
plot(time,Q_Triple6_error); hold on
plot(time,Q_Suzuki6_error)
plot(time,Q_Yoshida6_error)
xlabel('time','Interpreter','latex')


subplot(3,4,11)
plot(time,P_Triple6_error); hold on
plot(time,P_Suzuki6_error)
plot(time,P_Yoshida6_error)
xlabel('time','Interpreter','latex')


subplot(3,4,12)
plot(time,I_Triple6_error); hold on
plot(time,I_Suzuki6_error)
plot(time,I_Yoshida6_error)
xlabel('time','Interpreter','latex')


set(findobj(gcf,'type','axes'), ...
    'FontName'    , 'latex'   , ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'     , ...
    'YMinorGrid'  , 'off'     , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         )

set(findall(gcf,'type','legend'), ...
    'FontSize'      , 11        , ...
    'Interpreter'   , 'latex'   , ...
    'ItemTokenSize' , [10,5]    , ...
    'Location'      , 'best'    )


set(gcf,'Position',[150 200 1150 550])





% figure
% 
% set(0,'DefaultAxesFontSize', 12)
% 
% subplot(3,4,1)
% plot(time,H_Tao2_error); hold on
% plot(time,H_Strang_proj2_error)
% xlabel('time','Interpreter','latex')
% legend('Tao 2','Strang-proj 2')
% 
% title('$\frac{H-H_{0}}{H_{0}}$','interpreter','latex','FontSize', 16)
% 
% 
% subplot(3,4,2)
% plot(time,Q_Tao2_error); hold on
% plot(time,Q_Strang_proj2_error)
% xlabel('time','Interpreter','latex')
% 
% title('$\frac{Q-Q_{0}}{Q_{0}}$','interpreter','latex','FontSize', 16)
% 
% 
% subplot(3,4,3)
% plot(time,P_Tao2_error); hold on
% plot(time,P_Strang_proj2_error)
% xlabel('time','Interpreter','latex')
% 
% title('$\frac{P-P_{0}}{P_{0}}$','interpreter','latex','FontSize', 16)
% 
% 
% subplot(3,4,4)
% plot(time,I_Tao2_error); hold on
% plot(time,I_Strang_proj2_error)
% xlabel('time','Interpreter','latex')
% 
% title('$\frac{I-I_{0}}{I_{0}}$','interpreter','latex','FontSize', 16)
% 
% 
% subplot(3,4,5)
% plot(time,H_Tao4_error); hold on
% plot(time,H_Strang_proj4_error)
% xlabel('time','Interpreter','latex')
% legend('Tao 4','Strang-proj 4')
% 
% 
% subplot(3,4,6)
% plot(time,Q_Tao4_error); hold on
% plot(time,Q_Strang_proj4_error)
% xlabel('time','Interpreter','latex')
% 
% 
% subplot(3,4,7)
% plot(time,P_Tao4_error); hold on
% plot(time,P_Strang_proj4_error)
% xlabel('time','Interpreter','latex')
% 
% 
% subplot(3,4,8)
% plot(time,I_Tao4_error); hold on
% plot(time,I_Strang_proj4_error)
% xlabel('time','Interpreter','latex')
% 
% 
% subplot(3,4,9)
% plot(time,H_Tao6_error); hold on
% plot(time,H_Strang_proj6_error)
% xlabel('time','Interpreter','latex')
% legend('Tao 6','Strang-proj 6')
% 
% 
% subplot(3,4,10)
% plot(time,Q_Tao6_error); hold on
% plot(time,Q_Strang_proj6_error)
% xlabel('time','Interpreter','latex')
% 
% 
% subplot(3,4,11)
% plot(time,P_Tao6_error); hold on
% plot(time,P_Strang_proj6_error)
% xlabel('time','Interpreter','latex')
% 
% 
% subplot(3,4,12)
% plot(time,I_Tao6_error); hold on
% plot(time,I_Strang_proj6_error)
% xlabel('time','Interpreter','latex')
% 
% 
% set(findobj(gcf,'type','axes'), ...
%     'FontName'    , 'latex'   , ...
%     'Box'         , 'off'     , ...
%     'TickDir'     , 'out'     , ...
%     'TickLength'  , [.02 .02] , ...
%     'XMinorTick'  , 'on'      , ...
%     'YMinorTick'  , 'on'      , ...
%     'YGrid'       , 'off'     , ...
%     'YMinorGrid'  , 'off'     , ...
%     'XColor'      , [.3 .3 .3], ...
%     'YColor'      , [.3 .3 .3], ...
%     'LineWidth'   , 1         )
% 
% set(findall(gcf,'type','legend'), ...
%     'FontSize'      , 11        , ...
%     'Interpreter'   , 'latex'   , ...
%     'ItemTokenSize' , [10,5]    , ...
%     'Location'      , 'best'    )
% 
% 
% set(gcf,'Position',[250 250 1150 550])

% 
% filename = sprintf('10_vortices_%g.mat', abs(log10(tol)));
% path = [pwd, '\', filename];
% save(path,'omega','tol','time_range','dt','table_cols','time','H_invs','Q_invs','P_invs','I_invs')
% 
% 
% 



