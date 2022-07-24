clear all
close all
clc

%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

omega = 20;
time_range = 1000;
dt = 0.1;
tol = 1e-13;
Max_NW_Iter = 100;
dim = 1;


itr = time_range/dt;


q_int = -3;
p_int = 0;


%--------------------------------------------------------------------------
% Hamiltonian and gradiants of the system
%--------------------------------------------------------------------------
syms q [1 dim] real
syms p [1 dim] real

H =@(q,p) (q.^2+1).*(p.^2+1)/2;

gradq = @(q,p) q.*(p.^2+1);
gradp = @(q,p) p.*(q.^2+1);

gradqq = @(q,p) (p.^2+1);
gradpp = @(q,p) (q.^2+1);
gradqp = @(q,p) 2*q.*p;       


%--------------------------------------------------------------------------
% Geometric integrators
%--------------------------------------------------------------------------


Tao_args    = {gradq, gradp, q_int, p_int, dt, itr, omega};
Strang_args = {gradq, gradp, q_int, p_int, dt, itr};
semi_args   = {gradq, gradp, q_int, p_int, dt, itr, tol, Max_NW_Iter};
Gauss_args  = {gradq, gradp, gradqq, gradqp, gradpp, q_int, p_int, dt, itr, tol, Max_NW_Iter};

semi = semiexplicit;

[ q_1,  p_1,                   Defect_1,  time_1] = Tao_method.base2  (Tao_args{:});
[ q_2,  p_2,                              time_2] = Strang_proj.base2 (Strang_args{:});  
[ q_3,  p_3,  NW_3,  Error_3,  Defect_3,  time_3] = semi.base2        (semi_args{:});  
[ q_4,  p_4,  NW_4,  Error_4,             time_4] = IMid2             (Gauss_args{:});

[ q_5,  p_5,                   Defect_5,  time_5] = Tao_method.TJ4    (Tao_args{:});
[ q_6,  p_6,                              time_6] = Strang_proj.TJ4   (Strang_args{:}); 
[ q_7,  p_7,  NW_7,  Error_7,  Defect_7,  time_7] = semi.TJ4          (semi_args{:});   
[ q_8,  p_8,  NW_8,  Error_8,  Defect_8,  time_8] = semi.Su4          (semi_args{:});   
[ q_9,  p_9,  NW_9,  Error_9,             time_9] = IRK4              (Gauss_args{:});  
% [ q_9,  p_9,                              time_9] = RK4               (Strang_args{:}); 

[q_10, p_10,                  Defect_10, time_10] = Tao_method.TJ6    (Tao_args{:});
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



z_1  = [arrayfun(@(x) {q_1(x,:)},1:dim), arrayfun(@(x) {p_1(x,:)},1:dim)];
z_2  = [arrayfun(@(x) {q_2(x,:)},1:dim), arrayfun(@(x) {p_2(x,:)},1:dim)];
z_3  = [arrayfun(@(x) {q_3(x,:)},1:dim), arrayfun(@(x) {p_3(x,:)},1:dim)];
z_4  = [arrayfun(@(x) {q_4(x,:)},1:dim), arrayfun(@(x) {p_4(x,:)},1:dim)];
z_5  = [arrayfun(@(x) {q_5(x,:)},1:dim), arrayfun(@(x) {p_5(x,:)},1:dim)];
z_6  = [arrayfun(@(x) {q_6(x,:)},1:dim), arrayfun(@(x) {p_6(x,:)},1:dim)];
z_7  = [arrayfun(@(x) {q_7(x,:)},1:dim), arrayfun(@(x) {p_7(x,:)},1:dim)];
z_8  = [arrayfun(@(x) {q_8(x,:)},1:dim), arrayfun(@(x) {p_8(x,:)},1:dim)];
z_9  = [arrayfun(@(x) {q_9(x,:)},1:dim), arrayfun(@(x) {p_9(x,:)},1:dim)];
z_10 = [arrayfun(@(x) {q_10(x,:)},1:dim), arrayfun(@(x) {p_10(x,:)},1:dim)];
z_11 = [arrayfun(@(x) {q_11(x,:)},1:dim), arrayfun(@(x) {p_11(x,:)},1:dim)];
z_12 = [arrayfun(@(x) {q_12(x,:)},1:dim), arrayfun(@(x) {p_12(x,:)},1:dim)];
z_13 = [arrayfun(@(x) {q_13(x,:)},1:dim), arrayfun(@(x) {p_13(x,:)},1:dim)];
z_14 = [arrayfun(@(x) {q_14(x,:)},1:dim), arrayfun(@(x) {p_14(x,:)},1:dim)];

H_Tao2_error         = (H(z_1{:}) - H_int)/H_int;
H_Strang_proj2_error = (H(z_2{:}) - H_int)/H_int;
H_semi2_error        = (H(z_3{:}) - H_int)/H_int;
H_mid_error          = (H(z_4{:}) - H_int)/H_int;
H_Tao4_error         = (H(z_5{:}) - H_int)/H_int;
H_Strang_proj4_error = (H(z_6{:}) - H_int)/H_int;
H_Triple4_error      = (H(z_7{:}) - H_int)/H_int;
H_Suzuki4_error      = (H(z_8{:}) - H_int)/H_int;
H_IRK4_error         = (H(z_9{:}) - H_int)/H_int;
H_Tao6_error         = (H(z_10{:}) - H_int)/H_int;
H_Strang_proj6_error = (H(z_11{:}) - H_int)/H_int;
H_Triple6_error      = (H(z_12{:}) - H_int)/H_int;
H_Suzuki6_error      = (H(z_13{:}) - H_int)/H_int;
H_Yoshida6_error     = (H(z_14{:}) - H_int)/H_int;



%--------------------------------------------------------------------------
% Figures
%--------------------------------------------------------------------------
time = 0:dt:dt*itr;


figure
set(0,'DefaultAxesFontSize', 12)

set(gcf,'Position',[300 250 600 600])

sgtitle('$\frac{H-H_{0}}{H_{0}}$','interpreter','latex','FontSize', 16)

subplot(3,2,1)
plot(time,H_Tao2_error); hold on
plot(time,H_Strang_proj2_error)
xlabel('time','interpreter','latex')
legend('Tao 2','Strang-proj 2','Location','best')


subplot(3,2,2)
plot(time,H_Strang_proj2_error); hold on
plot(time,H_semi2_error,'--')
plot(time,H_mid_error)
xlabel('time','Interpreter','latex')
legend('Strang-proj 2','semiexplicit 2','midpoint','Location','best')


subplot(3,2,3)
plot(time,H_Tao4_error); hold on
plot(time,H_Strang_proj4_error)
xlabel('time','interpreter','latex')
legend('Tao 4','Strang-proj 4','Location','best')


subplot(3,2,4)
plot(time,H_Strang_proj4_error) ; hold on
plot(time,H_Triple4_error,'--')
plot(time,H_IRK4_error)
xlabel('time','Interpreter','latex')
legend('Strang-proj 4','semiexplicit 4','IRK4','Location','best')


subplot(3,2,5)
plot(time,H_Tao6_error); hold on
plot(time,H_Strang_proj6_error)
xlabel('time','interpreter','latex')
legend('Tao 6','Strang-proj 6','Location','best')


subplot(3,2,6)
plot(time,H_Strang_proj6_error) ; hold on
plot(time,H_Triple6_error,'--')
xlabel('time','Interpreter','latex')
legend('Strang-proj 6','semiexplicit 6','Location','best')


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


set(gcf,'Position',[100 50 800 600])







