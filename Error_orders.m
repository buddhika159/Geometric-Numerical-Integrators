clear all
close all
clc

omega = 20;
time = 100;
dt = 0.1./(1:5:100);
tol = 1e-15;
Max_NW_Iter = 100;
dim = 1;
name = 'Exact-problem';

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


z_int = num2cell([q_int;p_int]);

H_int = H(z_int{:});





for k=1:length(dt)
    itr = time/dt(k);
    
    q_1 = [];	q_2 = [];	q_3 = [];
    p_1 = [];	p_2 = [];	p_3 = [];
    
    q_4 = [];	q_5 = [];	q_6 = [];   
    p_4 = [];	p_5 = [];	p_6 = [];   
    
    q_7 = [];	q_8 = [];	q_9 = [];	q_10 = [];	q_11 = [];	q_12 = [];
    p_7 = [];	p_8 = [];	p_9 = [];	p_10 = [];	p_11 = [];	p_12 = [];
    
    q_13 = [];	q_14 = [];
    p_13 = [];	p_14 = [];
    
    q_1(:,1) = q_int;	q_2(:,1) = q_int;	q_3(:,1) = q_int;
    p_1(:,1) = p_int;	p_2(:,1) = p_int;	p_3(:,1) = p_int; 
    

    q_4(:,1) = q_int;	q_5(:,1) = q_int;	q_6(:,1) = q_int;    
    p_4(:,1) = p_int;	p_5(:,1) = p_int;	p_6(:,1) = p_int;    

    q_7(:,1) = q_int;	q_8(:,1) = q_int;	q_9(:,1) = q_int;	q_10(:,1) = q_int;	q_11(:,1) = q_int;	q_12(:,1) = q_int;
    p_7(:,1) = p_int;	p_8(:,1) = p_int;	p_9(:,1) = p_int;	p_10(:,1) = p_int;	p_11(:,1) = p_int;	p_12(:,1) = p_int;
    
    q_13(:,1) = q_int;	q_14(:,1) = q_int;
    p_13(:,1) = p_int;	p_14(:,1) = p_int;

    for i=1:itr


        [q_1_temp, p_1_temp, ~]         = Tao_method.base2  (gradq, gradp, q_1(:,i), p_1(:,i), dt(k), 1, omega);
        [q_2_temp, p_2_temp, ~]         = Tao_method.TJ4    (gradq, gradp, q_2(:,i), p_2(:,i), dt(k), 1, omega);
        [q_3_temp, p_3_temp, ~]         = Tao_method.TJ6    (gradq, gradp, q_3(:,i), p_3(:,i), dt(k), 1, omega);
        [q_4_temp, p_4_temp, ~]         = Strang_proj.base2 (gradq, gradp, q_4(:,i), p_4(:,i), dt(k), 1); 
        [q_5_temp, p_5_temp, ~]         = Strang_proj.TJ4   (gradq, gradp, q_5(:,i), p_5(:,i), dt(k), 1);   
        [q_6_temp, p_6_temp, ~]         = Strang_proj.TJ6   (gradq, gradp, q_6(:,i), p_6(:,i), dt(k), 1);
        [q_7_temp, p_7_temp, ~, ~, ~]   = semiexplicit.base2(gradq, gradp, q_7(:,i), p_7(:,i), dt(k), 1, tol, Max_NW_Iter); 
        [q_8_temp, p_8_temp, ~, ~, ~]   = semiexplicit.TJ4  (gradq, gradp, q_8(:,i), p_8(:,i), dt(k), 1, tol, Max_NW_Iter);   
        [q_9_temp, p_9_temp, ~, ~, ~]   = semiexplicit.TJ6  (gradq, gradp, q_9(:,i), p_9(:,i), dt(k), 1, tol, Max_NW_Iter);
        [q_10_temp, p_10_temp, ~, ~, ~] = semiexplicit.Su4  (gradq, gradp, q_10(:,i), p_10(:,i), dt(k), 1, tol, Max_NW_Iter);
        [q_11_temp, p_11_temp, ~, ~, ~] = semiexplicit.Su6  (gradq, gradp, q_11(:,i), p_11(:,i), dt(k), 1, tol, Max_NW_Iter);
        [q_12_temp, p_12_temp, ~, ~, ~] = semiexplicit.Yo6  (gradq, gradp, q_12(:,i), p_12(:,i), dt(k), 1, tol, Max_NW_Iter);   
        [q_13_temp, p_13_temp, ~, ~, ~] = IMid2             (gradq, gradp, gradqq, gradqp, gradpp, q_13(:,i), p_13(:,i), dt(k), 1, tol, Max_NW_Iter);   
        [q_14_temp, p_14_temp, ~, ~, ~] = IRK4              (gradq, gradp, gradqq, gradqp, gradpp, q_14(:,i), p_14(:,i), dt(k), 1, tol, Max_NW_Iter);   
        
        
        q_1(:,i+1) = q_1_temp(:,end);      p_1(:,i+1) = p_1_temp(:,end);
        q_2(:,i+1) = q_2_temp(:,end);      p_2(:,i+1) = p_2_temp(:,end);
        q_3(:,i+1) = q_3_temp(:,end);      p_3(:,i+1) = p_3_temp(:,end);
        q_4(:,i+1) = q_4_temp(:,end);      p_4(:,i+1) = p_4_temp(:,end);
        q_5(:,i+1) = q_5_temp(:,end);      p_5(:,i+1) = p_5_temp(:,end);
        q_6(:,i+1) = q_6_temp(:,end);      p_6(:,i+1) = p_6_temp(:,end);
        q_7(:,i+1) = q_7_temp(:,end);      p_7(:,i+1) = p_7_temp(:,end);
        q_8(:,i+1) = q_8_temp(:,end);      p_8(:,i+1) = p_8_temp(:,end);
        q_9(:,i+1) = q_9_temp(:,end);      p_9(:,i+1) = p_9_temp(:,end);
        q_10(:,i+1) = q_10_temp(:,end);    p_10(:,i+1) = p_10_temp(:,end);
        q_11(:,i+1) = q_11_temp(:,end);    p_11(:,i+1) = p_11_temp(:,end);
        q_12(:,i+1) = q_12_temp(:,end);    p_12(:,i+1) = p_12_temp(:,end);
        q_13(:,i+1) = q_13_temp(:,end);    p_13(:,i+1) = p_13_temp(:,end);
        q_14(:,i+1) = q_14_temp(:,end);    p_14(:,i+1) = p_14_temp(:,end);
        
        
        
        
        


    end

    H_Tao2_error(k)         = max(abs((H(q_1,p_1) - H_int)/H_int));
    H_Tao4_error(k)         = max(abs((H(q_2,p_2) - H_int)/H_int));
    H_Tao6_error(k)         = max(abs((H(q_3,p_3) - H_int)/H_int));
    H_Strang_proj2_error(k) = max(abs((H(q_4,p_4) - H_int)/H_int));
    H_Strang_proj4_error(k) = max(abs((H(q_5,p_5) - H_int)/H_int));
    H_Strang_proj6_error(k) = max(abs((H(q_6,p_6) - H_int)/H_int));
    H_semi2_error(k)        = max(abs((H(q_7,p_7) - H_int)/H_int));
    H_Triple4_error(k)      = max(abs((H(q_8,p_8) - H_int)/H_int));
    H_Triple6_error(k)      = max(abs((H(q_9,p_9) - H_int)/H_int));
    H_Suzuki4_error(k)      = max(abs((H(q_10,p_10) - H_int)/H_int));
    H_Suzuki6_error(k)      = max(abs((H(q_11,p_11) - H_int)/H_int));
    H_Yoshida6_error(k)     = max(abs((H(q_12,p_12) - H_int)/H_int));
    H_mid_error(k)          = max(abs((H(q_13,p_13) - H_int)/H_int));
    H_IRK4_error(k)         = max(abs((H(q_14,p_14) - H_int)/H_int));

    
    k
end





for j = 2:length(dt)
    slop1(j) = (log(H_Tao2_error(j-1))         - log(H_Tao2_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop2(j) = (log(H_Tao4_error(j-1))         - log(H_Tao4_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop3(j) = (log(H_Tao6_error(j-1))         - log(H_Tao6_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop4(j) = (log(H_Strang_proj2_error(j-1)) - log(H_Strang_proj2_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop5(j) = (log(H_Strang_proj4_error(j-1)) - log(H_Strang_proj4_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop6(j) = (log(H_Strang_proj6_error(j-1)) - log(H_Strang_proj6_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop7(j) = (log(H_semi2_error(j-1))        - log(H_semi2_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop8(j) = (log(H_Triple4_error(j-1))      - log(H_Triple4_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop9(j) = (log(H_Triple6_error(j-1))      - log(H_Triple6_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop10(j)= (log(H_Suzuki4_error(j-1))      - log(H_Suzuki4_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop11(j)= (log(H_Suzuki6_error(j-1))      - log(H_Suzuki6_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop12(j)= (log(H_Yoshida6_error(j-1))     - log(H_Yoshida6_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop13(j)= (log(H_mid_error(j-1))          - log(H_mid_error(j)))/(log(dt(j-1)) - log(dt(j)));
    slop14(j)= (log(H_IRK4_error(j-1))         - log(H_IRK4_error(j)))/(log(dt(j-1)) - log(dt(j)));
    
end

[dt' slop1' slop2' slop3' slop4' slop5' slop6' slop7' slop8' slop9' slop10' slop11' slop12' slop13' slop14']



figure(20)

set(gcf,'Position',[400   100   770   634])
set(0,'DefaultAxesFontSize', 12)

% mJet = jet(8);
% Markers = {'-+','-o','-*','-x','-v','-d','-^','-s','->','-<'};

loglog(dt,H_Tao2_error,        '-^','MarkerSize',8,'Color','#0072BD');hold on
loglog(dt,H_Strang_proj2_error,'-^','MarkerSize',8,'Color','#EDB120','LineWidth',2)
loglog(dt,H_semi2_error,       '-^','MarkerSize',8,'Color','#FF0000')
loglog(dt,H_mid_error,         '-^','MarkerSize',8,'Color','#FF00FF')

loglog(dt,H_Tao4_error,        '-*','MarkerSize',8,'Color','#0072BD')
loglog(dt,H_Strang_proj4_error,'-*','MarkerSize',8,'Color','#EDB120','LineWidth',2)
loglog(dt,H_Triple4_error,     '-*','MarkerSize',8,'Color','#FF0000')
loglog(dt,H_Suzuki4_error,     '-*','MarkerSize',8,'Color','#000000')
loglog(dt,H_IRK4_error,        '-*','MarkerSize',8,'Color','#A2142F')

loglog(dt,H_Tao6_error,        '-o','MarkerSize',8,'Color','#0072BD')
loglog(dt,H_Strang_proj6_error,'-o','MarkerSize',8,'Color','#EDB120','LineWidth',2)
loglog(dt,H_Triple6_error,     '-o','MarkerSize',8,'Color','#FF0000')
loglog(dt,H_Suzuki6_error,     '-o','MarkerSize',8,'Color','#000000')
loglog(dt,H_Yoshida6_error,    '-o','MarkerSize',8,'Color','#77AC30')

   
legend('Tao 2','Strang-proj 2','semiexplicit 2','Midpoint', ...
       'Tao 4','Strang-proj 4','semiexplicit 4','semiexplicit-S 4','IRK4', ...
       'Tao 6','Strang-proj 6','semiexplicit 6','semiexplicit-S 6','semiexplicit-Y 6');

xlabel('$\Delta t$','interpreter','latex','FontSize', 20)
ylabel('$\max\left|\frac{H-H_{0}}{H_{0}}\right|$','interpreter','latex','FontSize', 20)

% axis equal
set(gca, ...
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
    'LineWidth'   , 1         );


set(findall(gcf,'type','legend'), ...
    'FontSize'      , 12        , ...
    'Interpreter'   , 'latex'   , ...
    'Location'      , 'southeast'    )








% saveas(figure(10),[name,' tol=',num2str(tol),'.eps'])
% print -depsc exact_orders

% 
% 
% filename = ['C:\Users\banda\Dropbox\UT Dallas\modified\New folder (2)\Error orders\Exact_order_tol=',num2str(tol),'.mat'];
% save(filename,'dt','omega','tol','H_semi2_error','H_Suzuki4_error','H_Suzuki6_error','H_Yoshida6_error','H_Triple4_error','H_Triple6_error','H_mid_error','H_IRK4_error','H_Tao2_error','H_Tao4_error','H_Tao6_error')
% 
% 
% 


















