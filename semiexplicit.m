classdef semiexplicit

    methods (Static)
        function [qB, pB, xB, yB] = strang_split(gradq, gradp, qB, pB, xB, yB, dt)
            
            xA = xB + dt/2 * gradp(qB,yB);
            pA = pB - dt/2 * gradq(qB,yB);    
            qB = qB + dt * gradp(xA,pA);
            yB = yB - dt * gradq(xA,pA);
            xB = xA + dt/2 * gradp(qB,yB);
            pB = pA - dt/2 * gradq(qB,yB);

        end
        
        function [q, p, NW_itrs, Max_Error, Max_Defect, method_time] = base2(gradq, gradp, q_int, p_int, dt, itr, tol, Max_NW_Iter)
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
            z = zeros(4*dim,itr+1);               
            Defects = zeros(itr,1);
            Errors = zeros(itr,1);
            NW_itrs = 0; 
            
            z(:,1) = [q_int;p_int;q_int;p_int];
            
            A =@(z) [z(0*dim+1:1*dim) - z(2*dim+1:3*dim)
                     z(1*dim+1:2*dim) - z(3*dim+1:4*dim)];
            
            
            %--------------------------------------------------------------------------
            % 2nd order semiExplicit
            %--------------------------------------------------------------------------
            tic;
            for i = 1:itr
                mu = zeros(2*dim,1);
                z_int = z(:,i);
                
                error = 1;
                NW_i = 0;
                while (error > tol && NW_i < Max_NW_Iter)
                    z_in = z_int + [mu;-mu];
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, z_in(0*dim+1:1*dim), z_in(1*dim+1:2*dim), z_in(2*dim+1:3*dim), z_in(3*dim+1:4*dim), dt);
                    
                    proj = [qB;pB;xB;yB] + [mu;-mu];
                    
                    % Newton itterations
                    
                    A_val = A(proj);
                    mu = mu - A_val/4;
                    
                    error = norm(-A_val/4);
                    NW_i = NW_i + 1;
                end
                
                z(:,i+1) = proj;                 
                Defects(i) = norm(z(0*dim+1:2*dim,i+1) - z(2*dim+1:4*dim,i+1)); 
                Errors(i) = error; 
                NW_itrs = NW_itrs + NW_i;
                
                
            end
            method_time = toc;
            
            q = z(0*dim+1:1*dim,:);
            p = z(1*dim+1:2*dim,:);                                                    
            
            Max_Defect = max(Defects);
            Max_Error = max(Errors);             
            
            
            
            disp('2nd order semixplicit done')
            
        end
        
        
        function [q, p, NW_itrs, Max_Error, Max_Defect, method_time] = TJ4(gradq, gradp, q_int, p_int, dt, itr, tol, Max_NW_Iter)
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
            
            l = 2;
            TJgam_4(1) = 1/(2-2^(1/(l+1)));
            TJgam_4(2) = 1 - 2*TJgam_4(1);
            TJgam_4(3) = 1/(2-2^(1/(l+1)));
            
            
            dim = length(q_int);
            
            z = zeros(4*dim,itr+1);               
            Defects = zeros(itr,1); 
            Errors = zeros(itr,1);            
            NW_itrs = 0;             
            
            z(:,1) = [q_int;p_int;q_int;p_int];
            
            A =@(z) [z(0*dim+1:1*dim) - z(2*dim+1:3*dim)
                     z(1*dim+1:2*dim) - z(3*dim+1:4*dim)];
            
            
            %--------------------------------------------------------------------------
            % 4th order semiExplicit Triple jump
            %--------------------------------------------------------------------------
            tic;
            for i = 1:itr
                mu = zeros(2*dim,1);
                z_int = z(:,i);
                
                error = 1;
                NW_i = 0;
                while (error > tol && NW_i < Max_NW_Iter)
                    z_in = z_int + [mu;-mu];

                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, z_in(0*dim+1:1*dim), z_in(1*dim+1:2*dim), z_in(2*dim+1:3*dim), z_in(3*dim+1:4*dim), TJgam_4(1)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(2)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(3)*dt);
                    
                    proj = [qB;pB;xB;yB] + [mu;-mu];
                    
                    % Newton itterations
                    
                    A_val = A(proj);
                    mu = mu - A_val/4;
                    
                    error = norm(-A_val/4);
                    NW_i = NW_i + 1;
                end
                
                z(:,i+1) = proj;                 
                Defects(i) = norm(z(0*dim+1:2*dim,i+1) - z(2*dim+1:4*dim,i+1));    
                Errors(i) = error;             
                NW_itrs = NW_itrs + NW_i;
                
                
            end
            method_time = toc;
            
            q = z(0*dim+1:1*dim,:);
            p = z(1*dim+1:2*dim,:);                                                    
            
            Max_Defect = max(Defects);
            Max_Error = max(Errors);             
            
            
            disp('4th order semiexplicit Triple Jump done')
            
        end
        
        
        function [q, p, NW_itrs, Max_Error, Max_Defect, method_time] = TJ6(gradq, gradp, q_int, p_int, dt, itr, tol, Max_NW_Iter)
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
            
            l = 2;
            TJgam_4(1) = 1/(2-2^(1/(l+1)));
            TJgam_4(2) = 1 - 2*TJgam_4(1);
            TJgam_4(3) = 1/(2-2^(1/(l+1)));
            
            l = 4;
            TJgam_6(1) = 1/(2-2^(1/(l+1)));
            TJgam_6(2) = 1 - 2*TJgam_6(1);
            TJgam_6(3) = 1/(2-2^(1/(l+1)));
            
            
            dim = length(q_int);
            
            z = zeros(4*dim,itr+1);               
            Defects = zeros(itr,1); 
            Errors = zeros(itr,1);            
            NW_itrs = 0;             
            
            z(:,1) = [q_int;p_int;q_int;p_int];
            
            A =@(z) [z(0*dim+1:1*dim) - z(2*dim+1:3*dim)
                     z(1*dim+1:2*dim) - z(3*dim+1:4*dim)];
            
            
            %--------------------------------------------------------------------------
            % 6th order semiexplicit Triple jump
            %--------------------------------------------------------------------------
            tic;
            for i=1:itr
                mu = zeros(2*dim,1);
                z_int = z(:,i);
                
                error = 1;
                NW_i = 0;
                while (error > tol && NW_i < Max_NW_Iter)
                    z_in = z_int + [mu;-mu];
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, z_in(0*dim+1:1*dim), z_in(1*dim+1:2*dim), z_in(2*dim+1:3*dim), z_in(3*dim+1:4*dim), TJgam_4(1)*TJgam_6(1)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(2)*TJgam_6(1)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(3)*TJgam_6(1)*dt);
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(1)*TJgam_6(2)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(2)*TJgam_6(2)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(3)*TJgam_6(2)*dt);
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(1)*TJgam_6(3)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(2)*TJgam_6(3)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(3)*TJgam_6(3)*dt);
                    
                    
                    proj = [qB;pB;xB;yB] + [mu;-mu];
                    
                    % Newton itterations
                    
                    A_val = A(proj);
                    mu = mu - A_val/4;
                    
                    error = norm(-A_val/4);
                    NW_i = NW_i + 1;
                end
                
                z(:,i+1) = proj;                 
                Defects(i) = norm(z(0*dim+1:2*dim,i+1) - z(2*dim+1:4*dim,i+1));   
                Errors(i) = error;              
                NW_itrs = NW_itrs + NW_i;
                
                
            end
            method_time = toc;
            
            q = z(0*dim+1:1*dim,:);
            p = z(1*dim+1:2*dim,:);                                                    
            
            Max_Defect = max(Defects);
            Max_Error = max(Errors);             
            
            
            disp('6th order semiexplicit Triple Jump done')
            
            
        end
          
        
        function [q, p, NW_itrs, Max_Error, Max_Defect, method_time] = Su4(gradq, gradp, q_int, p_int, dt, itr, tol, Max_NW_Iter)
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
            
            l = 2;
            SUgam_4(1) = 1/(4-4^(1/(l+1)));
            SUgam_4(2) = SUgam_4(1);
            SUgam_4(3) = 1 - 4*SUgam_4(1);
            SUgam_4(4) = SUgam_4(1);
            SUgam_4(5) = SUgam_4(1);
            
            dim = length(q_int);
            
            z = zeros(4*dim,itr+1);               
            Defects = zeros(itr,1);
            Errors = zeros(itr,1);             
            NW_itrs = 0;             
            
            z(:,1) = [q_int;p_int;q_int;p_int];
            
            A =@(z) [z(0*dim+1:1*dim) - z(2*dim+1:3*dim)
                     z(1*dim+1:2*dim) - z(3*dim+1:4*dim)];
            
            
            %--------------------------------------------------------------------------
            % 4th order semiexplicit Suzuki
            %--------------------------------------------------------------------------
            tic;
            for i=1:itr
                mu = zeros(2*dim,1);
                z_int = z(:,i);
                
                error = 1;
                NW_i = 0;
                while (error > tol && NW_i < Max_NW_Iter)
                    z_in = z_int + [mu;-mu];
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, z_in(0*dim+1:1*dim), z_in(1*dim+1:2*dim), z_in(2*dim+1:3*dim), z_in(3*dim+1:4*dim), SUgam_4(1)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*dt);
                    
                    
                    proj = [qB;pB;xB;yB] + [mu;-mu];
                    
                    % Newton itterations
                    
                    A_val = A(proj);
                    mu = mu - A_val/4;
                    
                    error = norm(-A_val/4);
                    NW_i = NW_i + 1;
                end
                
                z(:,i+1) = proj;                 
                Defects(i) = norm(z(0*dim+1:2*dim,i+1) - z(2*dim+1:4*dim,i+1)); 
                Errors(i) = error;                
                NW_itrs = NW_itrs + NW_i;
                
                
            end
            method_time = toc;
            
            q = z(0*dim+1:1*dim,:);
            p = z(1*dim+1:2*dim,:);                                                    
            
            Max_Defect = max(Defects);  
            Max_Error = max(Errors);           
            
            
            disp('4th order semiexplicit Suzuki done')
            
        end
        
        
        function [q, p, NW_itrs, Max_Error, Max_Defect, method_time] = Su6(gradq, gradp, q_int, p_int, dt, itr, tol, Max_NW_Iter)
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
            
            l = 2;
            SUgam_4(1) = 1/(4-4^(1/(l+1)));
            SUgam_4(2) = SUgam_4(1);
            SUgam_4(3) = 1 - 4*SUgam_4(1);
            SUgam_4(4) = SUgam_4(1);
            SUgam_4(5) = SUgam_4(1);
            
            l = 4;
            SUgam_6(1) = 1/(4-4^(1/(l+1)));
            SUgam_6(2) = SUgam_6(1);
            SUgam_6(3) = 1 - 4*SUgam_6(1);
            SUgam_6(4) = SUgam_6(1);
            SUgam_6(5) = SUgam_6(1);
            
            dim = length(q_int);
            
            z = zeros(4*dim,itr+1);               
            Defects = zeros(itr,1); 
            Errors = zeros(itr,1);            
            NW_itrs = 0;             
            
            z(:,1) = [q_int;p_int;q_int;p_int];
            
            A =@(z) [z(0*dim+1:1*dim) - z(2*dim+1:3*dim)
                     z(1*dim+1:2*dim) - z(3*dim+1:4*dim)];
            
            
            %--------------------------------------------------------------------------
            % 6th order semiexplicit Suzuki
            %--------------------------------------------------------------------------
            tic;
            for i=1:itr
                mu = zeros(2*dim,1);
                z_int = z(:,i);
                
                error = 1;
                NW_i = 0;
                while (error > tol && NW_i < Max_NW_Iter)
                    z_in = z_int + [mu;-mu];
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, z_in(0*dim+1:1*dim), z_in(1*dim+1:2*dim), z_in(2*dim+1:3*dim), z_in(3*dim+1:4*dim), SUgam_4(1)*SUgam_6(1)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*SUgam_6(1)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*SUgam_6(1)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*SUgam_6(1)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*SUgam_6(1)*dt);
                    
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(1)*SUgam_6(2)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*SUgam_6(2)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*SUgam_6(2)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*SUgam_6(2)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*SUgam_6(2)*dt);
                    
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(1)*SUgam_6(3)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*SUgam_6(3)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*SUgam_6(3)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*SUgam_6(3)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*SUgam_6(3)*dt);
                    
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(1)*SUgam_6(4)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*SUgam_6(4)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*SUgam_6(4)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*SUgam_6(4)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*SUgam_6(4)*dt);
                    
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(1)*SUgam_6(5)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*SUgam_6(5)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*SUgam_6(5)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*SUgam_6(5)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*SUgam_6(5)*dt);
                    
                    
                    
                    
                    proj = [qB;pB;xB;yB] + [mu;-mu];
                    
                    % Newton itterations
                    
                    A_val = A(proj);
                    mu = mu - A_val/4;
                    
                    error = norm(-A_val/4);
                    NW_i = NW_i + 1;
                end
                
                z(:,i+1) = proj;                 
                Defects(i) = norm(z(0*dim+1:2*dim,i+1) - z(2*dim+1:4*dim,i+1));  
                Errors(i) = error;               
                NW_itrs = NW_itrs + NW_i;
                
                
            end
            method_time = toc;
            
            q = z(0*dim+1:1*dim,:);
            p = z(1*dim+1:2*dim,:);                                                    
            
            Max_Defect = max(Defects);
            Max_Error = max(Errors);             
            
            
            disp('6th order semiexplicit Suzuki done')
            
        end
        
        
        function [q, p, NW_itrs, Max_Error, Max_Defect, method_time] = Yo6(gradq, gradp, q_int, p_int, dt, itr, tol, Max_NW_Iter)
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
            
            YOgam_6(1) = 0.78451361047755726381949763;
            YOgam_6(2) = 0.23557321335935813368479318;
            YOgam_6(3) = -1.17767998417887100694641568;
            YOgam_6(4) = 1.31518632068391121888424973;
            YOgam_6(5) = YOgam_6(3);
            YOgam_6(6) = YOgam_6(2);
            YOgam_6(7) = YOgam_6(1);
            
            
            dim = length(q_int);
            
            z = zeros(4*dim,itr+1);               
            Defects = zeros(itr,1);
            Errors = zeros(itr,1);             
            NW_itrs = 0;             
            
            z(:,1) = [q_int;p_int;q_int;p_int];
            
            A =@(z) [z(0*dim+1:1*dim) - z(2*dim+1:3*dim)
                     z(1*dim+1:2*dim) - z(3*dim+1:4*dim)];
            
            
            %--------------------------------------------------------------------------
            % 6th order semiexplicit Yoshida
            %--------------------------------------------------------------------------
            tic;
            for i=1:itr
                mu = zeros(2*dim,1);
                z_int = z(:,i);
                
                error = 1;
                NW_i = 0;
                while (error > tol && NW_i < Max_NW_Iter)
                    z_in = z_int + [mu;-mu];
                    
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, z_in(0*dim+1:1*dim), z_in(1*dim+1:2*dim), z_in(2*dim+1:3*dim), z_in(3*dim+1:4*dim), YOgam_6(1)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(2)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(3)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(4)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(5)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(6)*dt);
                    [qB, pB, xB, yB] = semiexplicit.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(7)*dt);
                    
                    
                    proj = [qB;pB;xB;yB] + [mu;-mu];
                    
                    % Newton itterations
                    
                    A_val = A(proj);
                    mu = mu - A_val/4;
                    
                    error = norm(-A_val/4);
                    NW_i = NW_i + 1;
                end
                
                z(:,i+1) = proj;                 
                Defects(i) = norm(z(0*dim+1:2*dim,i+1) - z(2*dim+1:4*dim,i+1));   
                Errors(i) = error;              
                NW_itrs = NW_itrs + NW_i;
                
                
            end
            method_time = toc;
            
            q = z(0*dim+1:1*dim,:);
            p = z(1*dim+1:2*dim,:);                                                    
            
            Max_Defect = max(Defects);
            Max_Error = max(Errors);             
            
            
            disp('6th order semiexplicit Yoshida done')
            
        end
        
        
        
        
        
        
        
    end
    
end
















