classdef Tao_method
    
    methods (Static)
        
        function [qA, pA, xA, yA] = Tao_split(gradq, gradp, qA, pA, xA, yA, dt, omega)
       
            xA = xA + dt/2 * gradp(qA,yA);
            pA = pA - dt/2 * gradq(qA,yA);
            qB = qA + dt/2 * gradp(xA,pA);
            yB = yA - dt/2 * gradq(xA,pA);
            
            temp1 = (qB + xA + (qB - xA)*cos(2*omega*dt) + (pA - yB)*sin(2*omega*dt))/2;
            temp2 = (pA + yB + (pA - yB)*cos(2*omega*dt) - (qB - xA)*sin(2*omega*dt))/2;
            temp3 = (qB + xA - (qB - xA)*cos(2*omega*dt) - (pA - yB)*sin(2*omega*dt))/2;
            temp4 = (pA + yB - (pA - yB)*cos(2*omega*dt) + (qB - xA)*sin(2*omega*dt))/2;
            
            qA = temp1 + dt/2 * gradp(temp3,temp2);
            yA = temp4 - dt/2 * gradq(temp3,temp2);
            xA = temp3 + dt/2 * gradp(qA,yA);
            pA = temp2 - dt/2 * gradq(qA,yA);
                  
        end
       
        function [q, p, Max_Defect, method_time] = base2(gradq, gradp, q_int, p_int, dt, itr, omega)
            % INPUTS ------------------------------------------------------------------
            % gradq - Funciton handle for the gradiant vactor of the Hamiltonian w.r.t. q
            % gradp - Funciton handle for the gradiant vactor of the Hamiltonian w.r.t. p
            % q_int - Initial guess for q
            % p_int - Initial guess for p
            % dt    - step size
            % itr   - Number of iterations
            % omega - parameter of Tao's method
            
            % OUTPUTS -----------------------------------------------------------------
            % q - Approximates for q for the time frame 0:dt:dt*itr
            % p - Approximates for p for the time frame 0:dt:dt*itr
            
            
            
            
            
            %--------------------------------------------------------------------------
            % Initialize
            %--------------------------------------------------------------------------
            
            dim = length(q_int);
            
            q = zeros(dim,itr+1);
            p = zeros(dim,itr+1);
            x = zeros(dim,itr+1);
            y = zeros(dim,itr+1);
            
            Defects = zeros(itr,1);
            
            q(:,1) = q_int;
            p(:,1) = p_int;
            x(:,1) = q_int;
            y(:,1) = p_int;
            
            
            %--------------------------------------------------------------------------
            % 2nd order Tao
            %--------------------------------------------------------------------------
            tic
            for i=1:itr
                
                [q(:,i+1), p(:,i+1), x(:,i+1), y(:,i+1)] = Tao_method.Tao_split(gradq, gradp, q(:,i), p(:,i), x(:,i), y(:,i), dt, omega);
                
                Defects(i) = norm([q(:,i+1);p(:,i+1)] - [x(:,i+1);y(:,i+1)]);
                
            end
            method_time = toc;
            
            Max_Defect = max(Defects);
            
            disp('2nd order Tao done')
            
        end
        
        function [q, p, Max_Defect, method_time] = TJ4(gradq, gradp, q_int, p_int, dt, itr, omega)
            % INPUTS ------------------------------------------------------------------
            % gradq - Funciton handle for the gradiant vactor of the Hamiltonian w.r.t. q
            % gradp - Funciton handle for the gradiant vactor of the Hamiltonian w.r.t. p
            % q_int - Initial guess for q
            % p_int - Initial guess for p
            % dt    - step size
            % itr   - Number of iterations
            % omega - parameter of Tao's method
            
            % OUTPUTS -----------------------------------------------------------------
            % q - Approximates for q for the time frame 0:dt:dt*itr
            % p - Approximates for p for the time frame 0:dt:dt*itr
            
            %--------------------------------------------------------------------------
            % Initialize
            %--------------------------------------------------------------------------
            
            dim = length(q_int);
            
            
            l = 2;
            TJgam_4(1) = 1/(2-2^(1/(l+1)));
            TJgam_4(2) = 1 - 2*TJgam_4(1);
            TJgam_4(3) = 1/(2-2^(1/(l+1)));
            
            
            q = zeros(dim,itr+1);
            p = zeros(dim,itr+1);
            x = zeros(dim,itr+1);
            y = zeros(dim,itr+1);
            
            Defects = zeros(itr,1);
            
            q(:,1) = q_int;
            p(:,1) = p_int;
            x(:,1) = q_int;
            y(:,1) = p_int;
            
            
            %--------------------------------------------------------------------------
            % 4th order Tao Triple Jump
            %--------------------------------------------------------------------------
            
            tic;
            for i=1:itr
                
                [qA, pA, xA, yA] = Tao_method.Tao_split(gradq, gradp, q(:,i), p(:,i), x(:,i), y(:,i), TJgam_4(1)*dt, omega);
                [qA, pA, xA, yA] = Tao_method.Tao_split(gradq, gradp, qA, pA, xA, yA, TJgam_4(2)*dt, omega);
                [q(:,i+1), p(:,i+1), x(:,i+1), y(:,i+1)] = Tao_method.Tao_split(gradq, gradp, qA, pA, xA, yA, TJgam_4(3)*dt, omega);
                
                Defects(i) = norm([q(:,i+1);p(:,i+1)] - [x(:,i+1);y(:,i+1)]);
                
            end
            method_time = toc;
            
            Max_Defect = max(Defects);
            
            disp('4th order Tao Triple Jump done')
            
            
        end
        
        function [q, p, Max_Defect, method_time] = TJ6(gradq, gradp, q_int, p_int, dt, itr, omega)
            % INPUTS ------------------------------------------------------------------
            % gradq - Funciton handle for the gradiant vactor of the Hamiltonian w.r.t. q
            % gradp - Funciton handle for the gradiant vactor of the Hamiltonian w.r.t. p
            % q_int - Initial guess for q
            % p_int - Initial guess for p
            % dt    - step size
            % itr   - Number of iterations
            % omega - parameter of Tao's method
            
            % OUTPUTS -----------------------------------------------------------------
            % q - Approximates for q for the time frame 0:dt:dt*itr
            % p - Approximates for p for the time frame 0:dt:dt*itr
            
            %--------------------------------------------------------------------------
            % Initialize
            %--------------------------------------------------------------------------
            
            dim = length(q_int);
            
            
            pp = 2;
            TJgam_4(1) = 1/(2-2^(1/(pp+1)));
            TJgam_4(2) = 1 - 2*TJgam_4(1);
            TJgam_4(3) = 1/(2-2^(1/(pp+1)));
            
            pp = 4;
            TJgam_6(1) = 1/(2-2^(1/(pp+1)));
            TJgam_6(2) = 1 - 2*TJgam_6(1);
            TJgam_6(3) = 1/(2-2^(1/(pp+1)));
            
            
            q = zeros(dim,itr+1);
            p = zeros(dim,itr+1);
            x = zeros(dim,itr+1);
            y = zeros(dim,itr+1);
            
            Defects = zeros(itr,1);
            
            q(:,1) = q_int;
            p(:,1) = p_int;
            x(:,1) = q_int;
            y(:,1) = p_int;
            
            %--------------------------------------------------------------------------
            % 6th order Tao Triple Jump
            %--------------------------------------------------------------------------
            tic;
            for i=1:itr
                
                [qA, pA, xA, yA] = Tao_method.Tao_split(gradq, gradp, q(:,i), p(:,i), x(:,i), y(:,i), TJgam_4(1)*TJgam_6(1)*dt, omega);
                [qA, pA, xA, yA] = Tao_method.Tao_split(gradq, gradp, qA, pA, xA, yA, TJgam_4(2)*TJgam_6(1)*dt, omega);
                [qA, pA, xA, yA] = Tao_method.Tao_split(gradq, gradp, qA, pA, xA, yA, TJgam_4(3)*TJgam_6(1)*dt, omega);
                
                
                [qA, pA, xA, yA] = Tao_method.Tao_split(gradq, gradp, qA, pA, xA, yA, TJgam_4(1)*TJgam_6(2)*dt, omega);
                [qA, pA, xA, yA] = Tao_method.Tao_split(gradq, gradp, qA, pA, xA, yA, TJgam_4(2)*TJgam_6(2)*dt, omega);
                [qA, pA, xA, yA] = Tao_method.Tao_split(gradq, gradp, qA, pA, xA, yA, TJgam_4(3)*TJgam_6(2)*dt, omega);
                
                
                [qA, pA, xA, yA] = Tao_method.Tao_split(gradq, gradp, qA, pA, xA, yA, TJgam_4(1)*TJgam_6(3)*dt, omega);
                [qA, pA, xA, yA] = Tao_method.Tao_split(gradq, gradp, qA, pA, xA, yA, TJgam_4(2)*TJgam_6(3)*dt, omega);
                [q(:,i+1), p(:,i+1), x(:,i+1), y(:,i+1)] = Tao_method.Tao_split(gradq, gradp, qA, pA, xA, yA, TJgam_4(3)*TJgam_6(3)*dt, omega);
                
                Defects(i) = norm([q(:,i+1);p(:,i+1)] - [x(:,i+1);y(:,i+1)]);
                
            end
            method_time = toc;
            
            Max_Defect = max(Defects);        
            
            disp('6th order Tao Triple Jump done')
            
            
        end
        
        
        
        
        
        
        
        
        
        
        
    end
    
end
















