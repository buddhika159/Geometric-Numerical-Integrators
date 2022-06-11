classdef Strang

    methods (Static)
        function [qB, pB, xB, yB] = strang_split(gradq, gradp, qB, pB, xB, yB, dt)
            
            xA = xB + dt/2 * gradp(qB,yB);
            pA = pB - dt/2 * gradq(qB,yB);    
            qB = qB + dt * gradp(xA,pA);
            yB = yB - dt * gradq(xA,pA);
            xB = xA + dt/2 * gradp(qB,yB);
            pB = pA - dt/2 * gradq(qB,yB);

        end
        
        function [q, p, method_time] = base2(gradq, gradp, q_int, p_int, dt, itr)
            
            %--------------------------------------------------------------------------
            % Initialize
            %--------------------------------------------------------------------------
            
            dim = length(q_int);            
            q = zeros(dim,itr+1);    
            p = zeros(dim,itr+1);

            q(:,1) = q_int;
            p(:,1) = p_int;
            
            %--------------------------------------------------------------------------
            % 2nd order Strang_proj
            %--------------------------------------------------------------------------
            
            
            tic;
            for i = 1:itr
              
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, q(:,i), p(:,i), q(:,i), p(:,i), dt);
                    
                q(:,i+1) = qB;
                p(:,i+1) = pB;
                         
            end
            method_time = toc;
           
            
            
            
            disp('2nd order Strang done')
            
        end
        
        
        function [q, p, method_time] = TJ4(gradq, gradp, q_int, p_int, dt, itr)

            %--------------------------------------------------------------------------
            % Initialize
            %--------------------------------------------------------------------------
            
            l = 2;
            TJgam_4(1) = 1/(2-2^(1/(l+1)));
            TJgam_4(2) = 1 - 2*TJgam_4(1);
            TJgam_4(3) = 1/(2-2^(1/(l+1)));
            
            
            dim = length(q_int);            
            q = zeros(dim,itr+1);    
            p = zeros(dim,itr+1);

            q(:,1) = q_int;
            p(:,1) = p_int;
            
            
            %--------------------------------------------------------------------------
            % 4th order Strang_proj Triple jump
            %--------------------------------------------------------------------------
            tic;
            for i = 1:itr
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, q(:,i), p(:,i), q(:,i), p(:,i), TJgam_4(1)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(2)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(3)*dt);
                    
                q(:,i+1) = qB;
                p(:,i+1) = pB;
                      
            end
            method_time = toc;
            
            disp('4th order Strang Triple Jump done')
            
        end
        
        
        function [q, p, method_time] = TJ6(gradq, gradp, q_int, p_int, dt, itr)

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
            q = zeros(dim,itr+1);    
            p = zeros(dim,itr+1);

            q(:,1) = q_int;
            p(:,1) = p_int;
            
            
            %--------------------------------------------------------------------------
            % 6th order Strang_proj Triple jump
            %--------------------------------------------------------------------------
            tic;
            for i=1:itr
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, q(:,i), p(:,i), q(:,i), p(:,i), TJgam_4(1)*TJgam_6(1)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(2)*TJgam_6(1)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(3)*TJgam_6(1)*dt);
                    
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(1)*TJgam_6(2)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(2)*TJgam_6(2)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(3)*TJgam_6(2)*dt);
                    
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(1)*TJgam_6(3)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(2)*TJgam_6(3)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, TJgam_4(3)*TJgam_6(3)*dt);
                    
                q(:,i+1) = qB;
                p(:,i+1) = pB;
                      
            end
            method_time = toc;
            
            disp('6th order Strang Triple Jump done')
            
            
        end
        
        

        
        
        
        
        
        
        
        
    end
    
end
















