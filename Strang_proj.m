classdef Strang_proj

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
                    
                q(:,i+1) = (qB+xB)/2;
                p(:,i+1) = (pB+yB)/2;
                         
            end
            method_time = toc;
           
            
            
            
            disp('2nd order Strang_proj done')
            
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
                    
                q(:,i+1) = (qB+xB)/2;
                p(:,i+1) = (pB+yB)/2;
                      
            end
            method_time = toc;
            
            disp('4th order Strang_proj Triple Jump done')
            
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
                    
                q(:,i+1) = (qB+xB)/2;
                p(:,i+1) = (pB+yB)/2;
                      
            end
            method_time = toc;
            
            disp('6th order Strang_proj Triple Jump done')
            
            
        end
        
        
        function [q, p, method_time] = Su4(gradq, gradp, q_int, p_int, dt, itr)

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
            q = zeros(dim,itr+1);    
            p = zeros(dim,itr+1);

            q(:,1) = q_int;
            p(:,1) = p_int;
            
            
            %--------------------------------------------------------------------------
            % 4th order Strang_proj Triple jump
            %--------------------------------------------------------------------------
            tic;
            for i = 1:itr
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, q(:,i), p(:,i), q(:,i), p(:,i), SUgam_4(1)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*dt);
                    
                        
                q(:,i+1) = (qB+xB)/2;
                p(:,i+1) = (pB+yB)/2;
                      
            end
            method_time = toc;
            
            disp('4th order Strang_proj Suzuki done')
            
        end
        
        
        function [q, p, method_time] = Su6(gradq, gradp, q_int, p_int, dt, itr)

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
            q = zeros(dim,itr+1);    
            p = zeros(dim,itr+1);

            q(:,1) = q_int;
            p(:,1) = p_int;
            
            
            %--------------------------------------------------------------------------
            % 6th order Strang_proj Triple jump
            %--------------------------------------------------------------------------
            tic;
            for i=1:itr
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, q(:,i), p(:,i), q(:,i), p(:,i), SUgam_4(1)*SUgam_6(1)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*SUgam_6(1)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*SUgam_6(1)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*SUgam_6(1)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*SUgam_6(1)*dt);
                    
                    
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(1)*SUgam_6(2)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*SUgam_6(2)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*SUgam_6(2)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*SUgam_6(2)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*SUgam_6(2)*dt);
                    
                    
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(1)*SUgam_6(3)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*SUgam_6(3)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*SUgam_6(3)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*SUgam_6(3)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*SUgam_6(3)*dt);

                    
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(1)*SUgam_6(4)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*SUgam_6(4)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*SUgam_6(4)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*SUgam_6(4)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*SUgam_6(4)*dt);


                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(1)*SUgam_6(5)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(2)*SUgam_6(5)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(3)*SUgam_6(5)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(4)*SUgam_6(5)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, SUgam_4(5)*SUgam_6(5)*dt);
                    
                        
                q(:,i+1) = (qB+xB)/2;
                p(:,i+1) = (pB+yB)/2;
                      
            end
            method_time = toc;
            
            disp('6th order Strang_proj Suzuki done')
            
            
        end
        
        
        function [q, p, method_time] = Yo6(gradq, gradp, q_int, p_int, dt, itr)

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
            q = zeros(dim,itr+1);    
            p = zeros(dim,itr+1);

            q(:,1) = q_int;
            p(:,1) = p_int;
            
            
            %--------------------------------------------------------------------------
            % 6th order Strang_proj Triple jump
            %--------------------------------------------------------------------------
            tic;
            for i=1:itr
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, q(:,i), p(:,i), q(:,i), p(:,i), YOgam_6(1)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(2)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(3)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(4)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(5)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(6)*dt);
                [qB, pB, xB, yB] = Strang_proj.strang_split(gradq, gradp, qB, pB, xB, yB, YOgam_6(7)*dt);
                        
                        
                q(:,i+1) = (qB+xB)/2;
                p(:,i+1) = (pB+yB)/2;
                      
            end
            method_time = toc;
            
            disp('6th order Strang_proj Yoshida done')
            
            
        end
        
        
     
        
        
        
    end
    
end
















