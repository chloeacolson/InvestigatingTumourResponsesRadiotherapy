%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: - Example code to simulate different RT treatment protocols and 
%           plot the numerical solution
%           - Assume that the file in location 'file_path' contains the
%           parameters values that define the tumours of interest, i.e., 
%           the oxygen consumption rates, q1 and q3, the vascular volume,
%           V0, and the steady state tumour volume and oxygen
%           concentration,T0 and c0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% Import tumour information
    tumours = readtable('file_path');

% Define model parameters

    % Carrying capacity 
        K = 1;
    
    % Anoxic O2 threshold
        c_min = 1e-2;
    
    % Rate of oxygen release from vasculature
        g = 5;

    % proliferation rate of undamaged tumour cells
        k1 = 1e-2;
        q1 = tumours.q1;
        q3 = tumours.q3;
        q2 = k1*q3;

    % Cell oxygen consumption and proliferation rates
        q3sl = 1e-1*q3;
        q2sl = k1*q3sl;
        q1sl = 1e1*q1;

    % Cell death rates
        d1 = q2;
        d1sl = q2sl;

    % Number of RT fractions that correspond to the RT doses
    % 0, 1, 2, 3, 4, 5 Gy, respectively
        num_frac_RT = [0,80,40,26,20,16];
    
    % Time between fractions for 5x, 3x, 1x RT per week, respectively
    % 1 day = 1440 min, 2 days = 2880 min, 5 days = 7200 min
        time_btw_frac_RT = [1.44e3,2.88e3,7.2e3];
        
    % Duration of irradiation
        time_treat = 10;    
    
    % RT dose rates that correspond to the RT doses 0, 1, 2, 3, 4, 5 Gy,
    % respectively
        R = linspace(0,5,6)/(time_treat); 

    % Rates of sublethal (l1) and lethal (l2, l2sl) RT damage
        l1 = 10;
        l2 = 1;
        l2sl = 1;

    % RT damage repair rate
        mu = 5e-3;
        
    % Rate of mitotic catastrophe
        zeta = 5e-4;

    % Rate of dead cell clearance from tumour
        eta = 5e-5;

    % Solver options
        options = odeset('AbsTol',1e-10, 'RelTol',1e-10);

% Simulating a conventional RT schedule (5 x 2 Gy for 8 weeks)   

    % ID of the tumour we want to simulate treatment for
        j = 1;
        
    % RT dose rate/Number of RT fractions
        k = 3;
        
    % Number of fractions per week
        l = 1;

    % Define the initial conditions
        tt   = 0;
        TT_1 = tumours.T0(j);
        TT_2 = 0;
        TT_3 = 0;
        TT_4 = tumours.c0(j);
        V0   = tumours.V0(j);
    
    % Count the RT fraction number
        i = 1;
                
        while i < num_frac_RT(k)


            % Irradiation

            % Update time span
            t_1 = linspace(tt(end),tt(end) + time_treat,20); 

            % Update initial conditions
            IC = [TT_1(end), TT_2(end),TT_3(end),TT_4(end)];             

            % Solve the ODE model
            [tt_1,TT] = ode15s(@(t,T) RT(t, T, K, V0, q2(j), q1(j), q3(j),...
                                         g, q2sl(j), q1sl(j), q3sl(j),...
                                         l1, l2, l2sl, mu, eta, zeta,...
                                         d1(j), d1sl(j), c_min, R(k)),...
                                         t_1, IC, options);                  

            % Determing the break between the this fraction and the next
            if (mod(i,5) == 0)
                time_break = 2.88e3;
            else
                time_break = 0;
            end

            % Tumour growth between fractions

            % Update time span
            t_2 = linspace(tt_1(end),...
                tt_1(end) + time_btw_frac_RT(l) - time_treat + time_break, 500); 

            % Update initial condition
            IC_2 = [TT(end,1), TT(end,2), TT(end,3),TT(end,4)];  

            % Solve the model
            [tt_2,TTb] = ode45(@(t,T) RT(t, T, K, V0, q2(j), q1(j), q3(j),...
                                         g, q2sl(j), q1sl(j), q3sl(j),...
                                         l1, l2, l2sl, mu, eta, zeta,...
                                         d1(j), d1sl(j), c_min, R(1)),...
                                         t_2, IC_2, options); 

            % Appending the solutions for the simulated time periods
            tt = vertcat(tt,[tt_1;tt_2]);
            TT_1 = vertcat(TT_1, [TT(:,1);TTb(:,1)]);
            TT_2 = vertcat(TT_2, [TT(:,2);TTb(:,2)]);
            TT_3 = vertcat(TT_3, [TT(:,3);TTb(:,3)]);
            TT_4 = vertcat(TT_4, [TT(:,4);TTb(:,4)]);

            % Moving to the next fraction
            i = i+1;


        end

        % Final fraction

            % Irradiation

            % Update time span
            t_1 = linspace(tt(end),tt(end) + time_treat, 20);  

            % Update initial conditions
            IC = [TT_1(end), TT_2(end), TT_3(end),TT_4(end)];   

            % Solve the model
            [tt_1,TT] = ode15s(@(t,T) RT(t, T, K, V0, q2(j), q1(j), q3(j),...
                                         g, q2sl(j), q1sl(j), q3sl(j),...
                                         l1, l2, l2sl, mu, eta, zeta,...
                                         d1(j), d1sl(j), c_min, R(k)),...
                                         t_1, IC, options);

            % Post-RT tumour growth
            time_break_1 = 2.88e3; 
            time_break_2 = 1e6;    

            % Update time span
            t_2 = linspace(tt_1(end),...
                tt_1(end) + time_btw_frac_RT(l) - time_treat + time_break_1, 500); 

            % Update initial conditions
            IC_2 = [TT(end,1), TT(end,2), TT(end,3),TT(end,4)];                                        

            % Solve the ODE model
            [tt_2,TTb] = ode45(@(t,T) RT(t, T, K, V0, q2(j), q1(j), q3(j),...
                                         g, q2sl(j), q1sl(j), q3sl(j),...
                                         l1, l2, l2sl, mu, eta, zeta,...
                                         d1(j), d1sl(j), c_min, R(1)),...
                                         t_2, IC_2, options); 

            % Appending the simulated time periods
            tt = vertcat(tt,[tt_1;tt_2]);
            TT_1 = vertcat(TT_1, [TT(:,1);TTb(:,1)]);
            TT_2 = vertcat(TT_2, [TT(:,2);TTb(:,2)]);
            TT_3 = vertcat(TT_3, [TT(:,3);TTb(:,3)]);
            TT_4 = vertcat(TT_4, [TT(:,4);TTb(:,4)]);


            % Post-treatment tumour growth (comment this out for short-term
            % response only)

            % Update time span
            t_3 = linspace(tt(end),tt(end) + time_break_2, 500); 

            % Update initial conditions
            IC_3 = [TT_1(end), TT_2(end), TT_3(end),TT_4(end)];  

            % Solve the ODE model
            [tt_3,TTc] = ode45(@(t,T) RT(t, T, K, V0, q2(j), q1(j), q3(j),...
                                         g, q2sl(j), q1sl(j), q3sl(j),...
                                         l1, l2, l2sl, mu, eta, zeta,...
                                         d1(j), d1sl(j), c_min, R(1)),...
                                         t_3, IC_3, options); 

            % Appending the simulated time periods
            tt = vertcat(tt,tt_3);
            TT_1 = vertcat(TT_1, TTc(:,1));
            TT_2 = vertcat(TT_2, TTc(:,2));
            TT_3 = vertcat(TT_3, TTc(:,3));
            TT_4 = vertcat(TT_4, TTc(:,4));
   
% Plotting the numerical solution

    % Define figure
        p = figure
        hold on
    
    % Set axes properties
        set(gca,'FontSize',20,'FontName','Helvetica')
    
    % Plot the viable tumour cell volume
        plot(tt, TT_1+TT_2, 'LineWidth', 2.5,'color',"#00BFC4")
        
    % Plot the dead cell volume
        plot(tt, TT_3, 'LineWidth', 2.5,'color',"#C77CFF")
        
    % Plot the oxygen concentration
        plot(tt, TT_4, 'LineWidth', 2.5,'color',"#F8766D")
    
    % Legend
        legend('T+T_S','T_R','c','FontSize',16,...
           'FontName','Helvetica','location','east')
       
    % Define axis bounds (here set to the duration of the RT fractionation 
    % schedule)
        axis([0 80640 0 0.5])
    
    % x-axis label
        xlabel('time, t','FontSize',20,'FontName','Helvetica')
        
% Defining the ODE system
    function radiotherapy = RT(t, T, K, V_0, q_2, q_1, q_3, g, q_2sl, q_1sl,...
                               q_3sl, l1, l2, l2sl, mu, eta, zeta, d1, d1sl,...
                               cmin, R)

        DT = zeros(4,1);

        % Total tumour volume
        E = T(1)+T(2)+T(3)+V_0; 

        if T(4) >= cmin
            % Undamaged tumour cells  
            DT(1) =  q_2*T(4)*T(1)*(K - E) - l1*R*T(4)*T(1)...
                     - l2*R*T(4)*T(1) + mu*T(2);      

            % Sublethally damaged cells     
            DT(2) =  q_2sl*T(4)*T(2)*(K - E) + l1*R*T(4)*T(1)...
                     - mu*T(2) - (zeta + l2sl*R*T(4))*T(2); 

            % Lethally damaged cells            
            DT(3) =  l2*R*T(4)*T(1) + (zeta + l2sl*R*T(4))*T(2) - eta*T(3);   

           % Oxygen concentration
            DT(4) =  g*(1 - T(4))*V_0 - (q_1*T(1) + q_1sl*T(2))*T(4)...
                     - (q_3*T(1) + q_3sl*T(2))*(K -  E)*T(4);     

        elseif T(4) < cmin

            % Undamaged tumour cells
            DT(1) =  q_2*T(4)*T(1)*(K - E) - l1*R*T(4)*T(1)...
                     - l2*R*T(4)*T(1) + mu*T(2) - d1*(cmin - T(4))*T(1); 

            % Sublethally damaged cells
            DT(2) =  q_2sl*T(4)*T(2)*(K - E) + l1*R*T(4)*T(1) - mu*T(2)...
                     - (zeta + l2sl*R*T(4))*T(2) - d1sl*(cmin - T(4))*T(2); 

            % Lethally damaged cells
            DT(3) =  l2*R*T(4)*T(1) + (zeta + l2sl*R*T(4))*T(2) - eta*T(3);     

            % Oxygen concentration
            DT(4) =  g*(1 - T(4))*V_0 - (q_1*T(1) + q_1sl*T(2))*T(4)...
                     - (q_3*T(1) + q_3sl*T(2))*(K -  E)*T(4);                

        end

       radiotherapy = DT;
   
    end

