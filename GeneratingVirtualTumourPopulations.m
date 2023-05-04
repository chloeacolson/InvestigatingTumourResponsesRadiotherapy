%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Chloe Colson
% Date: 21 Sep 2022
% Function: Code to construct the nutrient limited (NL), space limited (SL)
%           and bistable (BS) virtual tumour populations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Define the model parameters

    % Anoxic O2 threshold
    
        c_min = 1e-2;
    
    % Rate of oxygen release from vasculature
    
        g = 5;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% 1. NL regime

    % Specify the seed for the random number generator
    
        rng(1)
    
    % Define the vascular volume
    
        V0 = 0.0005;
    
    % Define the threshold value of q1 (q1_threshold) above which 
    % (q1,q3,V0) corresponds to the NL regime for all q1 > q1_threshold 
    % and all q3 > 0
    
        q1_threshold = g*(1/cmin -1)*(V0/(1-V0));
        
    
    % Randomly select q1 and q3 from the uniform distributions
    % U(q1_threshold,10) and U(0.01,10), respectively
    
        q1 = q1_threshold + (10-q1_threshold)*rand(250,1); 
        q3 = 0.01 + (10-0.01)*rand(250,1);
    

    % Calculating steady state values for T and c for each (q1,q3,V0)
    % parameter set
    
        for i = 1:size(q1,1)
            c_SS(i,1) = (cmin*(q1(i) -3*q3(i)+(g/cmin+q3(i))*V0 + ...
                        sqrt((q1(i)-3*q3(i)+(g/cmin+q3(i))*V0)^2 +...
                        4*q3(i)*(2*(q1(i)-q3(i))+(g-q1(i)+q3(i))*V0))))...
                        /2*(2*(q1(i)-q3(i))+(g-q1(i)+q3(i))*V0);
            T_SS(i,1) = 2-V0-cmin/c_SS(i);
        end
    
    % Create table listing the parameter values for the N=250 NL virtual 
    % tumours
    
        ID = (1:250)';
        V0 = repelem(V0,250,1);
        
        NL_cohort = array2table(horzcat(ID,q1,q3,V0,c_SS,T_SS),...
                    'VariableNames', {'ID','q1','q3','V0','c0','T0'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% 2. SL regime

    
    % Specify the seed for the random number generator
    
        rng(1)
    
    % Define the vascular volume
    
        V0 = 0.005;
    
    % Define the threshold value of q1 (q1_threshold) such that (q1,q3,V0) 
    % corresponds to the SL or BS regimes for all q1 < q1_threshold 

        q1_threshold = g*(1/cmin -1)*(V0/(1-V0));   
    
    % Randomly select q1 values from the uniform distribution
    % U(0.01, q1_threshold) 
    
        q1 = 0.01 + (q1_threshold-0.01)*rand(250,1);


    % For each q1 value selected:
    % a. Define the threshold values q3 (q3_threshold) such that (q1,q3,V0)
    %    corresponds to the SL regime for all q3 <= q3_threshold
    % b. Randomly select a q3 value from the uniform distribution 
    %    U(0.01,q3_threhold)
    

        for i = 1:size(q1,1)
            q3_threshold = (-q1(i)+(3*g/cmin-2+q1(i))*V0-(g/cmin)*V0^2)...
                            /(-1+V0)^2 + (2*sqrt(g))/(-1+V0)^2* ...
                            (sqrt((-(2/cmin-1)*q1(i)*V0+(g*(2/(cmin^2)- ...
                            3/cmin+1)+q1(i)*(3/cmin-1))*V0^2 ...
                            -(g*(1/(cmin^2)-1/cmin)+q1(i)/cmin)*V0^3)));


            q3(i,1) = 0.01 + (q3_threshold-0.01)*rand(1,1);  
            
            % Ensure q3 < 10 (as q3_threshold can be greater than 10)
            
                while ((q3(i,1)) > 10)
                    q3(i,1) = 0.01 + (q3_threshold-0.01)*rand(1,1);
                end
        end

    % Calculating steady state values for T and c for each (q1,q3,V0)
    % parameter set
    
        for i = 1:size(q1,1)
            c_SS(i,1) = V0/(V0 + (q1(i)/g)*(1-V0));
            T_SS(i,1) = 1-V0;
        end 
    
    % Create table listing the parameter values for the N=250 NL virtual 
    % tumours 
    
        ID = (1:250)';
        V0 = repelem(V0,250,1);
        
        SL_cohort = array2table(horzcat(ID,q1,q3,V0,c_SS,T_SS),...
        'VariableNames', {'ID','q1','q3','V0','c0','T0'});
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% BS regime

    % Specify the seed for the random number generator
    
        rng(1)
    
    % Define the vascular volume
    
        V0 = 0.00275;
    
    % Define the threshold value of q1 (q1_threshold) such that (q1,q3,V0) 
    % corresponds to the SL or BS regimes for all q1 < q1_threshold 

        q1_threshold = g*(1/cmin -1)*(V0/(1-V0));   
    
    % Randomly select q1 values from the uniform distribution
    % U(0.01, q1_threshold) 
    
        q1 = 0.01 + (q1_threshold-0.01)*rand(250,1);


    % For each q1 value selected:
    % a. Define the threshold values q3 (q3_threshold) such that (q1,q3,V0)
    %    corresponds to the BS regime for all q3 > q3_threshold
    % b. Randomly select a q3 value from the uniform distribution 
    %    U(q3_threhold,10)
    
        for i = 1:size(q1,1)
            q3_threshold = (-q1(i)+(3*g/cmin-2+q1(i))*V0-(g/cmin)*V0^2)...
                            /(-1+V0)^2 + (2*sqrt(g))/(-1+V0)^2* ...
                            (sqrt((-(2/cmin-1)*q1(i)*V0+(g*(2/(cmin^2)- ...
                            3/cmin+1)+q1(i)*(3/cmin-1))*V0^2 ...
                            -(g*(1/(cmin^2)-1/cmin)+q1(i)/cmin)*V0^3)));


            q3(i,1) = q3_threshold + (10-q3_threshold)*rand(1,1);     
        end
       
    % Calculating steady state values for T and c for each (q1,q3,V0)
    % parameter set
    
        for i = 1:size(q1,1)
            c_SS(i,1) = (cmin*(q1(i) -3*q3(i)+(g/cmin+q3(i))*V0 + ...
                        sqrt((q1(i)-3*q3(i)+(g/cmin+q3(i))*V0)^2 +...
                        4*q3(i)*(2*(q1(i)-q3(i))+(g-q1(i)+q3(i))*V0))))...
                        /2*(2*(q1(i)-q3(i))+(g-q1(i)+q3(i))*V0);
            T_SS(i,1) = 2-V0-cmin/c_SS(i);
        end
    
    % Create table listing the parameter values for the N=250 NL virtual 
    % tumours
    
        ID = (1:250)';
        V0 = repelem(V0,250,1);
        
        BS_cohort = array2table(horzcat(ID,q1,q3,V0,c_SS,T_SS),...
        'VariableNames', {'ID','q1','q3','V0','c0','T0'});
