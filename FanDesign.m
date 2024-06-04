classdef FanDesign
    %FANDESIGN Correct Axial Fan preliminary design following
    %[kromerSoundEmissionLowpressure2017]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% Github Repo Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% A custom MATLAB® class definition containing "Class members" for
    %%%%% FanDesign: properties, methods, functions, events
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%  Description  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Blade element theory is used in the design methodology for low pressure
    % axial fans. The following code includes a free vortex and non-free vortex 
    % design, the latter using multiple controlled vortex approaches, influencing 
    % the loading distrubution. CV3 found suitiable for low noise designs%
    properties
        %% Roman Letters

        cm2             % Absolute Frame Meridonal Velocity

        d_bsi           % Blade section diamaters
        d_duct          % Diameter of Fan Duct/Shroud (INPUT)
        d_fan           % Diameter of Fan tip (INPUT)
        d_hub           % Diameter of Fan Hub (INPUT)
        k               % Intregration Constant k for Controlled Vortex Methods
        loading         % Blade Loading Profile Switch Case (INPUT)
        n               % Rotational Speed (INPUT)
        P_amb           % Pressure_Atmosphere (INPUT)
        P_b             % Blade Power
        profile_name    % x/y coordinates of airfoil from XFOIL (INPUT)
        r_bsi           % Blade section radii
        r_duct          % Radius of Fan Duct/Shroud
        r_fan           % Radius of Fan tip
        r_hub           % Radius of Fan Hub
        S               % Surface area between Hub and Tip of Fan,
        s_b             % Blade section spacing
        T               % Temperature (INPUT)
        tip_gap         % Tip Gap between blade and Duct
        Vdot            % Volume Flow Rate (INPUT)
        Vdot_d          % Actual Volumetric flow rate based on efficiency
        Vdot_dpre       % Preliminary Actual Volumetric flow rate based on efficiency

        Vel_U1          % Blade speed at inlet 
        Vel_U2          % Blade speed at exit
        Vel_c_m1        % Meridional velocity at inlet
        Vel_c_u1        % Circumferential velocity at inlet
        Vel_c_u2        % Circumferential velocity at exit
        Vel_c_m2        % Meridional velocity at exit
        Vel_C1          % Absolute velocity at inlet
        Vel_C2          % Absolute velocity at exit
        Vel_W1          % Relative velocity at inlet
        Vel_W2          % Relative velocity at exit
        Vel_w_m1        % Relative Meridional velocity at inlet
        Vel_w_m2        % Relative Meridional velocity at exit
        Vel_w_u1        % Relative Circumferential velocity at inlet
        Vel_w_u2        % Relative Circumferential velocity at exit
        Vel_w_mean      % Relative Mean velocity

        Vortex_a        % Control Coefficients for Controlled Vortex Loading distribution
        Vortex_b        % Control Coefficients for Controlled Vortex Loading distribution
        Vortex_d        % Control Coefficients for Controlled Vortex Loading distribution
        Vortex_e        % Control Coefficients for Controlled Vortex Loading distribution
        Vortex_m        % Control Coefficients for Controlled Vortex Loading distribution

        Y_bd            % Specific Blade Energy
        Y_bsd           % Specific Blade Energy at Blade Segments
        Y_tt
        zb              % Fan Blade number (INPUT)
        z_bf            % Blade section numbers (number of spanwise slices) (INPUT)

        %% Greek Letters
        alpha1_bf       % Angle of attack for each section (INPUT)
        alpha2_bf
        beta1_bf
        beta2_bf
        beta_bf_mean
        delta_bf        % Sweep angle for each section in Degrees (INPUT)
        dP_ts           % Total-to-Static Pressure Difference (INPUT)
        dP_tspre        % Preliminary Total-to-Static Pressure
        dP_tt           % Total-to-Total Pressure Difference
        dP_ttd          % Total-to-Total Pressure Difference Actual
        dP_ttd_cor      % Total-to-Total Pressure Difference Actual and corrected for blade sweep
        eta_h           % Hydraulic Efficiency (INPUT)
        eta_vol         % Volumetric Efficiency (INPUT)
        mu              % Viscosity
        nu_bf           % Dihedral angle for each blade section in Degrees  (INPUT)
        rho             % Density
        omega           % Rotational velocity (rad/s)

        nu              % Diameter Ratio
        phi             % Flow Coefficient

        psi_tt          % Total-to-Total Pressure Coefficient
        psi_ts          % Total-to-Static Pressure Coefficient
        eta_ts
        sigma_s         % Specific Speed
        delta_d         % Specific Diameter

    end

    methods
        function obj = FanDesign(inputArg1,inputArg2,inputArg3,inputArg4)
            %FANDESIGN Construction method for Fan Design class
            %   For robustness and condensed inputs, common parameters
            %   groups into inputArg1


            obj.Vdot = inputArg1(1);

            obj.d_fan = inputArg1(2);
            obj.d_hub = inputArg1(3);
            obj.d_duct = inputArg1(4);

            obj.r_fan = obj.d_fan/2;
            obj.r_hub = obj.d_hub/2;
            obj.r_duct = obj.d_duct/2;

            obj.n = inputArg1(5);

            obj.tip_gap = (obj.d_duct - obj.d_fan)/2;


            obj.zb = inputArg1(6);
            obj.eta_vol = inputArg1(7);
            obj.eta_h  = inputArg1(8);
            obj.z_bf = inputArg1(9);
            obj.loading = inputArg1(10); % Switch Case Statement 0=FV, 1=CV1, 2=CV2, 3=CV3

            %CHECK
            if obj.loading ~= 0 || obj.loading ~= 1 || obj.loading ~= 2 || obj.loading ~= 3
                % disp('ERROR - No valid blade loading distribution chosen (0=FV, 1=CV1, 2=CV2, 3=CV3)');
            end
            
            %%%%% inputArg1(2) %%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%
            obj.alpha1_bf = inputArg2(1,:);
            obj.delta_bf = inputArg2(2,:);
            obj.nu_bf    = inputArg2(3,:);

            if length(obj.alpha1_bf) ~= obj.z_bf
                disp('ERROR - Number of alpha angles not equal to section number');
            elseif length(obj.delta_bf) ~= obj.z_bf
                disp('ERROR - Number of sweep angles not equal to section number');
            elseif length(obj.nu_bf) ~= obj.z_bf
                disp('ERROR - Number of dihedral angles not equal to section number');
            end

            obj.profile_name = inputArg3;

            obj.Vortex_a = inputArg4(1);
            obj.Vortex_b = inputArg4(2);
            obj.Vortex_d = inputArg4(3);
            obj.Vortex_e = inputArg4(4);
            obj.Vortex_m = inputArg4(5);

            % Some basic parameters to calculate from Inputs 
            obj.T = inputArg1(11);
            obj.P_amb = inputArg1(12);
            obj.dP_ts = inputArg1(13);
            obj.S = (pi/4)*(obj.d_fan^2 - obj.d_hub^2);
            obj.omega = (2*pi*obj.n)/60;

            % Execute rest of Fan design methods here
            obj = thermodynamicProperties(obj);
            obj = segmentation(obj);


            % Depending on Blade Loading Profile use IF statement to decide
            % flow of programme

            % FREE VORTEX LOADING
            if obj.loading == 0
                obj = loadingtype0(obj);
            end

            if obj.loading == 1
                obj = loadingtype1(obj);
            end



        end

        function obj = thermodynamicProperties(obj)

            % Ideal Gas Density
            obj.rho = obj.P_amb / (287.052874*obj.T);

            % Sutherlands Law - Good for low pressure single component gases
            obj.mu = (1.716e-5)*((obj.T/273)^(3/2))*((273 + 111)/(obj.T + 111));


        end

        function obj = segmentation(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            for i = 1:obj.z_bf
                obj.d_bsi(i) = sqrt(obj.d_hub^2 + ((i-1)*(obj.d_fan^2 - obj.d_hub^2))/(obj.z_bf -1));
            end
            obj.r_bsi = obj.d_bsi/2;

            obj.s_b = (pi*obj.d_bsi)/obj.zb;
        end

        function obj = loadingtype0(obj)
            % LOADINGTYPE0 Execute fan design using Free vortex approach
            disp('Hello you are using FV Loading Criteria')

            obj.Vdot_d = obj.Vdot/obj.eta_vol;

            obj.cm2 = (4*obj.Vdot_d)/(pi*(obj.d_duct^2 - obj.d_hub^2));

            % Solution of total-to-total pressure difference for free
            % vortex design based on non-linear equation. Define a function
            % with extra parameters (defining the preliminary fan) and
            % solve using fzero (total-to-static dP) as guess.
            myfunc = @(x,A,B,C,D,E,F,G,H,I) -A + B*x - (pi/(C*D))*((x^2)/(E^2))*log(F/G) - 0.5*C*((4*H)/(pi*I))^2;

            A = obj.dP_ts;
            B = obj.eta_h;
            C = obj.rho;
            D = obj.S;
            E = obj.omega;
            F = obj.r_fan;
            G = obj.r_hub;
            H = obj.Vdot_d;
            I = obj.d_duct^2 - obj.d_hub^2;

            fun = @(x) myfunc(x,A,B,C,D,E,F,G,H,I);

            obj.dP_ttd = fzero(fun,obj.dP_ts); %,optimset('Display','iter'));

            obj.dP_tt = obj.eta_h*obj.dP_ttd;

            %% Blade Skew Correction

            if all(obj.delta_bf == obj.delta_bf(1))
                % Constant Sweep Angle
                obj.dP_ttd_cor = (1/((cosd(obj.delta_bf(1)))^0.62))*obj.dP_ttd;
            else
                disp('Run equation 4.19 using numerical integration')
            end

            % Specific Blade Energy
            obj.Y_bd = obj.dP_ttd/obj.rho;
            % Specific Blade Energy of each section
            obj.Y_bsd = obj.Y_bd*ones(1,obj.z_bf);


            %% Velocity Triangle Calculations

            for i = 1:obj.z_bf

                % Entry Triangle
                obj.Vel_U1(i) = 2*pi*obj.n*obj.r_bsi(i)/60;

                obj.Vel_c_m1(i) = obj.Vdot/(pi*(obj.r_fan^2 - obj.r_hub^2));
                obj.Vel_w_m1(i) = obj.Vdot/(pi*(obj.r_fan^2 - obj.r_hub^2));
                obj.Vel_c_u1(i) = tand(obj.alpha1_bf(i))*obj.Vel_c_m1(i);       %obj.Vel_c_m1(i)/tand(obj.alpha1_bf(i));
                obj.Vel_w_u1(i) = obj.Vel_U1(i) - obj.Vel_c_u1(i);
                obj.beta1_bf(i) = atand(obj.Vel_w_m1(i)/obj.Vel_w_u1(i));       %atand(obj.Vel_c_m1(i)/(obj.Vel_U1(i)-obj.Vel_c_u1(i)));

                obj.Vel_W1(i) = sqrt(obj.Vel_w_m1(i)^2 + obj.Vel_w_u1(i)^2);    %sqrt( obj.Vel_U1(i)^2 + obj.Vel_C1(i)^2 - 2*obj.Vel_U1(i)*obj.Vel_C1(i)*cosd(obj.alpha1_bf(i)));
                obj.Vel_C1(i) = sqrt(obj.Vel_c_u1(i)^2 + obj.Vel_c_m1(i)^2);

                % Exit Triangle
                obj.Vel_U2(i) = obj.Vel_U1(i);

                obj.Vel_c_m2(i) = obj.Vdot/(pi*(obj.r_fan^2 - obj.r_hub^2));
                obj.Vel_w_m2(i) = obj.Vdot/(pi*(obj.r_fan^2 - obj.r_hub^2));
                obj.Vel_c_u2(i) = obj.Vortex_b./obj.r_bsi(i);
                obj.alpha2_bf(i) = atand(obj.Vel_c_m2(i)./obj.Vel_c_u2(i));

                obj.Vel_C2(i) = sqrt( obj.Vel_c_m2(i)^2 + obj.Vel_c_u2(i)^2);
                obj.Vel_W2(i) = sqrt( obj.Vel_U2(i)^2 + obj.Vel_C2(i)^2 - 2*obj.Vel_U2(i)*obj.Vel_C2(i)*cosd(obj.alpha2_bf(i)));

                obj.Vel_w_u2(i) = sqrt( obj.Vel_W2(i)^2 - obj.Vel_w_m2(i)^2);

                obj.beta2_bf(i) = atand(obj.Vel_w_m2(i)/obj.Vel_w_u2(i));

                obj.Vel_w_mean(i) = 0.5*sqrt((obj.Vel_U1(i) + sqrt( obj.Vel_W2(i)^2 - obj.Vel_c_m1(i)^2))^2 + 4*obj.Vel_c_m1(i)^2);
                obj.beta_bf_mean(i) = atand((2*obj.Vel_c_m1(i))/(obj.Vel_U1(i) + sqrt(obj.Vel_W2(i)^2 - obj.Vel_c_m1(i)^2)));

            end

            %% Dimensionless Numbers
            obj.Y_tt = obj.Y_bd*obj.eta_h; % Total-to-Total Specific Energy

            obj.nu = obj.d_hub/obj.d_fan;   % Diameter Ratio

            obj.phi = (4*obj.Vdot)/(pi*(obj.d_fan^2 - obj.d_hub^2)*(obj.omega*obj.r_fan)); % flow Coefficient

            obj.psi_tt = 2*obj.Y_tt/(obj.omega*obj.r_fan)^2; % Total-to-Total Pressure Coefficient

            obj.psi_ts = (2*obj.dP_ts)/(obj.rho*((obj.omega*obj.r_fan)^2)); % Total-to-Static Pressure Coefficient

            obj.eta_ts = [];

            obj.sigma_s = (obj.psi_tt^(-0.75))*(obj.phi^(+0.5));

            obj.delta_d = (obj.psi_tt^(+0.25))*(obj.phi^(-0.5));

        end

        function obj = loadingtype1(obj)
            disp('Hello you are using CV1 Loading Criteria')
            % Controlled Vortex Distrubtion Type 1 rCu2 = ar + b
            obj.dP_tspre = 0;

            % while obj.dP_tspre > obj.dP_ts*1.01 || obj.dP_tspre <
            % obj.dP_ts*0.99 

            A = obj.Vortex_a;
            B = obj.Vortex_b;

            int_k = 0; %Set initial guess of integration constant k

            obj.Vdot_dpre = 0;

            obj.Vdot_d = obj.Vdot/obj.eta_vol;

            while obj.Vdot_dpre > obj.Vdot_d*1.001 || obj.Vdot_dpre < obj.Vdot_d*0.999
                
                % Controlled Vortex Type 1 (Eq 3.41)
                A = obj.omega*obj.eta_h*obj.Vortex_a.*obj.r_bsi;
                B = (obj.Vortex_a.^2).*(log(obj.r_bsi));
                C = (obj.Vortex_a.*obj.Vortex_b)./obj.r_bsi;

                D = sqrt( 2.*(A - B + C) + int_k);
                obj.Vel_c_m2 = D;
                obj.Vel_c_m2 = sqrt( 2.*(A - B + C) + int_k);

                obj.Vdot_dpre = trapz(obj.r_bsi,obj.Vel_c_m2.*obj.r_bsi)*2*pi;

                if obj.Vdot_dpre < obj.Vdot_d
                    int_k = int_k+0.01;
                elseif obj.Vdot_dpre > obj.Vdot_d
                    int_k = int_k-0.01;
                end
            end
            % Log Integration K into object 
            obj.k = int_k;

            obj.Vel_c_u2 = (obj.Vortex_a.*obj.r_bsi + obj.Vortex_b)./obj.r_bsi;            % Absolute circumferential Velocity at exit from vortex profile
            obj.Vel_U1 = 2.*pi.*obj.n.*obj.r_bsi./60;   % Blade Speed at exit = inlet 

            % Mass Average Velocity of Cm2 and Cu1
            func1 = (obj.Vel_c_u2.^2).*(obj.Vel_c_m2).*obj.rho.*2*pi.*obj.r_bsi;
            A = trapz(obj.r_bsi,func1)
            func2 = (obj.Vel_c_m2.^2).*(obj.Vel_c_m2).*obj.rho.*2*pi.*obj.r_bsi;
            B = trapz(obj.r_bsi,func2)

            % obj.dP_ttd = obj.dP_ts + A + B;

            obj.dP_ttd = obj.dP_ts + 0.5*obj.rho.*(A + B);
            

            % Blade Power
            obj.P_b = 2*pi*obj.rho*obj.omega.*trapz(obj.r_bsi,obj.Vel_c_m2.*obj.Vel_c_u2.*(obj.r_bsi.^2));
            % Specific Blade Energy
            obj.Y_bd = obj.P_b./(obj.rho.*obj.Vdot_dpre);
            % Specific Blade Energy for each section
            obj.Y_bsd = obj.omega.*obj.r_bsi.*obj.Vel_c_u2;

            % Blade Skew Correction
            % kk = (cosd(obj.delta_bf).^(0.62));
            func1 = obj.dP_ttd*(cosd(obj.delta_bf).^(0.62))*obj.rho.*obj.Vel_c_m2.*2.*pi.*obj.r_bsi;
            func2 = obj.dP_ttd*obj.rho.*obj.Vel_c_m2.*2.*pi.*obj.r_bsi;

            func1r = trapz(obj.r_bsi,func1);
            func2r = trapz(obj.r_bsi,func2);

            obj.dP_ttd_cor = (func2r/func1r)*obj.dP_ttd;

            %% Velocity Triangles

            % Entry Triangle
            obj.Vel_U1 = 2*pi*obj.n*obj.r_bsi/60;

            obj.Vel_c_m1 = obj.Vdot./(pi*(obj.r_fan.^2 - obj.r_hub.^2));%%%%%% cm1=cm2 (see p29 & Appendix B)
            obj.Vel_w_m1 = obj.Vdot./(pi*(obj.r_fan.^2 - obj.r_hub.^2));
            obj.Vel_c_u1 = tand(obj.alpha1_bf).*obj.Vel_c_m1;       %obj.Vel_c_m1(i)/tand(obj.alpha1_bf(i));
            obj.Vel_w_u1 = obj.Vel_U1 - obj.Vel_c_u1;
            obj.beta1_bf = atand(obj.Vel_w_m1./obj.Vel_w_u1);       %atand(obj.Vel_c_m1(i)/(obj.Vel_U1(i)-obj.Vel_c_u1(i)));

            obj.Vel_W1 = sqrt(obj.Vel_w_m1.^2 + obj.Vel_w_u1.^2);    %sqrt( obj.Vel_U1(i)^2 + obj.Vel_C1(i)^2 - 2*obj.Vel_U1(i)*obj.Vel_C1(i)*cosd(obj.alpha1_bf(i)));
            obj.Vel_C1 = sqrt(obj.Vel_c_u1.^2 + obj.Vel_c_m1.^2);

            % Exit Triangle
            obj.Vel_U2 = obj.Vel_U1;

            % obj.Vel_c_m2(i) = obj.Vdot/(pi*(obj.r_fan^2 - obj.r_hub^2));
            obj.Vel_w_m2 = obj.Vel_c_m2; %obj.Vdot/(pi*(obj.r_fan^2 - obj.r_hub^2));
            % obj.Vel_c_u2 = obj.Vortex_b./obj.r_bsi;
            obj.alpha2_bf = atand(obj.Vel_c_m2./obj.Vel_c_u2);

            obj.Vel_C2 = sqrt( obj.Vel_c_m2.^2 + obj.Vel_c_u2.^2);
            obj.Vel_W2 = sqrt( obj.Vel_U2.^2 + obj.Vel_C2.^2 - 2*obj.Vel_U2.*obj.Vel_C2.*cosd(obj.alpha2_bf));

            obj.Vel_w_u2 = sqrt( obj.Vel_W2.^2 - obj.Vel_w_m2.^2);

            obj.beta2_bf = atand(obj.Vel_w_m2./obj.Vel_w_u2);

            obj.Vel_w_mean = 0.5*sqrt((obj.Vel_U1 + sqrt( obj.Vel_W2.^2 - obj.Vel_c_m1.^2)).^2 + 4*obj.Vel_c_m1.^2);
            obj.beta_bf_mean = atand((2*obj.Vel_c_m1)./(obj.Vel_U1 + sqrt(obj.Vel_W2.^2 - obj.Vel_c_m1.^2)));

            %% Total-to-Static Pressure Check

            % Need to recalculate total-to-static pressure using calculated
            % values and update the loop.
            % Use the dP_tt = dP_tt + dynamic equation again, due to the
            % correction
            %

            obj.dP_tt = obj.dP_ttd_cor*obj.eta_h;
            obj.dP_tspre = obj.dP_ttd_cor*obj.eta_h - 0.5*obj.rho.*(A+B);
            Out = obj.dP_tspre;
            b = obj.Vortex_b;
        
            %
            %     if obj.dP_tspre > obj.dP_ts
            %
            %         obj.Vortex_b = obj.Vortex_b - 0.01;
            %
            %     elseif obj.dP_tspre < obj.dP_ts
            %
            %         obj.Vortex_b = obj.Vortex_b + 0.01;
            %
        
        %
        % end

        % while (CHECK parameter B)
        % loading distribution
        % while (CHECK A)
        % total pressure difference
        % blade Skew
        % Velocity Diagrams
        % Total-to-static pressure
        % end
        % end
        % dimensionless number s


    end


end
end

