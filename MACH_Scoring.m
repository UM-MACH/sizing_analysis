% ===========================================================
%                                         .:,:`           
%                                      `:,:               
%     ...       ...     `;;;          :,,     ;;;     ,;; 
%     ,,,:     ,,,,     ####        `,,:      ###     ;## 
%     ,,,,,   .,,,,     #####      `,,:       ###     ;## 
%     ,,,,,. `,,,,,    +##+##      ,,,        ###     ;## 
%     ,,,,,, ,,,,,,    ### ###    ,,,.        ###     ;## 
%     ,,, ,,,,, ,,,   .##` ###   ,,,,  #+     ########### 
%     ,,,  ,,,  ,,,   ###  `##;  *######+.,,,,########### 
%     ,,,  `,.  ,,,   ###   ###   ::,       .,,:.     :## 
%     ,,,   `   ,,,  ##########.  :,,         '`,,,,, ;## 
%     ,,,       ,,,  ###########   ,,,        ##@.,,,,, # 
%     ,,,       ,,, .##;     ###    ,,        ###  ,,,,,` 
%     ,,,       ,,: ###       ###    :,       ###  :,,,,, 
%                                   `:                 
% ===========================================================
%
% SUMMARY: Scoring analysis for the 2017-2018 DBF competition        

%% Initialization
% close all
clear all
% clc

%% =========================== Unit Conversions ============================= %%

in2m = 0.0254; %inches to meters
ft2m = in2m * 12; %feet to meters

%% =========================== Varibles ===================================== %%

AR = 4:8;
bombs = 4:8;
[AR, bombs] = meshgrid(AR, bombs);

g = 9.81; %m/s^2
num_wings = 1; %wings

weight_fuselage_intial = 10;  % fuse empty  Newtons

weight_fuselage = weight_fuselage_intial; % Newtons

%wing linear density - weight of the wing per inch
dens_lin_wing = 0.2875 * 0.0283495 * g / in2m;  % density N/m

thrust_to_weight = .7;
Takeoff_velocity = 12; %m/s
CD_0 = 0.06;

e = 0.80; %oswald efficiency
lap_length = 5000; % ft
airfoil_Cl_max = 1.4;

%historical data
air_density = 1.17; %air (kg/m^3)

%symbolic variables
v = sym('v','real');

wing_ref_area = 0.8; % inital guess for wing area (m^2)
weight_propulsion = 12; % inital guess for Propulsion System weight in Newtons



%% ================ Data for Propulsion System Weight ======================= %%

thrust_data = [0.5, 1, 2, 3, 4, 5, 15]*4.44822; % Thrust needed in Newtons
weight_data =[9, 12, 25.9, 37.5, 42, 55, 113]*0.28; % Weight of prop. sys. in N

prop_weight_coef = polyfit(thrust_data, weight_data, 1 );


%% ========================= Size Aircraft ================================= %%
% Iterative method to solve for wing_ref_area, AR, weight_propulsion, MTOW, weight_empty
% Note: this method is not very sophisticated and prone to error 

for i = 0:1000

    % calculate the wingspan from the aspect ratio and wing reference area
    span_wing = sqrt(AR .* wing_ref_area); %wingspan (m)

    % wing weight in Newtons 
    weight_wings = num_wings*(span_wing*dens_lin_wing);

    
    for k = 0:100

        weight_empty = weight_fuselage + weight_wings + weight_propulsion; %empty weight (N)
        MTOW = weight_empty + bombs * (0.0850486 * g); %Max Takeoff weight (N) (bombs weigh 3 oz 0.085 kg)
        
        thrust = thrust_to_weight*MTOW; %thrust (N)

        err = sum(sum(abs(polyval(prop_weight_coef, thrust) - weight_propulsion))); %sum of absolute error
        fprintf('err %f\n',err);

        if err < 1e-8
            fprintf('thrust converged after %d\n',k);
            break

        end
        
       weight_propulsion = weight_propulsion + ...
        ( polyval(prop_weight_coef, thrust) - weight_propulsion) ;

    end
    
   
    Cl_stall = airfoil_Cl_max * AR ./ (AR + 2); % finite wing correction 
    Cl_takeoff = Cl_stall/(1.1^2); % equation from 481
    wing_area_req = 2*MTOW./(Cl_takeoff .* air_density * (Takeoff_velocity^2)); %required wing area


    err_wing_ref_area = sum(sum(abs(wing_ref_area - wing_area_req)));
    fprintf('err  %d  %f\n ', i, err_wing_ref_area);

    if err_wing_ref_area < 1e-8
        fprintf('wing_ref_area converged after %d\n',i);
        break

    end

    wing_ref_area = wing_ref_area + 0.1*(wing_area_req - wing_ref_area);  %m^2


    % % AR = (span_wing.*0.0254).^2./wing_ref_area*2;
    % % EW = (Weight_Fuselage + 1 * 0.0283495 * passengers) + (dens_lin_wing* span_wing).*(span_wing./AR / 12 * 0.5 + 0.5) + P;
    % % MTOW = (EW+passengers*2*0.0283495)*g


    % if sum(sum( abs( (span_wing.*0.0254).^2./wing_ref_area*2 - AR) )) < 1e-2
    %     disp('break')
    %     break
    % end
    % i
end
% return
AR( AR<= 0) = 0.1
AR(isnan(AR)==1) = 0.1
MTOW( MTOW<= 0) = 1e6
MTOW(isnan(MTOW)==1) = 0.1
wing_ref_area( wing_ref_area<= 0) = inf
wing_ref_area( isnan(wing_ref_area)== 1) = 0.01


% %% ======================= Calculations ================================ %%


% calculate wing_ref_area based on throwing speed 

% span_wingpan = (AR .* wing_ref_area).^0.5;
% Chord_root = wing_ref_area./span_wingpan/(1+lambda)*2;

% Chord = sqrt(wing_ref_area./AR); %meters
% span_wingpan = wing_ref_area./Chord; %wingspan in meters
%     Aspect_R = span_wingpan./Chord


T_0 = thrust_to_weight*MTOW;
T_1 = -0.060;
T_2 = -0.015;


for j = [1:size(bombs,1)]
    j
    for i = [1:size(span_wing,2)]
        v_cruise(j,i) = double(max(vpasolve( v.^2*T_2 + v*T_1 + T_0(j,i) ==...
        0.5*air_density*v.^2*wing_ref_area(j,i).*...
        (CD_0 + 1./(pi*e*AR(j,i))*(2*MTOW(j,i)./(air_density   *v^2*wing_ref_area(j,i))).^2)...
        , v)));

        if imag(v_cruise(j,i)) || real(v_cruise(j,i)) < 0
            v_cruise(j,i) = 0;
%             disp('yes')
        end
        
    end
end

laps = (v_cruise)./(lap_length/3.28) * 60 * 10;
t_3laps = 3*(lap_length/3.28)./v_cruise*0.8;
RAC = weight_empty*0.224809.*span_wing;



%only account for well define aircraft
for j = [1:size(bombs,1)]
    for i = [1:size(span_wing,2)]
        if(j*3-45 > i)
            t_3laps(j,i) = 100000000;
            laps(j,i) = 0;
            
        end
    end
end        

v_cruise(v_cruise>25) = 25;


M1 = zeros(size(bombs));
M1(t_3laps/60 < 5) = 1;

M2 = 2*bombs./(t_3laps)./max(max(bombs./(t_3laps))); %mission 2 score
M3 =4*(bombs/2.*(bombs*2/2).*laps)./max(max((bombs/2.*(bombs*2/2).*laps))) + 2;

score = (M1 + M2 + M3)./RAC;
% score(M2 == 0) = 0;
% maximum = max(max(score));
% score = score./maximum;
% %% ========================== Plotting ================================= %%

figure
shading interp;
[ANALY2 , ANALY2] = contourf( bombs, span_wing, score);
set(ANALY2,'edgecolor','none');

title('Scoring Analysis','FontSize',23);


ylabel('b_s','FontSize',20);
xlabel('passengers','FontSize',20);
hold on 

max_score = max(max(score));
scatter(bombs(score== max_score), span_wing(score== max_score),250, 'ko','filled')

c = colorbar;
c.Label.String = 'Normalized Score';

fprintf('Max Score %.2f\n',max_score);
fprintf('RAC %.2f\n', RAC(score== max_score));
fprintf('passengers %.2f\n', bombs(score== max_score));
fprintf('span_wing [in]%.2f\n', span_wing(score== max_score));
fprintf('MTOW [N]%.2f\n', MTOW(score== max_score));
fprintf('wing_ref_area [m^2]%.2f\n', wing_ref_area(score== max_score));
fprintf('span_wingpan [in] %.2f\n', span_wing(score== max_score));
fprintf('Chord [in] %.2f\n', span_wing(score== max_score)/AR(score== max_score));
fprintf('AR %.2f\n', AR(score== max_score));
fprintf('V %.2f\n', v_cruise(score== max_score));
fprintf('Cl %.2f\n', Cl_takeoff(score== max_score));




figure
contourf( bombs, span_wing, RAC);
title('RAC','FontSize',23);
ylabel('b_s','FontSize',20);
xlabel('passengers','FontSize',20);
hold on 
c = colorbar;


% figure
% contourf( passengers, span_wing, (M1 + M2 + M3));
% title('(M1 + M2 + M3)','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, MTOW);
% title('MTOW','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, weight_empty);
% title('weight_empty','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, AR);
% title('AR','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, wing_ref_area);
% title('s','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, v_cruise);
% title('v','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;
