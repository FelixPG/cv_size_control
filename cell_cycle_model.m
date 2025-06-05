%% Numerical simulation code of the 2017 Chandler-Brown cell cycle model
clear all 

% Path where the sim workspace will be saved
path = pwd;

%% Simulation hyperparameters

% Define simulation length parameters; time/iteration = 1 minute
% Main parameters for the simulation, change below
start_phase_length = 2000;      % Number of iterations in simulation, default 2000, more than 1600 is fine.
start_phase_pop_cap = 10000;    % Maximum number of cells during steady state phase, default 10000, more than 6000 is fine.

control_pop_size = 1;           % 0 = like flask growth, 1 = like chemostat
generate_random_pop_start = 0;  % 0 = random seed using parameters below; 1 = load seed pop; 2 = load and seed with contraint
max_seed_size = 100;
min_seed_size = 40;
n_cells = 10;          %starting number of cells

pop_mu_i = 3e4;          %femtoliters conversion
pop_sigma_i = .17 * pop_mu_i;   %femtoliters conversion

%% Cell Parameters

% Taking a sensitivity vector to perturb parameters in the desired
% direction.
% Set epsilon to 0 to take baseline WT values or choose ppv=ppv_null
% ppv is a parameter perturbation vector in the 21 dimension of the
% parameter space.
% ppv below are identified from the parameter sensitivity analysis. 

epsilon = +0.5;

ppv_null = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

ppv_cv_full = [0.2,0.091,-0.12,0.16,0.38,-0.18,-0.18,-0.5,-0.35,0.5,0,0,0,0,0.12,0.06,0,-0.12,0.08,0,0.19]; % Gradient of CV Full distribution (mother and daughter)
ppv_cv_start = [0,0.27,-0.26,-0.15,-0.26,-0.11,0.26,-0.25,0.12,-0.16,0.65,-0.25,0.13,0,-0.063,0.14,0,0,0.2,-0.074,0.12]; % Gradient of the CV of the Daughter Start distribution

%ppv_random = 2*rand(21,1)-1; % Vector with components sampled from Unif in [-1,1[ 
%ppv_random = ppv_random/norm(ppv_random); % Normalize the parameter direction

% Choose which vector to use
ppv = ppv_null;

% Set the cell parameters for the simulation taking into account the ppv 
% perturbation if needed.

%%% DAUGHTERS %%%

% START size control daughters
% Slope (M), offset
G1_lambda = [1.85056126550248e-06 -0.0337256466753454]; % coeffvalues(min_fit1_act);
bud_prob_beta_1 = G1_lambda(1)*(1+epsilon*ppv(1));      % per minute; Slope
bud_prob_M0_1 = G1_lambda(2)*(1+epsilon*ppv(2));      % fl; Threshold

%Post-START G1 timer daughters
% Slope (M), offset
poststart_G1_timer = [-2.596953686853002e-04*(1+epsilon*ppv(3)) 34.316569344194875*(1+epsilon*ppv(4))];

% Division control, daughters
% Slope (M_Bud), Slope (M_START), offset
G2_lambda_2d = [(7.628399268975367e-06)*(1 + epsilon*ppv(5)) -2.983508662043306e-06*(1+epsilon*ppv(6)) -0.011854207243603*(1+epsilon*ppv(7))];%exp_2d_coeff;

% Constant growth rate
% for Daughters and Mothers in S/G2/M
all_fit_act = [0.004999129499166 0];
constant_growth_rate = (log(all_fit_act(1)*3 + 1)/3)*(1+epsilon*ppv(8)); %min_a_act; %per min

% Bud mass correction, daughters
% Growth in S/G2/M is in the bud, but some goes into the main cell body
bud_mass_correction = 1.369421416415018e+03*(1+epsilon*ppv(16));

% S/G2/M timer daughter cells
SG2M_timer = [-2.681657833310124e-04*(1+epsilon*ppv(17)) 1.094494905033565e+02*(1+epsilon*ppv(18))];%SG2M_timer_fit;

%%% MOTHERS %%%

% Post-START G1 timer mothers
% Slope (M), offset
mother_G1_timer = [(-5.429621092521212e-04)*(1+epsilon*ppv(9)) 47.73555293824496*(1+epsilon*ppv(10))];%mother_G1_fit;

% Division control, mothers
% Slope (M_bud), Slope (M_START), offset
mother_lambda = [6.610122886272357e-06*(1+epsilon*ppv(11)) -2.538505485426652e-06*(1+epsilon*ppv(12)) -0.006267310306677*(1+epsilon*ppv(13))];%exp_2d_coeff_mother;

% Mother growth in Post-START G1 is bilinear, NOT exponential!
mother_g1_gr_fit = [0.001953883911940*(1+epsilon*ppv(14)) 98.154228802411230*(1+epsilon*ppv(15))];

% S/G2/M Timer, mothers
mother_G2_timer = [-2.979402851076718e-04*(1+epsilon*ppv(19)) 1.049573409490653e+02*(1+epsilon*ppv(20))];%mother_G2_dur;

% Bud mass correction, mothers
% Growth in S/G2/M is in the bud, but some goes into the main cell body
mother_bud_mass_correction = 2.417144252429131e+03*(1+epsilon*ppv(21));%m_b_diff;

%%
% Build starting population

if generate_random_pop_start == 0

    mu = log((pop_mu_i^2)/sqrt(pop_sigma_i+pop_mu_i^2));    %values for log-normal
    sigma = sqrt(log(pop_sigma_i/(pop_mu_i^2)+1));  %values for log-normal
    cell_v = lognrnd(mu,sigma,n_cells,1);    %sample from log-normal volume dist
    cell_cycle = zeros(n_cells,1);  %begin with G1 cells; G1 = 0, G2 = 1
    cell_bud = zeros(n_cells,1);  %begin with no buds
    G2_counter = zeros(n_cells,1);  %cells count down G2
    mother_daughter = zeros(n_cells,1); %begin with all daughters; daughter = 0, mother = 1
    mother_G1_counter = zeros(n_cells,1); %all cells start in G1 as daughters
    post_start_G1 = zeros(n_cells,1);
    post_start_counter = zeros(n_cells,1);
    start_size = zeros(n_cells,1);
    G2_length_daughter = []; %time in G2 (if buds)
    mother_counter = [];

elseif generate_random_pop_start == 1

    load('Seed_Data_Size_Partitioning_Sim.mat');
    load('Starting_Parameters_Temp.mat');
    delete('Starting_Parameters_Temp.mat');
    seed_time = length(cell_v(1,:));
    choose_cells = randsample(length(cell_v(:,seed_time)),n_cells); %select cells from previous population
    cell_v = cell_v(choose_cells,seed_time);
    cell_cycle = cell_cycle(choose_cells,seed_time);  %begin with G1 cells; G1 = 0, G2 = 1
    cell_bud = cell_bud(choose_cells,seed_time);  %begin with no buds
    G2_counter = G2_counter(choose_cells,seed_time);  %cells count down G2
    mother_daughter = mother_daughter(choose_cells,seed_time); %begin with all daughters; daughter = 0, mother = 1
    mother_G1_counter = mother_G1_counter(choose_cells,seed_time); %all cells start in G1 as daughters
    post_start_G1 = post_start_G1(choose_cells,seed_time);
    post_start_counter = post_start_counter(choose_cells,seed_time);
    mother_bud_mass_defect = mother_bud_mass_defect(choose_cells);
    start_size = start_size(choose_cells);
    mother_counter = mother_counter(choose_cells);
    G2_length_mother = G2_length_mother(choose_cells);

end

%% Simulation

%Grow cells
%Cells first choose between bud, divide, grow; cells whose divide counters
%are greater than 0 will grow bud, else, divide; cells who are in G1 and
%don't get a bud draw will grow mother; cells with bud draws will initiate
%bud and grow bud

%Diagnostic Parameters
pop_vol_mean = zeros(start_phase_length,1);
pop_vol_sigma = zeros(start_phase_length,1);
mother_fraction = zeros(start_phase_length,1);
pop_cell_cycle = zeros(start_phase_length,1);

bud_prob = [];  %chance of budding (volume dependent)
bud = [];   %does this cell bud (if in G1)?

for i = 1:start_phase_length

    ['Current time: ' num2str(i)]

    % Update diagnostics
    pop_vol_mean(i) = mean(cell_v(:,i));
    pop_vol_sigma(i) = std(cell_v(:,i));
    cell_num(i) = length(cell_v(:,i));
    mother_fraction(i) = mean(mother_daughter(:,i));
    pop_cell_cycle(i) = mean(cell_cycle(:,i));

    ['# of cells: ',num2str(cell_num(i))]

    %Begin simulation - define global event stats
    bud_prob = bud_prob_beta_1.*cell_v(:,i) + bud_prob_M0_1; %piecewise linear fits from data
    bud_prob(bud_prob < 0) = 0;
    
    %bud_prob = 0.02*ones(length(cell_v(:,i))); % Constant timer rate in
    %pre-Start G1. No more size control.

    bud_prob = 1 - exp(-bud_prob*(1));
    
    bud = zeros(length(bud_prob),1);
    bud = transpose(bud);

    for k = 1:length(cell_v(:,i))   %for all cells
        growth_size = cell_v(k,i) + cell_bud(k,i);

        if mother_daughter(k,i) == 0 %if cells are newborn
            % Growth
            dv = (exp(constant_growth_rate) - 1) * growth_size;
            dv(dv < 0) = 0;
            if cell_cycle(k,i) == 0 %if cells are in G1
                if post_start_G1(k,i) == 0 %if cells are pre-Start
                    
                    % Deterministic START transition. pre-Start G1 becomes
                    % a deterministic sizer
                    %bud_thresh = 3.3e4;
                    %if cell_v(k,i) >= bud_thresh
                    %    bud_prob(k) = 1;
                    %else
                    %    bud_prob(k) = 0;
                    %end
                     
                    bud(k) = randsample([0 1], 1, true, [1-bud_prob(k), ...
                        bud_prob(k)]);

                    if bud(k) == 0
                        cell_v(k,i+1) = cell_v(k,i) + dv;
                        cell_cycle(k,i+1) = 0; %stays in G1
                        cell_bud(k,i+1) = 0;
                        mother_daughter(k,i+1) = mother_daughter(k,i);
                        G2_counter(k,i+1) = 0;
                        mother_G1_counter(k,i+1) = 0;
                        post_start_G1(k,i+1) = 0;
                        post_start_counter(k,i+1) = 0;
                        mother_counter(k) = 0;
                    elseif bud(k) == 1
                        cell_v(k,i+1) = cell_v(k,i) + dv;
                        cell_cycle(k,i+1) = 0; %stays in G1
                        cell_bud(k,i+1) = 0;
                        mother_daughter(k,i+1) = mother_daughter(k,i);
                        G2_counter(k,i+1) = 1;
                        mother_G1_counter(k,i+1) = 0;
                        post_start_G1(k,i+1) = 1;
                        post_start_counter(k,i+1) = ...
                            round(polyval(poststart_G1_timer,cell_v(k,i)));
                        if post_start_counter(k,i+1) < 0
                            post_start_counter(k,i+1) = 0;
                        end
                        start_size(k) = cell_v(k,i);
                        G2_length(k) = polyval(SG2M_timer,start_size(k));
                        if G2_length(k) < 0
                            G2_length(k) = 0;
                        end
                        mother_bud_mass_defect(k) = bud_mass_correction / G2_length(k);
                        mother_counter(k) = 0;
                    end
                elseif post_start_G1(k,i) == 1 %if cells are post-Start
                    if post_start_counter(k,i) > 0
                        cell_v(k,i+1) = cell_v(k,i) + dv;
                        cell_cycle(k,i+1) = 0; %stays in G1
                        cell_bud(k,i+1) = 0;
                        mother_daughter(k,i+1) = mother_daughter(k,i);
                        G2_counter(k,i+1) = G2_counter(k,i);
                        mother_G1_counter(k,i+1) = 0;
                        post_start_G1(k,i+1) = 1;
                        post_start_counter(k,i+1) = post_start_counter(k,i) - 1;
                        mother_counter(k) = 0;
                    elseif post_start_counter(k,i) == 0
                        cell_v(k,i+1) = cell_v(k,i); %bud grows, mother does not
                        cell_cycle(k,i+1) = 1; %change to G2
                        G2_counter(k,i+1) = G2_counter(k,i); %randomly select length of G2
                        cell_bud(k,i+1) = dv;
                        mother_daughter(k,i+1) = mother_daughter(k,i);
                        mother_G1_counter(k,i+1) = 0;
                        post_start_G1(k,i+1) = 0;
                        post_start_counter(k,i+1) = 0;
                        mother_counter(k) = 0;
                    end
                end
            elseif  cell_cycle(k,i) == 1 %if cells are in G2
                div_prob = G2_lambda_2d(1)*cell_bud(k,i) + ...
                    G2_lambda_2d(2)*start_size(k) + G2_lambda_2d(3);

                if div_prob < 0
                    div_prob = 0;
                end

                div_prob = 1 - exp(-div_prob*(1));

                G2_counter(k,i) = ...
                    randsample([0 1], 1, true, [div_prob, 1 - div_prob]);

                if G2_counter(k,i) > 0 %growth in G2
                    cell_v(k,i+1) = cell_v(k,i) + mother_bud_mass_defect(k); %bud grows, mother does not
                    cell_bud(k,i+1) = cell_bud(k,i) + ...
                        dv - mother_bud_mass_defect(k); %grow bud
                    G2_counter(k,i+1) = G2_counter(k,i) - 1; %shorten remaining G2
                    cell_cycle(k,i+1) = 1; %Stay in G2
                    mother_daughter(k,i+1) = mother_daughter(k,i);
                    mother_G1_counter(k,i+1) = 0;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    mother_counter(k) = mother_counter(k) + 1;
                elseif G2_counter(k,i) == 0
                    new_cell_num = length(cell_v(:,i))+1;
                    
                    % Division from 1st gen cells!
                    %div_frac = 2; %If you want to do symmetrical divisions
                    %div_frac_noise = sqrt(0.05)*randn; %with noise
                    cell_v(new_cell_num,i+1) = cell_bud(k,i);%cell_bud(k,i); %Make new daughter (cell_v(k,i)+cell_bud(k,i))/(div_frac+div_frac_noise)
                    cell_v(k,i+1) = cell_v(k,i);%cell_v(k,i); %No growth during cytokinesis (cell_v(k,i)+cell_bud(k,i))/(div_frac-div_frac_noise)

                    cell_bud(new_cell_num,i+1) = 0; %No bud on daughter
                    cell_cycle(new_cell_num,i+1) = 0; %New cell in G1
                    cell_bud(k,i+1) = 0; %Seperate bud from mother
                    cell_cycle(k,i+1) = 0; %Return mother to G1
                    mother_daughter(k,i+1) = 1; %Not first time mother
                    mother_daughter(new_cell_num,i+1) = 0; %New never budded daughter
                    mother_G1_counter(k,i+1) = ...
                        round(polyval(mother_G1_timer,cell_v(k,i)));
                    if mother_G1_counter(k,i+1) < 0
                        mother_G1_counter(k,i+1) = 0;
                    end
                    mother_G1_counter(new_cell_num,i+1) = 0;
                    G2_counter(k,i+1) = 0;
                    G2_counter(new_cell_num,i+1) = 0;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    post_start_G1(new_cell_num,i+1) = 0;
                    post_start_counter(new_cell_num,i+1) = 0;
                    mother_bud_mass_defect(new_cell_num) = 0;
                    mother_bud_mass_defect(k) = 0;
                    start_size(k) = cell_v(k,i+1);
                    start_size(new_cell_num) = 0;
                    mother_counter(k) = 0;
                    mother_counter(new_cell_num) = 0;
                    G2_length_mother(k) = ...
                        polyval(mother_G2_timer,cell_v(k,i));
                    G2_length_mother(new_cell_num) = 0;
                    M_D_parent_at_birth(new_cell_num) = mother_daughter(k,i);

                end

            end
        elseif mother_daughter(k,i) == 1 % if cells are 2nd gen
            growth_size = cell_v(k,i) + cell_bud(k,i);
            if cell_cycle(k,i) == 0
                dv = ((mother_g1_gr_fit(1)*growth_size + ...
                    mother_g1_gr_fit(2))*exp(mother_g1_gr_fit(1)) - ...
                    mother_g1_gr_fit(2))/mother_g1_gr_fit(1) - growth_size;
                dv(dv < 0) = 0;
                if mother_G1_counter(k,i) > 0 %if cell does not bud this time
                    cell_v(k,i+1) = cell_v(k,i) + dv;
                    cell_cycle(k,i+1) = 0; %stays in G1
                    cell_bud(k,i+1) = 0;
                    G2_counter(k,i+1) = 0;
                    mother_daughter(k,i+1) = mother_daughter(k,i);
                    mother_G1_counter(k,i+1) = mother_G1_counter(k,i) - 1;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    mother_counter(k) = mother_counter(k) + 1;

                elseif mother_G1_counter(k,i) == 0 %if cell buds
                    cell_v(k,i+1) = cell_v(k,i); %bud grows, mother does not
                    cell_cycle(k,i+1) = 1; %change to G2
                    G2_counter(k,i+1) = 1; %randomly select length of G2
                    cell_bud(k,i+1) = dv;
                    mother_daughter(k,i+1) = mother_daughter(k,i);
                    mother_G1_counter(k,i+1) = 0;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    mother_bud_mass_defect(k) = mother_bud_mass_correction / G2_length_mother(k);
                    mother_counter(k) = mother_counter(k) + 1;

                end

            elseif cell_cycle(k,i) == 1

                dv = (exp(constant_growth_rate) - 1) * growth_size;%dv = (exp(constant_growth_rate) - 1) * growth_size;
                dv(dv < 0) = 0;

                div_prob = mother_lambda(1)*cell_bud(k,i) + ...
                    mother_lambda(2)*start_size(k) + mother_lambda(3);              
                div_prob = 1 - exp(-div_prob*(1));
                div_prob(div_prob < 0) = 0;

                G2_counter(k,i) = ...
                    randsample([0 1], 1, true, [div_prob, 1 - div_prob]);

                if G2_counter(k,i) > 0 %growth in G2
                    cell_v(k,i+1) = cell_v(k,i) + mother_bud_mass_defect(k); %bud grows, mother does not
                    cell_bud(k,i+1) = cell_bud(k,i) + ...
                        dv - mother_bud_mass_defect(k); %grow bud
                    G2_counter(k,i+1) = G2_counter(k,i); %shorten remaining G2
                    cell_cycle(k,i+1) = 1; %Stay in G2
                    mother_daughter(k,i+1) = mother_daughter(k,i);
                    mother_G1_counter(k,i+1) = 0;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    mother_counter(k) = mother_counter(k) + 1;

                elseif G2_counter(k,i) == 0
                    new_cell_num = length(cell_v(:,i))+1;

                    %Division from 2nd gen cells or later
                    %div_frac = 2;
                    %div_frac_noise = sqrt(0.05)*randn;
                    cell_v(new_cell_num,i+1) = cell_bud(k,i);%(cell_v(k,i)+cell_bud(k,i))/(div_frac+div_frac_noise);%cell_bud(k,i); %Make new daughter
                    cell_v(k,i+1) = cell_v(k,i);%(cell_v(k,i)+cell_bud(k,i))/(div_frac-div_frac_noise);%cell_v(k,i); %No growth during cytokinesis

                    cell_bud(new_cell_num,i+1) = 0; %No bud on daughter
                    cell_cycle(new_cell_num,i+1) = 0; %New cell in G1
                    cell_bud(k,i+1) = 0; %Seperate bud from mother
                    cell_cycle(k,i+1) = 0; %Return mother to G1
                    mother_daughter(k,i+1) = 1; %Not first time mother
                    mother_daughter(new_cell_num,i+1) = 0; %New never budded daughter
                    mother_G1_counter(k,i+1) = ...
                        round(polyval(mother_G1_timer,cell_v(k,i)));
                    if mother_G1_counter(k,i+1) < 0
                        mother_G1_counter(k,i+1) = 0;
                    end
                    mother_G1_counter(new_cell_num,i+1) = 0;
                    G2_counter(k,i+1) = 0;
                    G2_counter(new_cell_num,i+1) = 0;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    post_start_G1(new_cell_num,i+1) = 0;
                    post_start_counter(new_cell_num,i+1) = 0;
                    mother_bud_mass_defect(new_cell_num) = 0;
                    mother_bud_mass_defect(k) = 0;
                    start_size(k) = cell_v(k,i+1);
                    start_size(new_cell_num) = 0;
                    mother_counter(k) = 0;
                    mother_counter(new_cell_num) = 0;
                    G2_length_mother(k) = ...
                        polyval(mother_G2_timer,cell_v(k,i));
                    if G2_length_mother(k) < 0;
                            G2_length_mother(k) = 0;
                        end
                    G2_length_mother(new_cell_num) = 0;

                    M_D_parent_at_birth(new_cell_num) = mother_daughter(k,i);

                end
            end
        end
    end

    if control_pop_size == 1
        select_cells = [];
        if length(cell_v(:,i+1)) > start_phase_pop_cap %caps population by randomly removing excess cells
            select_cells = randsample(length(cell_v(:,i+1)),start_phase_pop_cap);
            cell_v = cell_v(select_cells,:);
            cell_cycle = cell_cycle(select_cells,:);
            cell_bud = cell_bud(select_cells,:);
            G2_counter = G2_counter(select_cells,:);
            mother_G1_counter = mother_G1_counter(select_cells,:);
            mother_daughter = mother_daughter(select_cells,:);
            post_start_G1 = post_start_G1(select_cells,:);
            post_start_counter = post_start_counter(select_cells,:);
            mother_bud_mass_defect(select_cells);
            M_D_parent_at_birth = M_D_parent_at_birth(select_cells);
            start_size = start_size(select_cells);
            mother_counter = mother_counter(select_cells);
            G2_length_mother = G2_length_mother(select_cells);

        end
    end

    if length(find(cell_v(:,i+1) == 0)) > 0 % Negative growth is an error condition
        ['Negative Growth']
        return
    end
end

%%
% Data Processing

cell_cycle_change = [];
START_change = [];
birth_time = [];
first_daughter_bud = [];
first_daughter_cytokinesis = [];
G1_growth = [];
size_at_birth = [];
size_at_bud =  [];
cell_cycle_growth = [];
SG2M_growth = [];
first_daughter_start = [];
bud_size_at_cytokinesis = [];
volume_at_START = [];
n_data_cells = 0;
cell_cycle_length = [];
g1_length = [];
g2_length = [];
daughter_final_size = [];

for cell = 1:length(cell_v(:,1))
    cell_cycle_change(cell,:) = diff(cell_cycle(cell,:));
    START_change(cell,:) = diff(post_start_G1(cell,:));
    if isempty(find(cell_cycle_change(cell, :) == 1, 1)) == 0
        if isempty(find(cell_cycle_change(cell, :) == -1, 1)) == 0
            if isempty(find(START_change(cell,:) == 1, 1)) == 0
                n_data_cells = n_data_cells + 1;
                birth_time(n_data_cells) = find(cell_v(cell,:) > 0, 1);
                
                first_daughter_bud(n_data_cells) = find(cell_cycle_change(cell, :) == 1, 1);
                first_daughter_cytokinesis(n_data_cells) = find(cell_cycle_change(cell, :) == -1, 1);
                first_daughter_start(n_data_cells) = find(START_change(cell,:) == 1, 1);
                
                cell_cycle_length(n_data_cells) = first_daughter_cytokinesis(n_data_cells) - birth_time(n_data_cells);
                g1_length(n_data_cells) = first_daughter_bud(n_data_cells) - birth_time(n_data_cells);
                g2_length(n_data_cells) = first_daughter_cytokinesis(n_data_cells) - first_daughter_bud(n_data_cells);
                
                daughter_final_size(n_data_cells) = cell_v(cell,end)+cell_bud(cell,end);
                G1_growth(n_data_cells) = cell_v(cell,first_daughter_bud(n_data_cells)) - ...
                    cell_v(cell,birth_time(n_data_cells));

                size_at_birth(n_data_cells) = cell_v(cell,birth_time(n_data_cells));
                size_at_bud(n_data_cells) = cell_v(cell,first_daughter_bud(n_data_cells));

                cell_cycle_growth(n_data_cells) = cell_v(cell,first_daughter_cytokinesis(n_data_cells)) - ...
                    cell_v(cell,birth_time(n_data_cells)) + cell_bud(cell,first_daughter_cytokinesis(n_data_cells));

                SG2M_growth(n_data_cells) = cell_v(cell,first_daughter_cytokinesis(n_data_cells)) - ...
                    cell_v(cell,first_daughter_bud(n_data_cells)) + cell_bud(cell,first_daughter_cytokinesis(n_data_cells));

                bud_size_at_cytokinesis(n_data_cells) = cell_bud(cell, ...
                    first_daughter_cytokinesis(n_data_cells) - 1);
                volume_at_START(n_data_cells) = cell_v(cell,first_daughter_start(n_data_cells));
                
                M_D_at_birth_processed(n_data_cells) = M_D_parent_at_birth(cell);
               
            end
        end
    end
end

first_gen_mother_birth = [];
first_gen_mother_bud = [];
first_gen_mother_cytokinesis = [];
mother_birth_size = [];
mother_bud_size = [];
mother_cytokinesis_size = [];
mother_adder = [];
mother_final_size = [];

n_mother_cells = 0;
for cell = 1:length(cell_v(:,1))
    if length(find(cell_cycle_change(cell, :) == 1, 2)) > 1
        if length(find(cell_cycle_change(cell, :) == -1, 2)) > 1  
            n_mother_cells = n_mother_cells + 1;

            m_birth = find(cell_cycle_change(cell, :) == -1, 1);
            first_gen_mother_birth(n_mother_cells) = ...
                m_birth(1) + 1; % Time of birth for first gen mother cell
            m_bud = find(cell_cycle_change(cell, :) == 1, 2);
            first_gen_mother_bud(n_mother_cells) = ...
                m_bud(2); % Time of budding for first gen mother cell
            m_cytokinesis = find(cell_cycle_change(cell, :) == -1, 2);
            first_gen_mother_cytokinesis(n_mother_cells) = ...
                m_cytokinesis(2); % Time of cytokinesis for first gen mother cell

            mother_g1_length = first_gen_mother_bud - first_gen_mother_birth; % Post-START G1 length for first gen mother cell
            mother_g2_length = first_gen_mother_cytokinesis - first_gen_mother_bud; % G2 length first gen mother cell
            mother_cycle_length = first_gen_mother_cytokinesis - first_gen_mother_birth; % Cycle length for first gen mother cell

            mother_final_size(n_mother_cells) = cell_v(cell,end)+cell_bud(cell,end); % First gen mother cell size at end of sim
            mother_birth_size(n_mother_cells) = ...
                cell_v(cell,first_gen_mother_birth(n_mother_cells)) + ...
                cell_bud(cell,first_gen_mother_birth(n_mother_cells)); % Size at birth (START) for first gen mother cell
            mother_cytokinesis_size(n_mother_cells) = ...
                cell_v(cell,first_gen_mother_cytokinesis(n_mother_cells)) + ...
                cell_bud(cell,first_gen_mother_cytokinesis(n_mother_cells)); % Size at cytokinesis for first gen mother cell
            mother_adder(n_mother_cells) = ...
                mother_cytokinesis_size(n_mother_cells) - ...
                mother_birth_size(n_mother_cells); % Added volume over cycle for first gen mother cell
            mother_bud_size(n_mother_cells) = ...
                cell_v(cell,first_gen_mother_bud(n_mother_cells)) + ...
                cell_bud(cell,first_gen_mother_bud(n_mother_cells)); % Size at budding for first gen mother cell
        end
    end
end

%%%%% Filter cells after steady state only %%%%%%%
steady_state_time = 1000; %Time to start investigating cells

cells_after_steady_state = find(birth_time > steady_state_time); % Daughters born after steady_state_time
mothers_after_steady_state = find(first_gen_mother_birth > steady_state_time); % 1st gen. mothers born after steady_state_time

% Daughters
birth_time = birth_time(cells_after_steady_state);
first_daughter_bud = first_daughter_bud(cells_after_steady_state);
first_daughter_cytokinesis = first_daughter_cytokinesis(cells_after_steady_state);
cell_cycle_length = cell_cycle_length(cells_after_steady_state);
g1_length = g1_length(cells_after_steady_state);
g2_length = g2_length(cells_after_steady_state);
G1_growth = G1_growth(cells_after_steady_state);
size_at_birth = size_at_birth(cells_after_steady_state);
size_at_bud =  size_at_bud(cells_after_steady_state);
cell_cycle_growth = cell_cycle_growth(cells_after_steady_state);
SG2M_growth = SG2M_growth(cells_after_steady_state);
bud_size_at_cytokinesis = bud_size_at_cytokinesis(cells_after_steady_state);
volume_at_START = volume_at_START(cells_after_steady_state);
first_daughter_start = first_daughter_start(cells_after_steady_state);
M_D_at_birth_processed = M_D_at_birth_processed(cells_after_steady_state);
daughter_final_size = daughter_final_size(cells_after_steady_state);

% 1st gen. mothers
first_gen_mother_birth = first_gen_mother_birth(mothers_after_steady_state);
first_gen_mother_bud = first_gen_mother_bud(mothers_after_steady_state);
first_gen_mother_cytokinesis = first_gen_mother_cytokinesis(mothers_after_steady_state);
mother_g1_length = mother_g1_length(mothers_after_steady_state);
mother_g2_length = mother_g2_length(mothers_after_steady_state);
mother_cycle_length = mother_cycle_length(mothers_after_steady_state);
mother_final_size = mother_final_size(mothers_after_steady_state);
mother_birth_size = mother_birth_size(mothers_after_steady_state);
mother_bud_size = mother_bud_size(mothers_after_steady_state);
mother_cytokinesis_size = mother_cytokinesis_size(mothers_after_steady_state);
mother_adder = mother_adder(mothers_after_steady_state);

%%%%% Stats. for daughters from born from daughters or 1st gen mothers,
%%%%% Added by Felix 
size_at_birth_from_mother = [];
size_at_birth_from_daughter = [];
size_at_start_from_mother = [];
size_at_start_from_daughter = [];
size_at_bud_from_mother = [];
size_at_bud_from_daughter = [];
size_at_cytokinesis_from_mother = [];
size_at_cytokinesis_from_daughter = [];
size_at_birth_from_mother = [];
size_at_birth_from_daughter = [];

for i = 1:length(M_D_at_birth_processed)
    if M_D_at_birth_processed(i) == 1
        size_at_start_from_mother(i) = volume_at_START(i);
        size_at_birth_from_mother(i) = size_at_birth(i);
        size_at_bud_from_mother(i) = size_at_bud(i);
        size_at_cytokinesis_from_mother(i) = size_at_birth(i)+cell_cycle_growth(i);
    else
        size_at_start_from_daughter(i) = volume_at_START(i);
        size_at_birth_from_daughter(i) = size_at_birth(i);
        size_at_bud_from_daughter(i) = size_at_bud(i);
        size_at_cytokinesis_from_daughter(i) = size_at_birth(i)+cell_cycle_growth(i);
    end 
end

% Save the workspace. Analyze the run using the analyze_run.mlx notebook
save(strcat(path,'/sim_workspace.mat'))