clear; clc; clear all;
tic;
%% global variable definition
global model lambda GAP_vector record   % Global variable should also be defined in the sub-function.
lambda = 100; GAP_vector = [];    % global variable initialization
%% code information pre-statement
% This code uses single penalty factor % single model order.
% The optimization solver in this code is gurobi.
% This code generates random matrix using Gaussian random function: normrnd(mu, sigma, m, n).
% This code is suited to direct load control for ATDM: from zones to building cluster.
% This code make use of the function 'squeeze' to compress the dimension of matrix
%% import data
raw_data = load('welldata.mat');
raw_data=raw_data.welldata.data_2014_11;    % change directory for simplicity
% model = load('model_building.mat');
% model = model.model;    % change directory for simplicity
%% parameter setting 
num_sample = 1440;    % sample number
num_training = 0.75 * num_sample;    % training set number 
num_verifying = 0.25 * num_sample;    %verification set number
% building_id = [1 8 10 11 13 16 20]';     % building ID
building_id = [1 8 10 11 13 16 20]';
% building_id = [1]'; 
num_building = length(building_id);    % building number
model_order = 2;    % model order 
%% data pre_processing
% overall data acquisition (all buildings combined together)
model.overall_data.h_load = 0;
model.overall_data.tau_amb = [];
model.overall_data.radiation = [];
model.overall_data.tau_in = [];
for k=1:num_building
    temp_focused_building_data = raw_data(building_id(k,1)).building;
    model.individual_data(k).h_load = temp_focused_building_data(:,1);
    model.individual_data(k).tau_amb = temp_focused_building_data(:,2);
    model.individual_data(k).radiation = temp_focused_building_data(:,3);
    model.individual_data(k).tau_in = temp_focused_building_data(:,4:end);
    model.overall_data.h_load = model.overall_data.h_load + model.individual_data(k).h_load;
    model.overall_data.tau_in =  [model.overall_data.tau_in, model.individual_data(k).tau_in];
    clear temp_focused_building_data;    
end
model.overall_data.tau_amb = model.individual_data(1).tau_amb;
model.overall_data.radiation = model.individual_data(1).radiation;
% training data acquisition
model.overall_training_data.h_load = model.overall_data.h_load(1:num_training,:);
model.overall_training_data.tau_amb = model.overall_data.tau_amb(1:num_training,:);
model.overall_training_data.radiation = model.overall_data.radiation(1:num_training,:);
model.overall_training_data.tau_in = model.overall_data.tau_in(1:num_training,:);
% verification data acquisition
model.overall_verifying_data.h_load = model.overall_data.h_load(num_training+1:end,:);
model.overall_verifying_data.tau_amb = model.overall_data.tau_amb(num_training+1:end,:);
model.overall_verifying_data.radiation = model.overall_data.radiation(num_training+1:end,:);
model.overall_verifying_data.tau_in = model.overall_data.tau_in(num_training+1:end,:);
num_zone = size(model.overall_data.tau_in,2);    % number of zones in the building cluster
%% Perform the Privacy-Preserved Block Coordinated Descent Method Function
 PP_BCDM(model_order, num_training, num_zone);
%% Indicator calculation
% RMSE calculation
%   model.record_aggregate_I_full(i_penalty, i_order).indicator.RMSE = ...
%             sqrt(mean((model.nonlinear.verify.Tau_in_a_original(set_sample_training,:) - ...
%             model.nonlinear.verify.Tau_in_a_fore(set_sample_training,:)).^2));
%         model.record_aggregate_I_full(i_penalty, i_order).indicator.MAPE = ...
%             mean(abs((model.nonlinear.verify.Tau_in_a_original(set_sample_training,:) - ...
%             model.nonlinear.verify.Tau_in_a_fore(set_sample_training,:)) ./ ...
%             model.nonlinear.verify.Tau_in_a_fore(set_sample_training,:)))*100;
%         model.record_aggregate_I_full(i_penalty, i_order).indicator.R_square = ...
%             1 - ...
%             sum((model.nonlinear.verify.Tau_in_a_original(set_sample_training,:) - ...
%             model.nonlinear.verify.Tau_in_a_fore(set_sample_training,:)).^2) / ...
%             sum((model.nonlinear.verify.Tau_in_a_original(set_sample_training,:) - ...
%             mean(model.nonlinear.verify.Tau_in_a_original(set_sample_training,:))).^2);
model.training.aggregated_tau_in_ori = model.overall_training_data.tau_in*model.var.xi; % 1080*1
model.verifying.aggregated_tau_in_ori = model.overall_verifying_data.tau_in*model.var.xi; % 360*1
model.verifying.aggregated_tau_in_pre = zeros(num_verifying,1); %360*1
model.verifying.aggregated_tau_in_pre(1:model_order,:)=model.verifying.aggregated_tau_in_ori(1:model_order,:);
for i =model_order+1 :num_verifying % 逐点预测
    model.verifying.alpha_term = zeros(1,1);
    model.verifying.beta_term =zeros(1,1);
    model.verifying.gamma_term =zeros(1,1);
    model.verifying.theta_term = zeros(1,1);
    for j = 2:model_order+1
        model.verifying.alpha_term =  model.verifying.alpha_term + model.var.alpha(1,j) * model.verifying.aggregated_tau_in_pre(i-(j-1),1);
    end
    for j = 1:model_order+1
        model.verifying.beta_term = model.verifying.beta_term +model.var.beta(1,j) *model.overall_verifying_data.h_load(i-(j-1),1);
        model.verifying.gamma_term = model.verifying.gamma_term+model.var.gamma(1,j)*model.overall_verifying_data.tau_amb(i-(j-1),1);
        model.verifying.theta_term = model.verifying.theta_term+model.var.theta(1,j)*model.overall_verifying_data.radiation(i-(j-1),1);
    end
    model.verifying.aggregated_tau_in_pre(i,1)=model.verifying.alpha_term+model.verifying.beta_term+model.verifying.gamma_term+model.verifying.theta_term+model.var.occ(1+mod(i+num_training-1,48),1);
end

plot(model.verifying.aggregated_tau_in_ori);hold on; plot(model.verifying.aggregated_tau_in_pre);

model.verifying.error = model.verifying.aggregated_tau_in_ori-model.verifying.aggregated_tau_in_pre;
model.verifying.rmse = sqrt(mean(model.verifying.error.^2));
model.verifying.mpae = mean(abs(model.verifying.error./model.verifying.aggregated_tau_in_ori))*100;
model.verifying.r2= 1-sum(model.verifying.error.^2)/sum((model.verifying.aggregated_tau_in_ori-mean(model.verifying.aggregated_tau_in_ori)).^2);

 save('pp_data.mat','model');
 toc;
 disp(['operating time:', num2str(toc)]);
 
%% Privacy Preserved Block Coordinated Descent Method (PP-BCDM) for Aggregated Thermal Dynamic Model (ATDM)
function PP_BCDM (model_order, num_training,num_zone)
    %% global variable definition
    global model lambda GAP_vector record
    %% predefition & preprocess
    model.M_total = [];
    model.tau_in_hat_total =[];
    % solver related definition
    yalmip('clear');
    mysolver = 'gurobi'; 
    % iteration related definition
    Is_convergence = false;    % judge whether the iteration is converged; initial state: not converged
    num_iteration = 0;    % iteration number
    GAP = + inf;    % GAP = min{f_1 - f_2, (f_1-f_2)/f_2}. GAP is used to judge whether the iteration is terminated.
    GAP_vector = [];    % restore the GAP values in each iteration
    % ATDM parameter initialization
    model.var.xi_initial = 1 / num_zone * ones(num_zone, 1);    % aggregated coefficients' initialization: take its average value
    model.var.xi = model.var.xi_initial;    % In sub-problem one, xi is fixed. Because we first solve sub-problem one, so xi is fixed first.
    %% perform iterations
    while ~Is_convergence
        num_iteration = num_iteration + 1;
    M = normrnd(0.1, 0.1, num_zone, num_zone);    % generate transformation matrix M
%     M = eye(num_zone);
    model.M_total = [model.M_total;M];
    %% sub-problem one 
    % variable 
    model.var.alpha = sdpvar(1,model_order+1, 'full');
    model.var.beta = sdpvar(1,model_order+1, 'full');
    model.var.gamma = sdpvar(1,model_order+1, 'full');
    model.var.theta = sdpvar(1,model_order+1, 'full');
    model.var.occ = sdpvar(num_training-model_order,1, 'full');    % variable depicting occupant behavior
    model.var.residual = sdpvar(num_training-model_order,1, 'full');    % residual variable
    % constraints
    model.cons = [];
    model.cons = [model.cons model.var.alpha(1,1) == 1];    % alpha_0 = 1
    model.cons = [model.cons, model.var.occ(49 : num_training - model_order, 1) == model.var.occ(1 : num_training - model_order - 48, 1)];    % constraints for periodic occupant behaviour
    % data randomizatoin
    % data randomization for tau_in
    r_matrix = normrnd(0, 1, num_zone, num_zone, num_training);   % r_matrix for encryption
%     r_matrix = zeros(num_zone,num_zone,num_training);
    model.randomized_tau_times_xi = [];    % tau * xi + \sum{r}
    for i = 1 : num_zone
        r_added_to_each_zone = zeros(num_training, 1);
        if i == 1 r_added_to_each_zone = sum(r_matrix(1, 2 : end, :),2); 
        else if i == num_zone r_added_to_each_zone = -sum(r_matrix(1 : end-1, num_zone, :),1); 
            else
                r_added_to_each_zone = sum(r_matrix(i, i+1 : end, :), 2) - sum(r_matrix(1 : i-1, i, :), 1);
            end
        end
        model.randomized_tau_times_xi(:,i) = model.overall_training_data.tau_in(:, i) * model.var.xi(i) + r_added_to_each_zone(:);
    end
    model.aggregated_tau_in = sum(model.randomized_tau_times_xi,2);    % You can examine that this aggregated tau_in is equal to the one without randomized operatioin. 
    % optimizaton objective
    model.var.alpha_term = zeros(num_training - model_order, 1);    % alpha term in regression model
    model.var.beta_term = zeros(num_training - model_order, 1);    % beta term in regression model
    model.var.gamma_term = zeros(num_training - model_order, 1);    % gamma term in regression model
    model.var.theta_term = zeros(num_training - model_order, 1);    % theta term in regression model
    for i = 2 : model_order + 1    % alpha term calculation
        model.var.alpha_term = model.var.alpha_term + model.var.alpha(1, i) * ...
            model.aggregated_tau_in(model_order + 1 - (i - 1) : end - (i - 1), :);
    end
    for i = 1: model_order + 1    % the ith item means the (i-1)th order term in ATDM
        model.var.beta_term = model.var.beta_term + model.var.beta(1, i) * ...    % beta term calculation 
            model.overall_training_data.h_load(model_order + 1 - (i - 1) : end - (i - 1),:);
        model.var.gamma_term = model.var.gamma_term + model.var.gamma(1, i) * ...    % gamma term calculaton
            model.overall_training_data.tau_amb(model_order + 1 -(i -1) : end - (i - 1), :);
        model.var.theta_term = model.var.theta_term + model.var.theta(1, i) * ...    % theta term calculation 
            model.overall_training_data.radiation(model_order + 1 -(i -1) : end - (i - 1), :);
    end
    model.cons = [model.cons, model.var.residual ==  model.aggregated_tau_in(model_order + 1: end) - ...
        model.var.alpha_term - model.var.beta_term - model.var.gamma_term - model.var.theta_term - model.var.occ];    % residual calculation. What is more, don't forget the model.var.occ term.
    model.obj_1 = norm(model.var.residual(:), 2)^2 + lambda * norm(model.var.xi, 2)^2;    % set 2-norm residuals as the objective function
    % solve sub-problem one
    model.settings = sdpsettings('solver', mysolver, 'verbose', 0, 'savesolverinput', 1, 'savesolveroutput',1);
    model.solution_1 = optimize(model.cons, model.obj_1, model.settings);
    % record the solution
    model.var.alpha = myFun_GetValue(model.var.alpha);    % transform sdpvar into value (batch operation)
    model.obj_1 = value(model.obj_1);    % transform sdpvar into value(sinlge operation)
    %% sub-problem two
    % dada randomization for M
    u_matrix = normrnd(0, 1, num_zone, num_zone, num_zone);   % u_matrix for encryption
%     u_matrix= zeros(num_zone,num_zone,num_zone);
    for i = 1 : num_zone
        u_added_to_each_zone = zeros(num_zone, 1);
        if i == 1 u_added_to_each_zone = sum(u_matrix(1, 2 : end, :),2); 
        else if i == num_zone u_added_to_each_zone = -sum(u_matrix(1 : end-1, num_zone, :),1); 
            else
                u_added_to_each_zone = sum(u_matrix(i, i+1 : end, :), 2) - sum(u_matrix(1 : i-1, i, :), 1);
            end
        end
       aggregated_M(i, :) = M(i, :) + u_added_to_each_zone(:)';
    end
    model.aggregated_M = sum(aggregated_M, 1);    % get the aggregated_M
%     model.aggregated_M = sum(M,1);
    % data randomization for tau_in_hat_bar
    model.tau_in_hat = model.overall_training_data .tau_in( model_order + 1 : num_training, :);    % the definition of tau_in_hat 
    for i = 1 : model_order    % calculate tau_in_hat
        model.tau_in_hat = model.tau_in_hat - model.var.alpha(i+1) * ...
            model.overall_training_data.tau_in( model_order + 1 - i: num_training - i, :);
    end
    model.tau_in_hat_total=[model.tau_in_hat_total;model.tau_in_hat];
    model.aggregated_tau_in_hat_bar = zeros (num_training - model_order , num_zone);    % define aggregated tau_in_hat_bar
    N_matrix = normrnd(0 ,1, num_zone, num_zone, num_training - model_order, num_zone);    % N_matrix for encrption
%     N_matrix = zeros(num_zone,num_zone,num_training-model_order,num_zone);
    for i = 1 : num_zone
        N_added_to_each_zone = zeros(num_training - model_order, num_zone);
        if i == 1 N_added_to_each_zone = sum(N_matrix(1, 2 : end, :, :),2); 
        else if i == num_zone N_added_to_each_zone = -sum(N_matrix(1 : end-1, num_zone, :, :),1); 
            else
                N_added_to_each_zone = sum(N_matrix(i, i+1 : end, :, :), 2) - sum(N_matrix(1 : i-1, i, :, :), 1);
            end
        end
        model.aggregated_tau_in_hat_bar = model.aggregated_tau_in_hat_bar + ...
            model.tau_in_hat(:, i) * M(i, :) + squeeze(N_added_to_each_zone);
    end
    % variable
    % model.var.alpha directly takes the value of model.var.alpha (after myFun_GetValue operation)
    model.var.xi_bar = sdpvar(num_zone, 1, 'full');    % variable xi after transformation
    model.var.beta = sdpvar(1, model_order+1, 'full');
    model.var.gamma = sdpvar(1, model_order+1, 'full');
    model.var.theta = sdpvar(1, model_order+1, 'full');
    model.var.occ = sdpvar(num_training-model_order, 1, 'full');
    model.var.residual = sdpvar(num_training-model_order, 1, 'full');
    % constraints
    model.cons = []; 
    model.cons = [model.cons, model.var.occ(49 : num_training - model_order, 1) == model.var.occ(1 : num_training - model_order - 48, 1)];    % constraints for periodic occupant behaviour
    model.cons = [model.cons, model.aggregated_M * model.var.xi_bar == 1];    % the sum of aggregated coefficient is 1
%     model.cons = [model.cons, D * M model.var.xi_bar >= 0 ];    % aggregated coeffient is non-negative
    % optimization objective
    model.var.alpha_term = zeros(num_training - model_order, 1);    % alpha term in regression model
    model.var.beta_term = zeros(num_training - model_order, 1);    % beta term in regression model
    model.var.gamma_term = zeros(num_training - model_order, 1);    % gamma term in regression model
    model.var.theta_term = zeros(num_training - model_order, 1);    % theta term in regressino model
    model.var.alpha_term = model.aggregated_tau_in_hat_bar * model.var.xi_bar;    % alpha term calculation
    for i = 1: model_order + 1    % the ith item means the (i-1)th order term in ATDM
        model.var.beta_term = model.var.beta_term + model.var.beta(1, i) * ...    % beta term calculation 
            model.overall_training_data.h_load(model_order + 1 - (i - 1) : end - (i - 1),:);
        model.var.gamma_term = model.var.gamma_term + model.var.gamma(1, i) * ...    % gamma term calculaton
            model.overall_training_data.tau_amb(model_order + 1 -(i -1) : end - (i - 1), :);
        model.var.theta_term = model.var.theta_term + model.var.theta(1, i) * ...    % theta term calculation 
            model.overall_training_data.radiation(model_order + 1 -(i -1) : end - (i - 1), :);
    end
    model.cons = [model.cons, model.var.residual ==  model.var.alpha_term - ...
        model.var.beta_term - model.var.gamma_term - model.var.theta_term - model.var.occ];    % residual calculation. Note that the model.var.occ term should not be ignored.
    % solve sub-problem two
    model.settings = sdpsettings('solver', mysolver, 'verbose', 0, 'savesolverinput', 1, 'savesolveroutput',1);
    model.obj_2 = norm(model.var.residual(:), 2)^2 + lambda * norm(M * model.var.xi_bar, 2)^2;    % set 2-norm residuals as the objective function
    model.solution_2 = optimize(model.cons, model.obj_2, model.settings);
    % record the solution
    model.var = myFun_GetValue(model.var);    % transform sdpvar into value(batch operation)
    model.var.xi = M * value(model.var.xi_bar);    % recover aggregated coefficient xi. This means the BLA allocates xi_bar to each zone, and each zone recover their own aggregated coefficient xi.
    model.obj_2 = value(model.obj_2);    % transform sdpvar into value(sinlge operation)
    %% iteration termination condition
    GAP = min([abs(model.obj_1 - model.obj_2), abs(model.obj_1 - model.obj_2) / model.obj_2])
    if GAP < 10e-6
        Is_convergence = true;    % When GAP < delta, the iteration converged. Set the Is_convergence flag as true.
    end    
    GAP_vector = [GAP_vector, GAP];
   end
end
