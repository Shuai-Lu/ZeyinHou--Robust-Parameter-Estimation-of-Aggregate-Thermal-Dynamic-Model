clear; clc; 
%% global variable definition
global model lambda GAP_vector   % Global variable should also be defined in the sub-function.
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
building_id = [1 8 10 11 13 16]';     % building ID
% building_id = [1 8 10 11 13 16 20]';     % building ID
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
%% Privacy preservation process
W = normrnd(0.1,0.01,num_zone,num_zone);
% tau_in times W^T
p_matrix = normrnd(0, 1, num_zone, num_zone, num_training, num_zone); 
    model.tau_in_times_W_summation = zeros(num_training,num_zone); 
    for i = 1 : num_zone
        p_added_to_each_zone = zeros(num_training, num_zone);
        if i == 1 
            p_added_to_each_zone = sum(p_matrix(1, 2 : end, : , :),2); 
        elseif i == num_zone 
            p_added_to_each_zone = -sum(p_matrix(1 : end-1, num_zone, :, :),1); 
        else
                p_added_to_each_zone = sum(p_matrix(i, i+1 : end, :, :), 2) ...
                    - sum(p_matrix(1 : i-1, i, :, :), 1);
        end
        model.tau_in_times_W_summation = model.tau_in_times_W_summation +...
            model.overall_training_data.tau_in(:,i) * W(:,i)' +...
            squeeze(p_added_to_each_zone);
    end
    % WW^T
q_matrix = normrnd(0, 1, num_zone, num_zone, num_zone,num_zone); 
model.W_times_W_trans_summation = zeros(num_zone,num_zone);
 for i = 1 : num_zone
        q_added_to_each_zone = zeros(num_zone, num_zone);
        if i == 1 
            q_added_to_each_zone = sum(q_matrix(1, 2 : end, : , :),2); 
        elseif i == num_zone 
            q_added_to_each_zone = -sum(q_matrix(1 : end-1, num_zone, :, :),1); 
        else
                q_added_to_each_zone = sum(q_matrix(i, i+1 : end, :, :), 2) ...
                    - sum(q_matrix(1 : i-1, i, :, :), 1);
        end
        model.W_times_W_trans_summation = model.W_times_W_trans_summation +...
            W(:,i) * W(:,i)' +...
            squeeze(q_added_to_each_zone);
 end
 % 1^T times W^T
r_matrix = normrnd(0, 1, num_zone, num_zone, num_zone,1); 
model.ones_times_W_summation = zeros(num_zone,1);
 for i = 1 : num_zone
        r_added_to_each_zone = zeros(num_zone, 1);
        if i == 1 
            r_added_to_each_zone = sum(r_matrix(1, 2 : end, : , :),2); 
        elseif i == num_zone 
            r_added_to_each_zone = -sum(r_matrix(1 : end-1, num_zone, :, :),1); 
        else
                r_added_to_each_zone = sum(r_matrix(i, i+1 : end, :, :), 2) ...
                    - sum(r_matrix(1 : i-1, i, :, :), 1);
        end
        model.ones_times_W_summation = model.ones_times_W_summation +...
            W(:,i) +...
            squeeze(r_added_to_each_zone);
 end
%% training data disturbance
noise_ratio = load("noise_ratio_for_tau_in.mat");
noise_ratio = noise_ratio.ratio_m;
model.overall_training_data.tau_in  = model.overall_training_data.tau_in .* noise_ratio(:,1:num_zone);
model.overall_training_data.h_load = model.overall_training_data.h_load  .* noise_ratio(:,end-2);
model.overall_training_data.tau_amb = model.overall_training_data.tau_amb .* noise_ratio(:,end-1);
model.overall_training_data.radiation = model.overall_training_data.radiation .* noise_ratio(:,end);
%% Perform the HME
% 计算初始残差r_0
 tic;
 omega = ones(num_training-model_order,1);
 BCD(model_order, num_training, num_zone,omega);
 r_0 = model.var.residual; 
 k=0;
 delta_HME = 0.7;
 r=r_0;
 thre = [];
 while k<1000
     omega = ones(num_training-model_order,1);
     for t = 1:num_training-model_order
         if abs(r(t,1))<=delta_HME
             omega(t,1)= 1;
         else
             omega(t,1) = delta_HME/abs(r(t,1));
         end
     end
     BCD(model_order, num_training, num_zone,omega);
     r_new = model.var.residual;
     thre = [thre, norm(r_new-r,2)^2];
     if thre(end) < 10^-6
         break;
     end
     r = r_new;
     k=k+1;
 end
 toc;


%% Indicator calculation
model.training.aggregated_tau_in_ori = model.overall_training_data.tau_in*model.var.xi; % 1080*1
model.verifying.aggregated_tau_in_ori = model.overall_verifying_data.tau_in*model.var.xi; % 360*1
model.verifying.aggregated_tau_in_pre = zeros(num_verifying,1); % 360*1
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

 save('npp_data.mat','model');
 disp(['operating time:', num2str(toc)]);
 





%% Privacy Preserved Block Coordinated Descent Method (PP-BCDM) for Aggregated Thermal Dynamic Model (ATDM)
function BCD (model_order, num_training,num_zone,omega)
    %% global variable definition
    global model lambda GAP_vector
    %% predefition & preprocess
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
    model.aggregated_tau_in = model.overall_training_data.tau_in * model.var.xi;
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
    model.obj_1 = norm(model.var.residual(:).*sqrt(omega), 2)^2 + lambda * norm(model.var.xi, 2)^2;    % set 2-norm residuals as the objective function
    % solve sub-problem one
    model.settings = sdpsettings('solver', mysolver, 'verbose', 0, 'savesolverinput', 1, 'savesolveroutput',1);
    model.solution_1 = optimize(model.cons, model.obj_1, model.settings);
    % record the solution
    model.var.alpha = myFun_GetValue(model.var.alpha);    % transform sdpvar into value (batch operation)
    model.obj_1 = value(model.obj_1);    % transform sdpvar into value(sinlge operation)
    %% sub-problem two
    model.tau_in_hat = model.overall_training_data.tau_in( model_order + 1 : num_training, :);    % the definition of tau_in_hat 
    for i = 1 : model_order    % calculate tau_in_hat
        model.tau_in_hat = model.tau_in_hat - model.var.alpha(i+1) * ...
            model.overall_training_data.tau_in( model_order + 1 - i: num_training - i, :);
    end
    % variable
    % model.var.alpha directly takes the value of model.var.alpha (after myFun_GetValue operation)
    model.var.xi = sdpvar(num_zone, 1, 'full');    % variable xi after transformation
    model.var.beta = sdpvar(1, model_order+1, 'full');
    model.var.gamma = sdpvar(1, model_order+1, 'full');
    model.var.theta = sdpvar(1, model_order+1, 'full');
    model.var.occ = sdpvar(num_training-model_order, 1, 'full');
    model.var.residual = sdpvar(num_training-model_order, 1, 'full');
    % constraints
    model.cons = []; 
    model.cons = [model.cons, model.var.occ(49 : num_training - model_order, 1) == model.var.occ(1 : num_training - model_order - 48, 1)];    % constraints for periodic occupant behaviour
    model.cons = [model.cons, sum(model.var.xi) == 1];    % the sum of aggregated coefficient is 1
    % optimization objective
    model.var.alpha_term = zeros(num_training - model_order, 1);    % alpha term in regression model
    model.var.beta_term = zeros(num_training - model_order, 1);    % beta term in regression model
    model.var.gamma_term = zeros(num_training - model_order, 1);    % gamma term in regression model
    model.var.theta_term = zeros(num_training - model_order, 1);    % theta term in regressino model
    model.var.alpha_term = model.tau_in_hat * model.var.xi;    % alpha term calculation
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
    model.obj_2 = norm(model.var.residual(:).*sqrt(omega), 2)^2 + lambda * norm(model.var.xi, 2)^2;    % set 2-norm residuals as the objective function
    model.solution_2 = optimize(model.cons, model.obj_2, model.settings);
    % record the solution
    model.var = myFun_GetValue(model.var);    % transform sdpvar into value(batch operation)
    model.var.xi = value(model.var.xi);    % recover aggregated coefficient xi. This means the BLA allocates xi_bar to each zone, and each zone recover their own aggregated coefficient xi.
    model.obj_2 = value(model.obj_2);    % transform sdpvar into value(sinlge operation)
    %% iteration termination condition
    GAP = min([abs(model.obj_1 - model.obj_2), abs(model.obj_1 - model.obj_2) / model.obj_2])
    if GAP < 10e-6
        Is_convergence = true;    % When GAP < delta, the iteration converged. Set the Is_convergence flag as true.
    end    
    GAP_vector = [GAP_vector, GAP];
   end
end
