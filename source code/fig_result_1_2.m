% clear;clc;
%% data import
pp_data = load('pp_data.mat'); 
pp_model = pp_data.model;
% xi = var.xi; alpha = var.alpha; beta = var.beta; gamma = var.gamma; theta = var.theta; 
% tau_agg = pp_data.pp_model.data_all.Tau_in * xi;
% tau_pre = zeros(1440,1);
pp_model_order = 2;
% tau_pre(1:pp_model_order)=tau_agg(1:pp_model_order);
% for i =pp_model_order+1 : 1440 % 逐点预测
%     alpha_term = zeros(1,1);
%     beta_term =zeros(1,1);
%     gamma_term =zeros(1,1);
%     theta_term = zeros(1,1);
%     for j = 2:pp_model_order+1
%        alpha_term =  alpha_term + alpha(1,j) * tau_pre(i-(j-1),1);
%     end
%     for j = 1:pp_model_order+1
%         beta_term = beta_term + beta(1,j) * data.h_sum(i-(j-1),1);
%         gamma_term = gamma_term + gamma(1,j)* data.Tau_amb(i-(j-1),1);
%         theta_term = theta_term + theta(1,j)* data.radiation(i-(j-1),1);
%     end
%     tau_pre(i,1)=alpha_term + beta_term + gamma_term + var.occ(1+mod(i-1,48),1); 
% end
t=1:1440;t=t';


% tau_agg = pp_model.overall_data.tau_in*pp_model.var.xi; % 1440*1
% tau_pre(1:pp_model_order,:)=tau_agg(1:pp_model_order,:);
% for i =pp_model_order+1 : 1440 % 逐点预测
%     alpha_term = zeros(1,1);
%     beta_term =zeros(1,1);
%     gamma_term =zeros(1,1);
%     theta_term = zeros(1,1);
%     for j = 2:pp_model_order+1
%         alpha_term =  alpha_term + pp_model.var.alpha(1,j) * tau_pre(i-(j-1),1);
%     end
%     for j = 1:pp_model_order+1
%         beta_term = beta_term +pp_model.var.beta(1,j) *pp_model.overall_data.h_load(i-(j-1),1);
%         gamma_term = gamma_term+pp_model.var.gamma(1,j)*pp_model.overall_data.tau_amb(i-(j-1),1);
%         theta_term = theta_term+pp_model.var.theta(1,j)*pp_model.overall_data.radiation(i-(j-1),1);
%     end
%     tau_pre(i,1)=alpha_term+beta_term+gamma_term+theta_term+pp_model.var.occ(1+mod(i-pp_model_order,48),1);
% end

tau_ori = [pp_model.training.aggregated_tau_in_ori; pp_model.verifying.aggregated_tau_in_ori];
tau_pre = [pp_model.training.aggregated_tau_in_pre; pp_model.verifying.aggregated_tau_in_pre];

%%
h_fig = figure();          % gcf: get current figure
h_axis = gca;              % gca: get current axis
%% set position & color
% position, color,
left = 10; bottom = 10; width = 18; height = 5;
% units:inches|centimeters|normalized|points|{pixels}|characters
set(h_fig, 'Units','centimeters', 'position', [left, bottom, width, height], 'color', 'w');
%% Remove the blank edge
set(gca,'LooseInset',get(gca,'TightInset'));
%% Setting color
cmap = brewermap(8,'Set1');
h_fig.Colormap = cmap;
h_axis.Colormap = cmap;
colororder(cmap)
%% plot
% for i = 1:64
% plot(t/2,pp_model.overall_data.tau_in(:,i),'--', 'color',[0.75 0.75 0.75], 'LineWidth',0.5);
% hold on;
% end
% h0 = plot(t/2,pp_model.overall_data.tau_in,'--', 'color',[0.75 0.75 0.75], 'LineWidth',0.5);hold on;
% tau_meas = [pp_model.overall_training_data.tau_in;pp_model.overall_verifying_data.tau_in];
% h0 = plot(t/2,tau_meas,'--', 'color',[0.75 0.75 0.75], 'LineWidth',0.5);hold on;
h1 = plot(t/2,tau_ori,'-.r', 'LineWidth',2); hold on;
h2 = plot(t/2,tau_pre,'-b', 'LineWidth',2); hold on;
grid on;
set(h_axis, 'XTick', [0:50:750]);                        
axis([t(1)/2, t(end)/2, 10, 35])
%% axis
set(gca,'FontName','Times New Roman','FontSize',12);
%% Title & Label
xlabel('Time (h)','FontSize',12);
ylabel('℃');
legend([h1,h2], {'Aggregate state (real)', ...
    'Aggregate state (predicted)'}, ...
    'Orientation','horizontal', ...
    'NumColumns',2, ...
    'Location', 'north', ...
    'FontSize',12,'FontName','Times New Roman');
legend boxoff;

  % % line
    px=[0.71 0.71]; py=[0.23 0.75];
    annotation('line', px,py, 'color', [0 0 0], 'LineWidth', 0.5);
  % % arrow
    px=[0.70, 0.62]; py=[0.7 0.7];
    annotation('textarrow',px, py, 'LineWidth' ,0.5, ...
    'HeadStyle', 'cback2', 'HeadWidth', 6, 'HeadLength', 6);
  % % text
    px=[0.57, 0.57]; py=[0.11 0.11];
    annotation('textbox', [px,py], 'String', 'Training data', 'FitBoxToText', 'on', ...
    'LineStyle','none','FontSize',10,'FontName','Times New Roman');
% % arrow
    px=[0.72, 0.80]; py=[0.7 0.7];
    annotation('textarrow',px, py, 'LineWidth' ,0.5, ...
    'HeadStyle', 'cback2', 'HeadWidth', 6, 'HeadLength', 6);
% % text
    px=[0.72, 0.57]; py=[0.11 0.11];
    annotation('textbox', [px,py], 'String', 'Test data', 'FitBoxToText', 'on', ...
    'LineStyle','none','FontSize',10,'FontName','Times New Roman');



