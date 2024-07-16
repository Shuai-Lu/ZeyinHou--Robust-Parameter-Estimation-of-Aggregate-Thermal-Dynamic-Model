clear;clc;
%% data import
pp_data = load('pp_data.mat'); pp_var = pp_data.model.var;
npp_data =load('npp_data.mat'); npp_var = npp_data.model.var;
pp_xi = pp_var.xi;
npp_xi = npp_var.xi;
% delta_r = (pp_xi-npp_xi)./npp_xi;
delta = pp_xi-npp_xi;
k=1:64;



%% plot figure 1 
h_fig = figure(1);          % gcf: get current figure
h_axis = gca;              % gca: get current axis
%% set position & color
% position, color,
left = 10; bottom = 10; width = 18; height = 5;
% units:inches|centimeters|normalized|points|{pixels}|characters
set(h_fig, 'Units','centimeters', 'position', [left, bottom, width, height], 'color', 'w');
%% Remove the blank edge
% set(gca,'LooseInset',get(gca,'TightInset'));
%% Setting color
cmap = brewermap(8,'Set1');
h_fig.Colormap = cmap;
h_axis.Colormap = cmap;
colororder(cmap)
%% plot
set(gca,'LooseInset',[0.05,0.05,0.05,0.05]);
h_line(1) = plot(k,pp_xi,'o'); hold on;
h_line(2) = plot(k,npp_xi,'*');
grid on;
%% axis
set(gca,'FontName','Times New Roman','FontSize',12);
%% Title & Label
% title('Aggregation coefficient results','FontSize',12);
xlabel('Zone No.','FontSize',12);
ylabel('\xi');
ylim([min(pp_xi)-0.005,max(pp_xi)+0.015]);
legend(h_line(1:2), {'privacy-preser.', ...
    'non-privacy-preser.'}, ...
    'Orientation','horizontal', ...
    'NumColumns',2, ...
    'Location', 'north', ...
    'FontSize',12,'FontName','Times New Roman');
legend boxoff;
set(gca,'TickLabelInterpreter','latex');


%% plot figure 2
h_fig = figure(2);          % gcf: get current figure
h_axis = gca;              % gca: get current axis
%% set position & color
% position, color,
left = 10; bottom = 10; width = 18; height = 5;
% units:inches|centimeters|normalized|points|{pixels}|characters
set(h_fig, 'Units','centimeters', 'position', [left, bottom, width, height], 'color', 'w');
%% Remove the blank edge
% set(gca,'LooseInset',get(gca,'TightInset'));
%% Setting color
cmap = brewermap(8,'Set1');
h_fig.Colormap = cmap;
h_axis.Colormap = cmap;
colororder(cmap)
set(gca,'LooseInset',[0.05,0.05,0.05,0.05]);
h_line(3) = plot(k',(delta),'*'); hold on;
grid on;
% axis
set(gca,'FontName','Times New Roman','FontSize',12);
% Title & Label
xlabel('Zone No.','FontSize',12);
ylabel('%','FontSize',12);
ylim([-1*10e-9,1*10e-9]);


%% plot figure 3
pp_data = load('pp_data.mat'); 
pp_model = pp_data.model;
pp_model_order = 2;
t=1:1440;t=t';
tau_agg = pp_model.overall_data.tau_in*pp_model.var.xi; % 1440*1
tau_pre(1:pp_model_order,:)=tau_agg(1:pp_model_order,:);
for i =pp_model_order+1 : 1440 % 逐点预测
    alpha_term = zeros(1,1);
    beta_term =zeros(1,1);
    gamma_term =zeros(1,1);
    theta_term = zeros(1,1);
    for j = 2:pp_model_order+1
        alpha_term =  alpha_term + pp_model.var.alpha(1,j) * tau_pre(i-(j-1),1);
    end
    for j = 1:pp_model_order+1
        beta_term = beta_term +pp_model.var.beta(1,j) *pp_model.overall_data.h_load(i-(j-1),1);
        gamma_term = gamma_term+pp_model.var.gamma(1,j)*pp_model.overall_data.tau_amb(i-(j-1),1);
        theta_term = theta_term+pp_model.var.theta(1,j)*pp_model.overall_data.radiation(i-(j-1),1);
    end
    tau_pre(i,1)=alpha_term+beta_term+gamma_term+theta_term+pp_model.var.occ(1+mod(i-pp_model_order,48),1);
end
h_fig = figure(3);          % gcf: get current figure
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
h0 = plot(t/2,pp_model.overall_data.tau_in,'--', 'color',[0.75 0.75 0.75], 'LineWidth',0.5);hold on;
h1 = plot(t/2,tau_agg,'-.r', 'LineWidth',2); hold on;
h2 = plot(t/2,tau_pre,'-b', 'LineWidth',2); hold on;
grid on;
set(h_axis, 'XTick', [0:50:750]);                        
axis([t(1)/2, t(end)/2, 10, 35])
%% axis
set(gca,'FontName','Times New Roman','FontSize',12);
%% Title & Label
xlabel('Time(h)','FontSize',12);
ylabel('℃');
legend([h1,h2,h0(1)], {'Aggregate state (real)', ...
    'Aggregate state (predicted)','Real indoor temperatures of zones'}, ...
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
    px=[0.72, 0.80]; py=[0.7 0.7];
    annotation('textarrow',px, py, 'LineWidth' ,0.5, ...
    'HeadStyle', 'cback2', 'HeadWidth', 6, 'HeadLength', 6);
% % text
    px=[0.72, 0.57]; py=[0.11 0.11];
    annotation('textbox', [px,py], 'String', 'Test data', 'FitBoxToText', 'on', ...
    'LineStyle','none','FontSize',10,'FontName','Times New Roman');

disp('Figure 1 is the comparison of aggregation coefficient values (privacy-preserved vs. non-privacy-preserved algorithms).');
disp('Figure 2 is aggregation errors between privacy-preserved and non-privacy-preserved algorithm.');
disp('Figure 3 desribes the aggregation state and real indoor temperatures of all zones.');


%% display parameter values
disp(['------------------------------------------------------------------------']);
disp(['The parameter estimation results through the privacy-preserved algorithm']);
disp(['------------------------------------------------------------------------']);
disp(['alpha:    ',num2str(pp_var.alpha)]);
disp(['beta:    ',num2str(pp_var.beta)]);
disp(['gamma:   ',num2str(pp_var.gamma)]);
disp(['theta:   ',num2str(pp_var.theta)]);
disp(['------------------------------------------------------------------------']);
disp(['The parameter estimation results through the non-privacy-preserved algorithm']);
disp(['------------------------------------------------------------------------']);
disp(['alpha:    ',num2str(npp_var.alpha)]);
disp(['beta:    ',num2str(npp_var.beta)]);
disp(['gamma:   ',num2str(npp_var.gamma)]);
disp(['theta:   ',num2str(npp_var.theta)]);
disp(['------------------------------------------------------------------------']);
disp(['The computation results of indicators through the privacy-preserved algorithm']);
disp(['------------------------------------------------------------------------']);
disp(['RMSE:    ',num2str(pp_data.model.verifying.rmse)]);
disp(['MAPE:    ',num2str(pp_data.model.verifying.mpae)]);
disp(['R^2:     ',num2str(pp_data.model.verifying.r2)]);
disp(['------------------------------------------------------------------------']);
disp(['The computation results of indicators through the non-privacy-preserved algorithm']);
disp(['------------------------------------------------------------------------']);
disp(['RMSE:    ',num2str(npp_data.model.verifying.rmse)]);
disp(['MAPE:    ',num2str(npp_data.model.verifying.mpae)]);
disp(['R^2:     ',num2str(npp_data.model.verifying.r2)]);
disp(['------------------------------------------------------------------------']);




