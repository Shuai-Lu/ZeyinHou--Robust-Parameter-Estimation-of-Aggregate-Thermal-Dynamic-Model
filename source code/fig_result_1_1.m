clear;clc;
%% data import
pp_data = load('pp_data.mat'); pp_var = pp_data.model.var;
npp_data =load('npp_data.mat'); npp_var = npp_data.model.var;
pp_xi = pp_var.xi;
npp_xi = npp_var.xi;
% delta_r = (pp_xi-npp_xi)./npp_xi;
delta = pp_xi-npp_xi;
k=1:64;
%%
h_fig = figure(); % gcf: get current figure
h_axis = gca; % gca: get current axis
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
legend(h_line(1:2), {'PPHME.', ...
'HME.'}, ...
'Orientation','horizontal', ...
'NumColumns',2, ...
'Location', 'northeast', ...
'FontSize',12,'FontName','Times New Roman');
legend boxoff;
set(gca,'TickLabelInterpreter','latex');

%% plot
% set(gca,'LooseInset',[0.05,0.05,0.05,0.05]);
% h_line(3) = plot(k',(delta),'*'); hold on;
% grid on;
% % axis
% set(gca,'FontName','Times New Roman','FontSize',12);
% % Title & Label
% xlabel('Zone No.','FontSize',12);
% ylabel('%','FontSize',12);
% ylim([-1*10e-7,1*10e-7]);