clc; clear; close all; format compact;
dataName = 'Yang7'; 
[option, para1, para2, para3, para4, para5] = OptSetting(dataName);
runtest(dataName, @mu1, @dmu1, para1, option);
runtest(dataName, @mu2, @dmu2, para2, option);
runtest(dataName, @mu3, @dmu3, para3, option);
runtest(dataName, @mu4, @dmu4, para4, option);
runtest(dataName, @mu5, @dmu5, para5, option);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- Sub-functions  ------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = mu1(para, x)
% definition of the description mu_1 with
% three parameters mu_max, alpha and x_opt.
% input:
%	para : vector
%	parameters in the description
%	x : vector
% 	measured irradiance
mumax = para(1); al = para(2); xopt = para(3);
y = mumax*x ./ ( x + mumax/al * (x./xopt - 1).^2 );
end
function y = dmu1(para, x)
% definition of the partial derivative of
% mu_1 w.r.t. each parameters.
% input:
%	para : vector
%	parameters in the description
%	x : vector
% 	measured irradiance
mumax = para(1); al = para(2); xopt = para(3);
dmumax = al^2*xopt^4*x.^2 ./ ...
    ( mumax*(xopt-x).^2 + al*xopt^2*x ).^2;
dal = (mumax*xopt)^2*x .* (xopt-x).^2 ./ ...
    ( mumax*(xopt-x).^2 + al*xopt^2*x ).^2;
dxopt = -2*mumax^2*al*xopt*x.^2 .* (xopt-x) ./ ...
    ( mumax*(xopt-x).^2 + al*xopt^2*x ).^2;
y = [dmumax, dal, dxopt];
end

function y = mu2(para, x)
% definition of the description mu_2 with
% three parameters mu^bar, K_x and K_i.
% input:
%	para : vector
%	parameters in the description
%	x : vector
% 	measured irradiance
mubar = para(1); Kx = para(2); Ki = para(3);
y = mubar*x ./ ( x + Kx + x.^2/Ki );
end
function y = dmu2(para, x)
% definition of the partial derivative of
% mu_2 w.r.t. each parameters.
% input:
%	para : vector
%	parameters in the description
%	x : vector
% 	measured irradiance
mubar = para(1); Kx = para(2); Ki = para(3);
dmubar = x ./ ( x + Kx + x.^2/Ki );
dKx = -mubar*x ./ ( x + Kx + x.^2/Ki ).^2;
dKi = mubar*x.^3 ./ ( Ki*x + Ki*Kx + x.^2 ).^2;
y = [dmubar, dKx, dKi];
end

function y = mu3(para, x)
% definition of the description mu_3 with
% three parameters a, b and c.
% input:
%	para : vector
%	parameters in the description
%	x : vector
% 	measured irradiance
a = para(1); b = para(2); c = para(3);
y = x ./ ( a*x.^2 + b*x + c );
end
function y = dmu3(para, x)
% definition of the partial derivative of
% mu_3 w.r.t. each parameters.
% input:
%	para : vector
%	parameters in the description
%	x : vector
% 	measured irradiance
a = para(1); b = para(2); c = para(3);
da = -x.^3 ./ (a*x.^2 + b*x +c ).^2;
db = -x.^2 ./ (a*x.^2 + b*x +c ).^2;
dc = -x ./ (a*x.^2 + b*x +c ).^2;
y = [da, db, dc];
end

function y = mu4(para, x)
% definition of the description mu_4 with
% three parameters mu^bar, K_x and K_i.
% input:
%	para : vector
%	parameters in the description
%	x : vector
% 	measured irradiance
mubar = para(1); Kx = para(2); Ki = para(3);
y = mubar*Ki*x ./ (x+Kx) ./ (x+Ki);
end
function y = dmu4(para, x)
% definition of the partial derivative of
% mu_4 w.r.t. each parameters.
% input:
%	para : vector
%	parameters in the description
%	x : vector
% 	measured irradiance
mubar = para(1); Kx = para(2); Ki = para(3);
dmubar = Ki*x ./ (x+Kx) ./ (x+Ki);
dKx = -mubar*Ki*x ./ (x+Kx).^2 ./ (x+Ki);
dKi = mubar*x.^2 ./ (x+Kx) ./ (x+Ki).^2;
y = [dmubar, dKx, dKi];
end

function y = mu5(para, x)
% definition of the description mu_5 with
% two parameters gamma_max and x^star. 
% input:
%	para : vector
%	parameters in the description
%	x : vector
% 	measured irradiance
gammax = para(1); xstar = para(2);
y = 4*gammax*( xstar*x ./ (x+xstar).^2 );
end
function y = dmu5(para, x)
% definition of the partial derivative of
% mu_5 w.r.t. each parameters.
% input:
%	para : vector
%	parameters in the description
%	x : vector
% 	measured irradiance
gammax = para(1); xstar = para(2);
dgammax = 4*xstar*x ./ (x+xstar).^2;
dxstar = 4*gammax*x .* (xstar-x) ./ (x+xstar).^3;
y = [dgammax, dxstar];
end

function SSE = computeSSE(mu, x, y)
% compute SSE value for one description mu
% input:
%	mu : function
%	one of the description
%	x : vector
% 	measured irradiance
%	y : vector
%	measured growth rate
SSE = sum((mu(x) - y).^2);
end

function AICc = computeAICc(x, SSE, para)
% compute AICc value for one description mu
% input:
%	x : vector
% 	measured irradiance
%	SSE : float
%	SSE value for mu
%	para : vector
%	optimal parameters of mu
AICc = length(x)*log(SSE/length(x)) + 2*(length(para)+1) ...
    + 2*(length(para)+1)*(length(para)+2) ...
    / (length(x)-length(para)-2);
end

function BIC = computeBIC(x, SSE, para)
% compute BIC value for one description mu
% input:
%	x : vector
% 	measured irradiance
%	SSE : float
%	SSE value for mu
%	para : vector
%	optimal parameters of mu
BIC = length(x)*log(SSE/length(x)) + (length(para)+1)*log(length(x));
end

function PEMAC = computePEMAC(x, SSE, para, e)
% compute PEMAC value for one description mu
% input:
%	x : vector
% 	measured irradiance
%	SSE : float
%	SSE value for mu
%	para : vector
%	optimal parameters of mu
%	e : float
%	error propagation of mu
PEMAC = 2*length(x)*log(mean(e)) + length(x)*log(SSE/length(x)) ...
    + 2*(length(para)+1)*(length(para)+2) ...
    / (length(x)-length(para)-2);
end

function [FIM, s] = computeFIM(x, y, mu, dmu, para)
% compute Fisher information matrix and the
% residual mean square for one description mu
% input:
%	x : vector
% 	measured irradiance
%	y : vector
% 	measured growth rate
%	mu : function
%	one of the description
%	dmu : function
%	partial derivative of mu
%	para : vector
%	optimal parameters of mu
Q = diag(max( max(y)/50, 0.1*y ).^2)^-1;
FIM = dmu(para, x)'*Q*dmu(para, x); 
s = sqrt( ( (mu(para, x)-y)'*Q*(mu(para, x)-y) ) / ...
    ( length(x)-length(para) ) );
end

function [e, sgm, R] = computeSensitivity(x, dmu, para, FIM, s, z) 
% compute error propagation, standard derivation
% and the correlation matrix for one description
% input:
%	x : vector
% 	measured irradiance
%	dmu : function
%	partial derivative of mu
%	para : vector
%	optimal parameters of mu
%	FIM : matrix
%	Fisher information matrix
%	s : float
%	residual mean square
%	z : vector
%	variable in error propagation
V = FIM^-1; sd = sqrt(diag(V)); R = V./(sd*sd');
sgm = s*sqrt( length(x)*diag(V) );
e = sqrt(sum( dmu(para, z).^2.*(sgm.^2)', 2 ));
end

function plotPIcurve(varname, funcInfo, x, y, mu, dmu, para, FIM, s)
% plot the prediction interval of one description
% input:
%	varname : string
%	name of the data
%	x : vector
% 	measured irradiance
%	y : vector
% 	measured growth rate
%	mu : function
%	one of the description
%	dmu : function
%	partial derivative of mu
%	para : vector
%	optimal parameters of mu
%	FIM : matrix
%	Fisher information matrix
%	s : float
%	residual mean square
xt = 1:1:1600; %grid for illustration
sh = 3600; %covert second to hour
% sd = 3600*24; %convert second to day
e = computeSensitivity(x, dmu, para, FIM, s, xt');
%significant level value
if length(x) == 12
    sgflvl = 2.201;
elseif length(x) == 13
    sgflvl = 2.179;
elseif length(x) == 22
    sgflvl = 2.080;
end
figure
plot(x, y*sh, '+', xt,  mu(para, xt)*sh, '-', ...
    xt, ( mu(para, xt) - sgflvl*e' )*sh, '--', ...
    xt, ( mu(para, xt) + sgflvl*e' )*sh, '--', ...
    'LineWidth',2, 'MarkerSize', 10);
xlabel('$x$ ($$\mu$$mol photons m$$^{-2} $$s$$^{-1}$$)', ...
    'interpreter', 'latex'); xlim([-2 1600]);
% set ylabel for different data
if strcmp(varname, 'Anning_PcHL')
    ylabel('P$^C$ (h$$^{-1}$$)', 'interpreter', 'latex');
else %data Yang
    ylabel('P$_n$ ($$\mu$$mol O$_2$ mg$$^{-1} $$Chla h$$^{-1}$$)', ...
    'interpreter', 'latex');
end
legend('Mesurement', 'Estimation', 'Location', 'SouthEast');
set(gca, 'linewidth', 2); set(gca, 'FontSize', 20);
% save figure
figname = strcat('PI_', varname, '_', funcInfo.function, '.png');
if isfile(['fig/' figname]) == 0
    saveas(gca, ['fig/' figname]);
    disp(['Figure saved as: ', figname]);
end
end

function plotSensitivityCurve(varname, funcInfo, x, dmu, para, FIM, s)
% plot sensitivity curve of one description
% input:
%	varname : string
%	name of the data
%	x : vector
% 	measured irradiance
%	dmu : function
%	partial derivative of mu
%	para : vector
%	optimal parameters of mu
%	FIM : matrix
%	Fisher information matrix
%	s : float
%	residual mean square
xt = 1:1:1600; %grid for illustration
[e, sgm] = computeSensitivity(x, dmu, para, FIM, s, xt');
v = dmu(para, xt');
figure
plot(xt, (sgm'.*v./e).^2, ...
    'LineWidth',2, 'MarkerSize', 10);
xlabel('$x$ ($$\mu$$mol photons m$$^{-2} $$s$$^{-1}$$)', ...
    'interpreter', 'latex'); xlim([0 1600]); ylim([0 1]);
% set legend for different descriptions
if strcmp(funcInfo.function, 'mu1')
    legend({'$$ \mu_{\max}$$', '$$ \alpha$$', '$$ x_{opt}$$'}, ...
        'Location', 'southeast', 'interpreter', 'latex');
elseif strcmp(funcInfo.function, 'mu2')
    legend({'$$ \bar{\mu}$$', '$$ K_x$$', '$$ K_i$$'}, ...
        'Location', 'southeast', 'interpreter', 'latex');
elseif strcmp(funcInfo.function, 'mu3')
    legend({'$$ a$$', '$$ b$$', '$$ c$$'}, ...
       'Location', 'southeast', 'interpreter', 'latex');
elseif strcmp(funcInfo.function, 'mu4')
    legend({'$$ \bar{\mu}$$', '$$ K_x$$', '$$ K_i$$'}, ...
        'Location', 'southeast', 'interpreter', 'latex');
elseif strcmp(funcInfo.function, 'mu5')
    legend({'$$ \gamma_{\max}$$', '$$ x^{\star}$$'}, ...
        'Location', 'southeast', 'interpreter', 'latex');
end
set(gca, 'linewidth', 2); set(gca, 'FontSize', 20);
% save figure
figname = strcat('Sensitivity_', varname, '_', funcInfo.function, '.png');
if isfile(['fig/' figname]) == 0
    saveas(gca, ['fig/' figname]);
    disp(['Figure saved as: ', figname]);
end
end

function [x, y] = readData(varname)
% In Anning, they considered net growth rate
% whereas they considered growth rate in Yang, 
% thus the first meusurement is negative, and
% we convert it to net growth rate
% input:
%	varname : string
%	name of the data
originalPath = pwd; cd("data/");
filename = [varname '.txt']; 
data = readmatrix(filename); cd(originalPath);
x = data(:, 1);
%check whether net growth rate in data 
if data(1, 2) < 0 
    y = (data(:, 2) - data(1, 2))/3600; %convert hour to second
else
    y = data(:, 2)/3600; %convert hour to second
end
end

function [option, para1, para2, para3, para4, para5] = OptSetting(varname)
% Some results are very sensitive to the initial
% value in the optimization, as we are searching
% for reasonable local optimum
% input:
%	varname : string
%	name of the data
option = optimset('TolFun', 1e-12, 'TolX', 1e-12, ...
    'MaxFunEvals', 1e10, 'MaxIter', 1e8);
if strcmp(varname, 'Yang5') %data Yang5
    para1 = [1, 1, 1]; para2 = [-10, -10, -1]; 
    para3 = [0, -1, 1000]; para4 = [100, 10, 10];
elseif strcmp(varname, 'Anning_PcHL') %data Anning_PcHL
    para1 = [0, 0, 100]; para2 = [0, 1, 1]; 
    para3 = [0, 0, 1]; para4 = [1, 1, 1];
else %data Yang1-4, Yang6-7
    para1 = [0, 1, 1]; para2 = [0, 1, 1]; 
    para3 = [0, 0, 0]; para4 = [1, 1, 1];
end
para5 = [0, 2];
end

function runtest(varname, mu, dmu, ParaInit, option)
% run all requested tests in the paper
% input:
%	varname : string
%	name of the data
%	mu : function
%	one of the description
%	dmu : function
%	partial derivative of mu
%	ParaInit : vector
%	initial value for optimization
[x, y] = readData(varname); funcInfo = functions(mu);
disp(['Run the test for dataset ', varname, ...
    ' with description ', funcInfo.function]);
% minimize the value of SSE to find best parameter values
[ParaOpt, SSE] = fminsearch(@(para)computeSSE(@(x)mu(para, x), ...
    x, y), ParaInit, option);
disp('Optimal parameter values :')
fprintf(1, '%.10e\n', ParaOpt); fprintf(1, 'SSE = %.4e\n', SSE);
[FIM, s] = computeFIM(x, y, mu, dmu, ParaOpt);
disp('The correlation matrix is')
if strcmp(funcInfo.function, 'mu4') && abs(det(FIM)) < 1e-12
    PEMAC = NaN;
    fprintf(2, 'Fisher information matrix is singular !!\n');
else
    [e, ~, R] = computeSensitivity(x, dmu, ParaOpt, FIM, s, x);
    disp(R);
    PEMAC = computePEMAC(x, SSE, ParaOpt, e);
    plotPIcurve(varname, funcInfo, x, y, mu, dmu, ParaOpt, FIM, s)
    plotSensitivityCurve(varname, funcInfo, x, dmu, ParaOpt, FIM, s)
end
AICc = computeAICc(x, SSE, ParaOpt);
BIC = computeBIC(x, SSE, ParaOpt);
fprintf(1, 'AICc= %.4f |BIC= %.4f |PEMAC= %.4f\n', AICc, BIC, PEMAC);
disp('---------------------------')
end