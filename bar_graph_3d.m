clear; clc; close all;

xlsxfile = fullfile('/Users/jshine/My Drive/Python/fourier_roughness','fourier_parameters.xlsx');

A = readmatrix(xlsxfile,'Sheet','Amplitude');
phi = readmatrix(xlsxfile,'Sheet','Phase Shift');
N = readmatrix(xlsxfile,'Sheet','N');
M = readmatrix(xlsxfile,'Sheet','M');
lambda_k = readmatrix(xlsxfile,'Sheet','Wavelength');

fig = figure();
fig.Units = 'normalized';
fig.Position = [0, 0, 0.55, 0.95];

h = bar3(A,0.8);
set(h,'EdgeColor','None');
view(45, 55)
colormap('jet')
for i = 1:length(h)
     zdata = get(h(i),'Zdata');
     set(h(i),'Cdata',zdata)
end

xticks(linspace(1,2*N+1,5));
yticks(linspace(1,2*M+1,5));
% zticks(linspace(0,lambda_k/2,8));
xticklabels(string(round(linspace(-N,N,5))));
yticklabels(string(round(linspace(-M,M,5))));
xlabel('N')
ylabel('M')
zlabel('A_{n,m} [mils]')

ax = gca;
ax.FontSize = 32;
ax.FontName = 'Times';
ax.XTickLabelRotation = 0;
ax.YTickLabelRotation = 0;

saveas(fig,fullfile('/Users/jshine/My Drive/Python/fourier_roughness', 'Ogive Forebody OI Amplitude Coefficients.png'))
