%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelling post-stack data from the reflectivity series             %%%%
% d = s(t) * r(t) = G r(t)                                           %%%%
%--------------------------------------------------------------------%%%%
% Last update: May, 2020                                             %%%%
% Version: 0.0                                                       %%%%
% Any problem please report to: sergioluiz.pesquisa@gmail.com        %%%%
% Product: sergioluizedu/MbPsiRefEnt2020                             %%%%
%                                                                    %%%%
% Federal University of Rio Grande do Norte (UFRN)                   %%%%
% Programa de Pós-Graduação em Física (PPGF-UFRN)                    %%%%
%--------------------------------------------------------------------%%%%
% Copyright (C) 2020 Sérgio Luiz Eduardo (sergioluizufrn@ufrn.edu.br)%%%%
% This program is free software: you can redistribute it and/or modify%%%
% it under the terms of the GNU General Public License as published  %%%%
% by the Free Software Foundation, either version 3 of the License,  %%%%
% or (at your option) any later version.                             %%%%
%                                                                    %%%%
% This program is distributed in the hope that it will be useful, but%%%%
% WITHOUT ANY WARRANTY; without even the implied warranty of         %%%%
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.               %%%%
%                                                                    %%%%
% See the GNU General Public License for more details.               %%%%
%                                                                    %%%%
% See http://www.gnu.org/licenses/.                                  %%%% 
%--------------------------------------------------------------------%%%%
% Purpose:                                                           %%%%
% This code performs the modelling of Post-Stack data from a         %%%%
% reflectivity series model.                                         %%%%
%--------------------------------------------------------------------%%%%

%% setup
clear all; clc; close all; addpath('./core/');

%% Acquisition parameters
freqPICO=55;                % Peak Ricker wavelet frequency
snr=30;                     % Gaussian noise with a signal-to-noise ratio (SNR)
nt=2001;                    % Time samples
nx=1701;                    % Lateral distance samples (x-axis)
nz=1751;                    % Depth samples (z-axis)
dt=1e-3;                    % t-Discretization in sec
dx=10.;                     % x-Discretization in meters
dz=2.;                      % z-Discretization in meters
janT=1:1801; janX=100:500;  % Parameters for cutting the Marmousi2 model

%% Impedance model in time-depth domain 
fileID = fopen('./models/Marmousi2Ztime_sergioMod_dt_.001_dx_10_nt_2001_nx_1701.bin','r');
Z=fread(fileID,[nt nx],'float32');  % Read model
fclose(fileID);
Z=Z(janT,janX);                     % Cut model
Z=resizem(Z,[501 length(janX)]);    % Resample model

[nt,nx]=size(Z);                    % New dimensions
x=0:dx:(nx-1)*dx; t=(0:nt-1)*dt;    % Spatial and temporal coordinates
%% Obtaining the reflectivity model from the acoustic impedance model
% r = 0.5 d/dt ln(Z)
ref=.5.*[diff(log(Z)); zeros(1,nx)];

%% Seismic Source: Ricker wavelet
s=ricker(freqPICO,dt,nt);

%% Build operator  convolution (two-dimensional)
G = convmtx2(s,nt,1); G=G(1:nt,:);

%% Observed data
dobs=G*ref;                 % noiseless data
dobs_n=awgn(dobs,snr);      % noisy data


%% save observed data and acquisition parameters
mkdir data
dlmwrite('./data/dobs.dat',dobs)
dlmwrite(['./data/dobs_snr_' num2str(snr) '.dat'],dobs_n)

save('./data/AcquisitionParameters.mat','s','G','nt','dt','t','janT','janX','nx','nz','dx','dz')

%% plots
close all
figure('units','normalized','outerposition',[0 0 1 1])
subplot(231); 
imagesc(x./1000,t,dobs,[-.01 .01]); title('Noiseless data')
xlabel('Distance (km)'); ylabel('Time (sec)')

subplot(232); 
imagesc(x./1000,t,dobs_n,[-.01 .01]); title(['Noisy data - SNR = ' num2str(snr)])
xlabel('Distance (km)'); ylabel('Time (sec)')

subplot(233); 
imagesc(x./1000,t,dobs-dobs_n,[-.01 .01]); title('Noise')
xlabel('Distance (km)'); ylabel('Time (sec)')

subplot(234); 
histogram(dobs-dobs_n,'normalization','probability'); title('Noise distribution')
xlabel('Noise'); ylabel('Probability distribution')

colormap bone

subplot(235); 
plot(dobs(:),dobs_n(:),'ks'); 
xlabel('Noiseless data'); 
ylabel('Noisy data')
y=linspace(min([dobs(:);dobs_n(:)]),max([dobs(:);dobs_n(:)]),length(dobs(:)));
hold on; plot(y,y,'r-','linewidth',2)
axis([min([dobs(:);dobs_n(:)]) max([dobs(:);dobs_n(:)]) min([dobs(:);dobs_n(:)]) max([dobs(:);dobs_n(:)])])
title(['Pearson Correlation: R = ' num2str(corr(dobs(:),dobs_n(:)))])

% Compute spectrum
L=nt; Fs=1/dt;

Ys = fft(s);
Y = fft(dobs);
Yn = fft(dobs_n);

P2s = abs(Ys/L);
P2 = abs(Y/L);
P2n = abs(Yn/L);

P1s = P2s(1:L/2+1);
P1 = P2(1:L/2+1);
P1n = P2n(1:L/2+1);

P1s(2:end-1) = 2*P1s(2:end-1);
P1(2:end-1) = 2*P1(2:end-1);
P1n(2:end-1) = 2*P1n(2:end-1);

f = Fs*(0:(L/2))/L;

subplot(236);
plot(f,P1,'linewidth',3) 
hold on; plot(f,P1s,'linewidth',3) 
hold on; plot(f,P1n,'linewidth',3) 
title('Single-Sided Amplitude Spectrum')
xlabel('f (Hz)')
ylabel('|P1(f)|')
axis([0 5*freqPICO -inf inf])
legend('Noiseless data','Seismic source','Noisy data')

saveas(gcf,'Observed_data','png')

%%
% well-positions in metters
p1=1500; p2=2000; p3=2500; p4=3000; p5=3500;
p1=p1/dx+1; p2=p2/dx+1; p3=p3/dx+1; p4=p4/dx+1; p5=p5/dx+1;


figure('units','normalized','outerposition',[0 0 1 1])
subplot(151); 
plot(dobs(:,p1),t,'linewidth',3); title(['Trace at ' num2str((p1-1)*dx/1000) ' km'])
hold on; plot(dobs_n(:,p1),t,'linewidth',1.5); 
xlabel('Distance (km)'); ylabel('Time (sec)')
ax = gca; ax.YDir = 'reverse';
legend('Noiseless','Noisy')


subplot(152); 
plot(dobs(:,p2),t,'linewidth',3); title(['Trace at ' num2str((p2-1)*dx/1000) ' km'])
hold on; plot(dobs_n(:,p2),t,'linewidth',1.5); 
xlabel('Distance (km)'); ylabel('Time (sec)')
ax = gca; ax.YDir = 'reverse';
legend('Noiseless','Noisy')


subplot(153); 
plot(dobs(:,p3),t,'linewidth',3); title(['Trace at ' num2str((p3-1)*dx/1000) ' km'])
hold on; plot(dobs_n(:,p3),t,'linewidth',1.5); 
xlabel('Distance (km)'); ylabel('Time (sec)')
ax = gca; ax.YDir = 'reverse';
legend('Noiseless','Noisy')


subplot(154); 
plot(dobs(:,p4),t,'linewidth',3); title(['Trace at ' num2str((p4-1)*dx/1000) ' km'])
hold on; plot(dobs_n(:,p4),t,'linewidth',1.5); 
xlabel('Distance (km)'); ylabel('Time (sec)')
ax = gca; ax.YDir = 'reverse';
legend('Noiseless','Noisy')


subplot(155); 
plot(dobs(:,p5),t,'linewidth',3); title(['Trace at ' num2str((p5-1)*dx/1000) ' km'])
hold on; plot(dobs_n(:,p5),t,'linewidth',1.5); 
xlabel('Distance (km)'); ylabel('Time (sec)')
ax = gca; ax.YDir = 'reverse';
legend('Noiseless','Noisy')


saveas(gcf,'Observed_data_traces','png')