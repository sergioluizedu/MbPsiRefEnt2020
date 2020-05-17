%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-Stack reflection data inversion based on Tsallis statistics   %%%%
% with on l-BFGS method                                              %%%%
%--------------------------------------------------------------------%%%%
% References:                                                        %%%%
%                                                                    %%%%
% de Lima, I.P. et. al, Tsallis Entropy, Likelihood, and the Robust  %%%%
% Seismic Inversion. Entropy 2020, 22, 464.                          %%%%
% https://doi.org/10.3390/e22040464                                  %%%%
%                                                                    %%%%
% da Silva, S. L. E. F. et. al, Robust full-waveform inversion using %%%% 
% q-statistics, Phys. A Stat. Mech. Appl., 2020, 548, 124473.        %%%%
% https://doi.org/10.1016/j.physa.2020.124473                        %%%%
%                                                                    %%%%
% Nocedal, J. et. al, Numerical Optimisation, 2nd ed.; Springer: New %%%%
% York, NY, USA, 2006.                                               %%%%
%                                                                    %%%%
% Byrd, R.H. et. al, A limited memory algorithm for bound constrained%%%%
% optimization. J. Sci. Comp. 1995, 16, 1190–1208.                   %%%%
% http://users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf   %%%%
%--------------------------------------------------------------------%%%%
% Last update: February, 2020                                        %%%%
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
% This code performs the Post-Stack inversion from reflectivity      %%%%
% series model, using l-BFGS method.
%--------------------------------------------------------------------%%%%

%% setup
clear all; clc; close all; addpath('./core/');

%% Inputs
iterMax=10;         % Maximum number of iterations
q=1;                % q-Parameter 
                    % q = 1         Convencional case; 
                    % 1 < q < 3     q-Gaussian case;
dir='./results/';   % Directory where the results will be saved
snr=30;             % snr=0 (Noiseless data)

%% Acquisition parameters
dtrue=dlmread('./data/dobs.dat'); 
if snr==0
    dobs=dtrue;
else
    dobs=dlmread(['./data/dobs_snr_' num2str(snr) '.dat']);
end
load('./data/AcquisitionParameters.mat');
% Validation of the q-value
if ((q<1) || (q>=3)); error('q-parameter invalid!'); end
%% Initial model
fileID = fopen(['./models/Initial_Marmousi2Ztime_sergioMod_dt_.001_dx_10_nt_2001_nx_1701_suav_250'],'r');
Zt0=fread(fileID,[2001 1701],'float32'); fclose(fileID);
Zt0=Zt0(janT,janX); Zt0=resizem(Zt0,[nt length(janX)]); 
ref0=.5.*[diff(log(Zt0)); zeros(1,nx)]; 
d0=G*ref0;      

x=0:dx:(nx-1)*dx; z=0:dz:(nz-1)*dz;

mkdir(dir)

%% misfit function
fh=@(m) misfit(m,dobs,G,nt,nx,q);


%% Inversion process
options.dir=dir; options.itermax=iterMax;
[ref_rec,phi] = mylbfgs(fh,ref0(:),options);
[ref_rec]=reshape(ref_rec,size(ref0));
dlmwrite([dir 'misfit.dat'],phi); dlmwrite([dir 'r_rec.dat'],ref_rec);


%% plot
close all
dmod=G*ref_rec;
figure('units','normalized','outerposition',[0 0 1 1])
subplot(241); 
imagesc(x./1000,t,d0,[-.01 .01]); title('Inital modelled data')
xlabel('Distance (km)'); ylabel('Time (sec)')

subplot(242); 
imagesc(x./1000,t,dmod,[-.01 .01]); title('Final Modelled data')
xlabel('Distance (km)'); ylabel('Time (sec)')

subplot(243); 
imagesc(x./1000,t,dobs,[-.01 .01]); title('Observed data')
xlabel('Distance (km)'); ylabel('Time (sec)')

subplot(244); 
imagesc(x,t,ref_rec,[-.01 .01]); title('Reflectivity: Inversion result')
xlabel('Distance (km)'); ylabel('Time (sec)')


subplot(245); 
imagesc(x./1000,t,dobs-dmod,[-.01 .01]); title('Data misfit')
xlabel('Distance (km)'); ylabel('Time (sec)')


colormap bone

subplot(246); 
plot(dobs(:),dmod(:),'ks'); 
xlabel('Observed data'); 
ylabel('Modelled dat')
y=linspace(min([dobs(:);dmod(:)]),max([dobs(:);dmod(:)]),length(dobs(:)));
hold on; plot(y,y,'r-','linewidth',2)
axis([min([dobs(:);dmod(:)]) max([dobs(:);dmod(:)]) min([dobs(:);dmod(:)]) max([dobs(:);dmod(:)])])
title(['Pearson Correlation: R = ' num2str(corr(dobs(:),dmod(:)))])

subplot(247); 
plot(phi./max(phi(:)),'ks-'); title('Convergence')
xlabel('Iteration'); ylabel('(Normalized) misfit function')
axis([0 length(phi)+1 -.08 1.08])

subplot(248); 
histogram(dtrue-dobs,'normalization','probability'); title('Hypothetical error distribution')
hold on;histogram(dtrue-dmod,'normalization','probability'); 
xlabel('Error'); ylabel('Probability distribution')
legend('True erro','Our error')

saveas(gcf,'Results','png')

%%
% well-positions in metters
p=[1500 2000 2500 3000 3500];
p=p./dx+1;

figure('units','normalized','outerposition',[0 0 1 1])
for ix=1:length(p)
    subplot(length(p),1,ix); 
    plot(t,dobs(:,p(ix)),'linewidth',3); title(['Trace at ' num2str((p(ix)-1)*dx/1000) ' km'])
    hold on; plot(t,dtrue(:,p(ix)),'linewidth',2); 
    hold on; plot(t,dmod(:,p(ix)),'linewidth',1.5); 
    ylabel('Displacement (km)'); xlabel('Time (sec)')
    ax = gca; ax.YDir = 'reverse';
    legend('Observed data','True data','Modelled data')
end

saveas(gcf,'Results_traces','png')