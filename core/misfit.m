%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Misfit function based on q-Gaussian distribution                   %%%%
%--------------------------------------------------------------------%%%%
% References:                                                        %%%%
%                                                                    %%%%
% da Silva, S. L. E. F. et. al, Robust full-waveform inversion using %%%% 
% q-statistics, Phys. A Stat. Mech. Appl., 2020, 548, 124473.        %%%%
% https://doi.org/10.1016/j.physa.2020.124473                        %%%%
%                                                                    %%%%
% f (r,q) = 1/(q-1) sum{ ln[1+(q-1)/(3-q)(dobs-dmod(r))^2]  }        %%%% 
%--------------------------------------------------------------------%%%%
% Use:                                                               %%%%
%   [f,g] = misfit(m,dobs,G,nt,nx,q)                                 %%%%
%                                                                    %%%%
% Input:                                                             %%%%
%       m       - Reflectivy model                                   %%%%
%       dobs    - Observed/True data                                 %%%%
%       G       - Convolution operator                               %%%%
%       nt, nx  - Number of time and spatial samples                 %%%%
%       q       - q-Parameter                                        %%%%
% Output:                                                            %%%%
%       f       - Misfit function value                              %%%%
%       g       - Gradient of the same size as the input vector m    %%%%
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

function [f,g]=misfit(m,dobs,G,nt,nx,q)

m=reshape(m,[nt,nx]);

dmod=G*m;

if q==1
    f=.5*norm(dmod(:)-dobs(:),2)^2;
    g=G'*(dmod-dobs);
else
    f=sum((1/(q-1)).*log(1+((q-1)/(3-q)).*(dmod(:)-dobs(:)).^2));
    g=((1/(q-1)).*G)'*(2*((q-1)/(3-q)))*(dmod-dobs)./(1+((q-1)/(3-q)).*(dmod-dobs).^2);
end

g=g(:);
end