function nu = cal_visco_LabData_wet(nx,nz,T,P,C0,rho_ol)
% cal_visco_LabData -- A function for calculating effective viscosity from
%                       resistivity under "wet" condition
% nx -- The number of horizontal grids in the profile (unit: km)
% nz -- The number of vertical grids in the profile (unit:km)
% T --  Temperature (unit:K)
% P -- Pressure (unit:GPa)
% C0 -- Water content (unit:ppm)
% rho -- resistivity

% by liman 2023-3-2 try to get the viscosity from resistivity
% Comments, bug reports and questions, please sent to:
% manli@zju.edu.cn.
% $Revision: 1.0

R = 8.314472;

n  = 3.5;%%
% C0 = 80*ones(size(C0));%%

ru    = 1.2;%%

ep    = 1e-15; 
Aw = 794.3282;%%%
Ew = 470. ;%%
Vw = 24.;%%

sigma_hydrous = 10^(-1.37);%% Gardes et al., 2014
H_hydrous = 89;%%
alpha_G=1.79;%%

%%


for i = 1:nz
    for j = 1:nx
    if (T(i,j)< 800)
        nu = 10^25;
    elseif(T(i,j) >= 800 )
        nu0(i,j)=ep^((1.0-n)/n)*1e6*(Aw*C0(i,j)^ru*exp(-(1000*(Ew+P(i,j)*Vw)/R/T(i,j))))^(-1./n);   %%karato et al., 2003
        rho0(i,j) =1/(sigma_hydrous*C0(i,j)*exp(-1000*(H_hydrous-alpha_G*C0(i,j)^(1./3.))/R/T(i,j))) ; % Gardes et al., 2014
        nu(i,j) = nu0(i,j)*(rho_ol(i,j)/rho0(i,j))^(5/3.5/3);
    else
        nu(i,j) = "NaN";
    end 
    end 
end

end
