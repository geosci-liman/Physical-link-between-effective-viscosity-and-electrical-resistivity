function nu = cal_visco_LabData(nx,nz,T,P,C0,rho)
% cal_visco_LabData -- A function for calculating effective viscosity from
%                       resistivity
% nx -- The number of horizontal grids in the profile (unit: km)
% nz -- The number of vetical grids in the profile (unit:km)
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
% C0 = 80*ones(size(C0));

ru    = 1.2;%%
ep    = 1e-15; 
Ad    = 1258925.412   ;
Aw = 794.3282;%%%
Ed = 510. ;
Ew = 470. ;%%
Vd = 14. ;
Vw = 24.;%%

sigma_hydrous = 10^(-1.37);%% Gardes et al., 2014
H_hydrous = 89;%%
alpha_G=1.79;%%

sigma_polaron = 10^2.34;%%
H_polaron = 144;%%

%%
nu = zeros(nz,nx);
nu0 = zeros(nz,nx);
rho0 = zeros(nz,nx);
for i = 1:nz
    for j = 1:nx
    if (T(i,j)< 800)
        nu = 10^25;
    elseif(T(i,j) >= 800 && T(i,j)<1300)
        nu0(i,j)=ep^((1.0-n)/n)*1e6*(Aw*C0(i,j)^ru*exp(-(1000*(Ew+P(i,j)*Vw)/R/T(i,j))))^(-1./n);   %%karato et al., 2003
        rho0(i,j) =1/(sigma_hydrous*C0(i,j)*exp(-1000*(H_hydrous-alpha_G*C0(i,j)^(1./3.))/R/T(i,j))) ; % Gardes et al., 2014
        nu(i,j) = nu0(i,j)*(rho(i,j)/rho0(i,j))^(5/3.5/3);
    elseif(T(i,j)>1300)
        nu0(i,j) = ep^((1.0-n)/n)*1e6*(Ad*exp(-(1000*(Ed+P(i,j)*Vd)/R/T(i,j))))^(-1./n); %%karato et al., 2003
        rho0(i,j) = 1.0/(sigma_polaron*exp(-1000.*H_polaron/R/T(i,j))) ;           %!Gardes et al., 2014
        nu(i,j) = nu0(i,j)*(rho(i,j)/rho0(i,j))^(2/3.5);
    else
        nu(i,j) = "NaN";
    end 
    end 
end

end