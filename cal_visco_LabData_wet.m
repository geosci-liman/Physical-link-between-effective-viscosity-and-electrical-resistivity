function nu = cal_visco_LabData_wet(nx,nz,T,P,C0,rho_ol)
R = 8.314472;
k = 1.380649*10^(-23) ; 

n  = 3.5;%%
% C0 = 80*ones(size(C0));%%

ru    = 1.2;%%
belta = 0.126 ;
ep    = 1e-15; epp=1e-15 ;%%
Ad    = 1258925.412   ;
Aw = 794.3282;%%%
Ed = 510. ;
Ew = 470. ;%%
Vd = 14. ;
Vw = 24.;%%

Ap = 10^2.6 ;
fo2 = 6.2e-7        ;
Qu2 = 449.0    ;
qs = 0.20     ;
np = 3.5    ;%%

delta_hp=0.91   ;
alpha = 0.09  ;
sigma_0 = 10^3.05 ;
re = 0.86;

sigma_p = 10^2.1   ;
qe = 0.05;
H_e = 154. ;

sigma_hydrous = 10^(-1.37);%% Gardes et al., 2014
H_hydrous = 89;%%
alpha_G=1.79;%%

sigma_polaron = 10^2.34;%%
H_polaron = 144;%%

%%


for i = 1:nz
    for j = 1:nx
    if (T(i,j)< 800)
        nu = 1;
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