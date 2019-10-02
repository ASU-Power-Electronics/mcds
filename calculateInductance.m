%% calculateInductance
% Calculates inductance of set of windings on the same magnetic core.  Uses
% Hurley's formulation (1994) of summation of images.  Accurate to ~10%* prior
% to inclusion of gamma terms; now accurate to arbitrary tolerance.
%
% *results below are pre-inclusion of gamma terms
% airTest = [25, 25, sqrt(15e-3*25e-3), sqrt(15e-3*25e-3), 250e-3, getGMD(10e-3, 10e-3)];
% gives 3.7011e-6, 28.1099e-6 (8e-6, 17e-6)
% coreTest = [25, 25, 15e-3, 25e-3, 15e-3, 25e-3, 10e-3, 10e-3, 10e-3, 10e-3, 250e-3, 10.8e-3, 75];
% gives 86.3393e-6, 40.4391e-6 (85e-6, 34e-6)
% for a total of 158.5894e-6 (144e-6)
% which is a 10.13% overestimate

function [M, Ml] = calculateInductance(Windings, Core, N_w, f_s, tol)
    global MU_0
    
    % core variable extraction
    length = Core.l_e;
    area = Core.A_e;
    mu_r = Core.mu_r;
    sigma_c = 1/Core.material.main.rho;
    b = sqrt(area/pi); % core radius, circular cross-section transformation
    
    omega = 2*pi*f_s;
    gamma_0 = sqrt(1j*omega*mu_r*MU_0*sigma_c);
    gammaTerm = real(besseli(1, gamma_0*b)/(gamma_0*b*besseli(0, gamma_0*b)));
    
    % Inductance matrices, self/mutual, and leakage
    M = zeros(N_w);
    Ml = M;
    
    for i = 1:N_w
        % variable extraction
        N1 = Windings{i}.N; % number of turns
        r1 = Windings{i}.diameter/2; % mean radius of coil
        rIn1 = r1 - (Windings{i}.N_L*Windings{i}.d_o)/2; % inner radius
        rOut1 = r1 + (Windings{i}.N_L*Windings{i}.d_o)/2; % outer radius
        h1 = rOut1 - rIn1; % radial height
        w1 = Windings{i}.tpl*Windings{i}.d_o*Windings{i}.NPW*Windings{i}.bifilar; % width (breadth) of coil
        A = sqrt(rIn1*rOut1); % geometric mean radius
        
        if N1 > 1
            GMDself = getGMD(h1, w1); % GMD of coil to itself
        else
            GMDself = 0;  % single filament has no GMD
        end
        
        [Lab, Lal] = airInductance(N1, N1, A, A, length, GMDself, tol);
        [Lcb, Lcl] = coreInductance(N1, N1, rIn1, rOut1, rIn1, rOut1, h1, w1, h1, w1, length, b, mu_r, sigma_c, omega, tol);
        
        M(i, i) = Lab + Lal + Lcb*(gammaTerm - 1/mu_r) + Lcl;
        Ml(i, i) = Lal + Lcl;
        
%         fprintf('\nCoil %d\n', i)
%         fprintf('Inductance without core = %g\n', Lab + Lal)
%         fprintf('Constant air term = %g\n', Lab)
%         fprintf('Leakage air term = %g\n', Lal)
%         fprintf('Inductance due to core = %g\n', Lcb + Lcl)
%         fprintf('Constant core term = %g\n', Lcb)
%         fprintf('Leakage core term = %g\n', Lcl)
%         fprintf('Self-inductance = %g\n', M(i, i))
%         fprintf('Self-leakage = %g\n', Ml(i, i))
        
        if i < N_w
            for j = i + 1:N_w
                % variable extraction
                N2 = Windings{j}.N;
                r2 = Windings{j}.diameter/2;
                rIn2 = r2 - (Windings{j}.N_L*Windings{j}.d_o)/2;
                rOut2 = r2 + (Windings{j}.N_L*Windings{j}.d_o)/2;
                h2 = rOut2 - rIn2;
                R = sqrt(rIn2*rOut2);
                w2 = Windings{j}.tpl*Windings{j}.d_o*Windings{j}.NPW*Windings{i}.bifilar;
                GMD = R - A;
                
                [Mab, Mal] = airInductance(N1, N2, A, R, length, GMD, tol);
                [Mcb, Mcl] = coreInductance(N1, N2, rIn1, rOut1, rIn2, rOut2, h1, w1, h2, w2, length, b, mu_r, sigma_c, omega, tol);
        
                M(i, j) = Mab + Mal + Mcb*(gammaTerm - 1/mu_r) + Mcl;
                M(j, i) = M(i, j);
                Ml(i, j) = Mal + Mcl;
                Ml(j, i) = Ml(i, j);
                
%                 fprintf('\nCoils %d and %d\n', i, j)
%                 fprintf('Inductance without core = %g\n', Mab + Mal)
%                 fprintf('Constant air term = %g\n', Mab)
%                 fprintf('Leakage air term = %g\n', Mal)
%                 fprintf('Inductance due to core = %g\n', Mcb + Mcl)
%                 fprintf('Constant core term = %g\n', Mcb)
%                 fprintf('Leakage core term = %g\n', Mcl)
%                 fprintf('Mutual inductance = %g\n', M(i, j))
%                 fprintf('Mutual leakage = %g\n', Ml(i, j))
                
            end
        end
    end
end