clear;

% defining constants
hbar = 6.626e-34 / (2*pi); % plank's constant
lambda = 804.3e-9; % wavelength
kr = (2*pi)/lambda;
kx = linspace(-5*kr, 5*kr, 100);
m = 87*0.33*(1.66e-27)*7e4; % mass
Er = ((hbar^2)*(kr^2))/(2*m);
delta = 0;
delta1 = (-2*Er)/hbar;
delta2 = (-2*Er)/hbar;
omega = (4.85*Er)/hbar;
epsilon = 0;
theta1 = 0;
theta2 = 0;

%Defining empty lists and arrays that will be filled
E1 = zeros(3, length(kx));
E2 = zeros(3, length(kx));
Eb1 = zeros(3, length(kx));
Eb2 = zeros(3, length(kx));
Erf = zeros(3, length(kx));
Ebrf = zeros(3, length(kx));
Edrf = zeros(3, length(kx));
Ebdrf = zeros(3, length(kx));
pa = [];
pa2 = [];
pa3 = [];
pa1d = [];
pa2d = [];
pa3d = [];
pa_int = [];

%finding the Hamiltonian's eigenvalues
for ii = 1:length(kx)
    % H1 = Raman w/o Detuning, H2 = Raman w/Detuning
    % lowercase b indicates basis state
    H1 = [(hbar/(2*m))*(kx(ii)+2*kr)^2 - delta, omega/2, 0;
      omega/2, ((hbar/(2*m))*(kx(ii)^2)) - epsilon, omega/2;
      0, omega/2, (hbar/(2*m))*(kx(ii)-2*kr)^2 + delta];
    H2 = [(hbar/(2*m))*(kx(ii)+2*kr)^2 - delta2, omega/2, 0;
      omega/2, ((hbar/(2*m))*(kx(ii)^2)) - epsilon, omega/2;
      0, omega/2, (hbar/(2*m))*(kx(ii)-2*kr)^2 + delta2];
    E1(:,ii) = eig(H1);
    E2(:,ii) = eig(H2);
    Hb1 = [(hbar/(2*m))*(kx(ii)+2*kr)^2 - delta, 0, 0;
      0, ((hbar/(2*m))*(kx(ii)^2)) - epsilon, 0;
      0, 0, (hbar/(2*m))*(kx(ii)-2*kr)^2 + delta];
    Hb2 = [(hbar/(2*m))*(kx(ii)+2*kr)^2 - delta2, 0, 0;
      0, ((hbar/(2*m))*(kx(ii)^2)) - epsilon, 0;
      0, 0, (hbar/(2*m))*(kx(ii)-2*kr)^2 + delta2];
    Eb1(:,ii) = eig(Hb1);
    Eb2(:,ii) = eig(Hb2);
    % Hrf = rf w/o detuning, Hdrf = rf w/detuning
    Hrf = [((hbar/(2*m))*kx(ii)^2)-delta1, (1/2)*omega*(exp(1i*theta1)), 0;
        (1/2)*omega*(exp(-1i*theta1)), ((hbar/(2*m))*kx(ii)^2)-epsilon, (1/2)*omega*(exp(1i*theta2));
        0, (1/2)*omega*(exp(-1i*theta2)), ((hbar/(2*m))*kx(ii)^2)+delta2];
    Erf(:,ii) = eig(Hrf);
    Hbrf = [((hbar/(2*m))*kx(ii)^2)-delta1, 0, 0;
        0, ((hbar/(2*m))*kx(ii)^2)-epsilon, 0;
        0, 0, ((hbar/(2*m))*kx(ii)^2)+delta2];
    Ebrf(:,ii) = eig(Hbrf);
    Hdrf = [((hbar/(2*m))*kx(ii)^2)-0, (1/2)*omega*(exp(1i*theta1)), 0;
        (1/2)*omega*(exp(-1i*theta1)), ((hbar/(2*m))*kx(ii)^2)-epsilon, (1/2)*omega*(exp(1i*theta2));
        0, (1/2)*omega*(exp(-1i*theta2)), ((hbar/(2*m))*kx(ii)^2)+0];
    Edrf(:,ii) = eig(Hdrf);
    Hbdrf = [((hbar/(2*m))*kx(ii)^2)-0, 0, 0;
        0, ((hbar/(2*m))*kx(ii)^2)-epsilon, 0;
        0, 0, ((hbar/(2*m))*kx(ii)^2)+0];
    Ebdrf(:,ii) = eig(Hbdrf);
end

%finding the minimum points on the lowest band and checking for doubles
%(due to symmetry)
checker1 = 10;
checker2 = 10;
checkerrf = 0;
checkerdrf = 0;
for ii = 1:length(kx)
    if checker1 > E1(1, ii)
        checker1 = E1(1, ii);
        loc1 = ii;
    end
    if checker1 == E1(1, ii)
        loc1b = ii;
    end
    if checker2 > E2(1, ii)
        checker2 = E2(1, ii);
        loc2 = ii;
    end
    if checker2 == E2(1, ii)
        loc2b = ii;
    end
    if checkerrf > Erf(1, ii)
        checkerrf = Erf(1, ii);
        locrf = ii;
    end
    if checkerrf == Erf(1, ii)
        locrfb = ii;
    end
     if checkerdrf > Edrf(1, ii)
        checkerdrf = Edrf(1, ii);
        locdrf = ii;
    end
    if checkerdrf == Edrf(1, ii)
        locdrfb = ii;
    end
end

% calculating the PA rate and defining variables that are only used here
kx_temp = kr;
theta1 = linspace(0, 2*pi, 30);
delta1 = 0;
delta2_detuned = delta1;
omega1 = omega; %if omega is set close to 100 it matches Felicia's output
omega2 = omega1;
for ii = 1:30
    Hrf2 = [((hbar/(2*m))*kx_temp^2)-delta1, (1/2)*omega1*(exp(1i*theta1(ii))), 0;
        (1/2)*omega1*(exp(-1i*theta1(ii))), ((hbar/(2*m))*kx_temp^2)-epsilon, (1/2)*omega2*(exp(1i*theta2));
        0, (1/2)*omega2*(exp(-1i*theta2)), ((hbar/(2*m))*kx_temp^2)+delta2_detuned];
    [Eigvec0,b] = eig(Hrf2);
    pa_temp = (conj(Eigvec0(2,1)^2))*((Eigvec0(2,1)^2)) + 4*((conj(Eigvec0(1,1)*Eigvec0(3,1)))*(Eigvec0(1,1)*Eigvec0(3,1))) - 4*real((Eigvec0(2,1))^2 * conj(Eigvec0(1,1)) * conj(Eigvec0(3,1)));
    pa_int_temp = (conj(Eigvec0(2,1)^2)*(Eigvec0(2,1)^2)) + 4*((conj(Eigvec0(1,1)*Eigvec0(3,1)))*(Eigvec0(1,1)*Eigvec0(3,1)));
    pa = cat(2, pa, pa_temp);
    pa_int = cat(2, pa_int, pa_int_temp);
end


%plotting the energy-momentum dispersion relation and pa rate outputs
figure;
plot(kx/kr, E1, 'LineWidth', 2.5); hold on;
plot(kx/kr, Eb1, "k--", 'LineWidth', 2.5);
plot(kx(loc1)/kr, checker1, 'k*', 'LineWidth', 2.5);
plot(kx(loc1b)/kr, checker1, 'k*', 'LineWidth', 2.5);
xlabel('$$\frac{k_x}{k_r}$$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Energy $$(E_r)$$', 'Interpreter', 'latex', 'FontSize', 18);
title('Energy-Momentum Dispersion Curves for $$\delta=0$$', 'Interpreter', 'latex', 'FontSize', 20);
ylim([-4 10]);
grid on;

figure;
plot(kx/kr, E2, 'LineWidth', 2.5); hold on;
plot(kx/kr, Eb2, "k--", 'LineWidth', 2.5);
plot(kx(loc2)/kr, checker2, 'k*', 'LineWidth', 2.5);
xlabel('$$\frac{k_x}{k_r}$$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Energy $$(E_r)$$', 'Interpreter', 'latex', 'FontSize', 18);
title('Energy-Momentum Dispersion Curves for $$\hbar * \delta =-2E_r$$', 'Interpreter', 'latex', 'FontSize', 20);
ylim([-4 10]);
grid on;

figure;
plot(kx/kr, Erf, 'LineWidth', 2.5); hold on;
plot(kx/kr, Ebrf, "k--", 'LineWidth', 2.5);
plot(kx(locrf)/kr, checkerrf, 'k*', 'LineWidth', 2.5);
plot(kx(locrfb)/kr, checkerrf, 'k*', 'LineWidth', 2.5);
xlabel('$$\frac{k_x}{k_r}$$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Energy $$(E_r)$$', 'Interpreter', 'latex', 'FontSize', 18);
title('Energy-Momentum Dispersion Curves for Detuned RF-Coupling Scheme', 'Interpreter', 'latex', 'FontSize', 18);
ylim([-5 10]);
grid on;

figure;
plot(kx/kr, Edrf, 'LineWidth', 2.5); hold on;
plot(kx/kr, Ebdrf, "k--", 'LineWidth', 2.5);
plot(kx(locdrf)/kr, checkerdrf, 'k*', 'LineWidth', 2.5);
plot(kx(locdrfb)/kr, checkerdrf, 'k*', 'LineWidth', 2.5);
xlabel('$$\frac{k_x}{k_r}$$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Energy $$(E_r)$$', 'Interpreter', 'latex', 'FontSize', 18);
title('Energy-Momentum Dispersion Curves for RF-Coupling Scheme', 'Interpreter', 'latex', 'FontSize', 18);
ylim([-5 20]);
grid on;

figure;
plot(theta1 - theta2, pa, 'LineWidth', 2.5); hold on;
plot(theta1 - theta2, pa_int, 'LineWidth', 2.5);
xlabel('$$\theta_1 - \theta_2$$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$$\frac{k_{sup}}{k_{0,0}}$$', 'Interpreter', 'latex', 'FontSize', 18);
title('RF Normalized PA Rate Dependence on $$\theta_1 - \theta_2$$', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interference','No Interference', 'Interpreter', 'latex')
xticks(0:(pi/4):2*pi);
xlim([0 2*pi]);
grid on;

% making all of my variables symbolic so the eigenstates and PA rate can be
% solved analytically
kx = 0;
delta1 = 0;
delta2 = 0;
syms hbar
syms m 
syms omega1
syms omega2
syms theta1
syms epsilon
syms theta2
%analytically defining the matrix
Hrf_analytic = sym([0, (1/2)*omega1*(exp(1i*theta1)), 0;
        (1/2)*omega1*(exp(-1i*theta1)), 0, (1/2)*omega2*(exp(1i*theta2));
        0, (1/2)*omega2*(exp(-1i*theta2)), 0]);
[Erf_analytic, a] = eig(Hrf_analytic);
disp('omega 1 != omega 2');
disp('Eigenvalues:');
disp(a);
disp('bare state representation:');
bsr = Erf_analytic(1,1) * a(1,1) + Erf_analytic(2,1) * a(2,2) + Erf_analytic(3,1) * a(3,3);
disp(bsr);
disp('Eigenstates:');
disp(Erf_analytic);
pa_analytic = (conj(Erf_analytic(2,1))^2)*((Erf_analytic(2,1)^2)) + 4*((conj(Erf_analytic(1,1)*Erf_analytic(3,1)))*(Erf_analytic(1,1)*Erf_analytic(3,1))) - 4*real(Erf_analytic(2,1))^2 * Erf_analytic(1,1) * Erf_analytic(3,1);
disp('PA Rate:');
disp(pa_analytic);
pa_a = rewrite(pa_analytic, "sincos");
disp(pa_a);
% so it looks like the first and second expressions depend on theta 1 - theta 2
% because those two values are included in those expressions

%analytically defining the matrix for omega 1 = omega 2
syms omega
Hrf_analytic = sym([0, (1/2)*omega*(exp(1i*theta1)), 0;
        (1/2)*omega*(exp(-1i*theta1)), 0, (1/2)*omega*(exp(1i*theta2));
        0, (1/2)*omega*(exp(-1i*theta2)), 0]);
[Erf_analytic, a] = eig(Hrf_analytic);
disp('omega 1 = omega 2');
disp('Eigenvalues:');
disp(a);
disp('bare state representation:');
bsr = Erf_analytic(1,1) * a(1,1) + Erf_analytic(2,1) * a(2,2) + Erf_analytic(3,1) * a(3,3);
disp(bsr);
disp('Eigenstates:');
disp(Erf_analytic);
pa_analytic = (conj(Erf_analytic(2,1))^2)*((Erf_analytic(2,1)^2)) + 4*((conj(Erf_analytic(1,1)*Erf_analytic(3,1)))*(Erf_analytic(1,1)*Erf_analytic(3,1))) - 4*real(Erf_analytic(2,1))^2 * Erf_analytic(1,1) * Erf_analytic(3,1);
disp('PA Rate:');
disp(pa_analytic);
pa_a = rewrite(pa_analytic, "sincos");
disp(pa_a);