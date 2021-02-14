function [x, y, t, psi, psire, psiim, psimod, v] = ...
    sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)
%
% Inputs
%
%   tmax:   Maximum integration time
%   level:  Discretization level
%   lamdba: dt/dx
%   idtype: Selects initial data type
%   idpar:  Vector of initial data parameters
%   vtype:  Selects potential type
%   vpar:   Vector of potential parameters
%
% Outputs
%
%   x:      Vector of x coordinates [nx]
%   y:      Vector of y coordinates [ny]
%   t:      Vector of t coordinates [nt]
%   psi:    Array of computed psi values [nt x nx x ny]
%   psire:  Array of computed psi_re values [nt x nx x ny]
%   psiim:  Array of computed psi_im values [nt x nx x ny]
%   psimod: Array of computed sqrt(psi psi*) values [nt x nx x ny]
%   v:      Array of potential values [nx x ny]
% =======================================================================

nx = 2^level + 1;
ny = 2^level + 1;
dx = 2^(-level);
dy = dx;
dt = lambda * dx;
nt = round(tmax/dt) + 1;

% Initialize outputs
x = linspace(0.0, 1.0, nx);
y = linspace(0.0, 1.0, ny);
t = (0:nt-1) * dt;

psi    = zeros(nt, nx, ny);
psire  = zeros(nt, nx, ny);
psiim  = zeros(nt, nx, ny);
psimod = zeros(nt, nx, ny);

v = zeros(nx, ny);

% Select and set initial data
if idtype == 0
    % Exact family
    mx = idpar(1);
    my = idpar(2);
    for n_x = 1:nx
        for n_y = 1:ny
            psi(1, n_x, n_y) = sin(mx*pi*x(n_x)) * sin(my*pi*y(n_y));
        end
    end
elseif idtype == 1
    % Boosted Gaussian
    x0 = idpar(1);
    y0 = idpar(2);
    delx = idpar(3);
    dely = idpar(4);
    px = idpar(5);
    py = idpar(6);
    for n_x = 1:nx
        for n_y = 1:ny
            psi(1, n_x, n_y) = exp(1i*px*x(n_x)) * exp(1i*py*y(n_y)) ...
                * exp(-((x(n_x)-x0)^2/delx^2 + (y(n_y)-y0)^2/dely^2));
        end
    end
end

% Update other solutions to reflect initial data
psire(1, :, :)  = real(psi(1, :, :));
psiim(1, :, :)  = imag(psi(1, :, :));
psimod(1, :, :) = abs(psi(1, :, :));

% Select and set potential (only non-zero if vtype == 1
if vtype == 1
    % Rectangluar barrier or well
    xmin = vpar(1);
    xmax = vpar(2);
    ymin = vpar(3);
    ymax = vpar(4);
    Vc = vpar(5);
    v(xmin:xmax, ymin:ymax) = Vc;
    
elseif vtype == 2
    % Double slit
    x1 = vpar(1);
    x2 = vpar(2);
    x3 = vpar(3);
    x4 = vpar(4);
    Vc = vpar(5);
    
    j_prime = (ny-1)/4 + 1;
    
    % Set v to Vc for j = j_prime, j = j_prime+1
    v(:, j_prime)   = Vc;
    v(:, j_prime+1) = Vc;
    % Set slit openings to 0
    v(x1:x2, j_prime)   = 0;
    v(x3:x4, j_prime+1) = 0;
end


% Initialize x storage for sparse matrix and RHS
dl_x = zeros(nx, 1);
d_x  = zeros(nx, 1);
du_x = zeros(nx, 1);

% Set up tridiagonal system (LHS coefficents)
dl_x = -((1i*dt) / (2*dx^2)) .* ones(nx, 1);
d_x  =  (1 + (1i*dt)/dx^2) .* ones(nx, 1);
du_x = dl_x;

% Fix up x boundary cases
d_x(1) = 1.0;
du_x(2) = 0.0;
dl_x(nx-1) = 0.0;
d_x(nx) = 1.0;

% Define sparse matrix A_x
A_x = spdiags([dl_x d_x du_x], -1:1, nx, nx);

% Compute solution
for n = 1:nt - 1
    % Stage 1: For each j=2,3,..ny-1 solve for psi^{n+0.5}_{i,j}, i=1,2,..n
    % allocate intermediate storage
    halfpsi = zeros(nx, ny);
    for h = 2:ny-1
        f1a = zeros(nx, 1);
        f1b = zeros(nx, 1);
        yyoperator = (psi(n, 2:nx-1, h+1) - 2.*psi(n, 2:nx-1, h) + ... 
            psi(n, 2:nx-1, h-1))./dy^2;
        f1a(2:nx-1) = psi(n, 2:nx-1, h) + ((1i*dt)/2) .* yyoperator - ...
            ((1i*dt)/2) .* v(2:nx-1, h).'.*psi(n, 2:nx-1, h);
        
        f1b(2:nx-1) = f1a(2:nx-1) + ((1i*dt)/(2*dx^2)) .* ...
            (f1a(3:nx) - 2.*f1a(2:nx-1) + f1a(1:nx-2));
        
        % solve and update psi
        halfpsi(:, h) = A_x \ f1b;
        
    end
    
    % Stage 2: For each i=2,3,..nx-1 solve for psi^{n+1}_{i,j}, j=1,2,...n
    for u = 2:nx-1
        dl_y = zeros(nx, 1);
        d_y  = zeros(nx, 1);
        du_y = zeros(nx, 1);
        f_y  = zeros(nx, 1);

        % Create new tridiagonal system
        dl_y = (-0.5i*dt / dy^2) .* ones(nx, 1);
        d_y  = (1 + 1i*dt/dy^2 + 1i.*dt/2.*v(u, :).') .*ones(nx, 1);
        du_y = (-0.5i*dt / dy^2) .* ones(nx, 1);
        
        %BCs
        d_y(1) = 1.0;
        du_y(2) = 0.0;
        dl_y(nx-1) = 0.0;
        d_y(nx) = 1.0;
        
        A_y = spdiags([dl_y d_y du_y], -1:1, nx, nx);
        
        f_y(2:nx-1) = halfpsi(u, 2:nx-1);
        f_y(1)  = 0.0;
        f_y(nx) = 0.0;
        
        % Solve
        psi(n+1, u, :) = A_y \ f_y;
    end
    
    % Update remaining output vectors
    psire(n+1, :, :) = real(psi(n+1, :, :));
    psiim(n+1, :, :) = imag(psi(n+1, :, :));
    psimod(n+1, :, :) = abs(psi(n+1, :, :));
end
end

