% % Convergence tests of 2D Schrodinger Eqn using ADI Scheme
% % Simply uncomment the upper or lower half of this script to run the
% % convergence test or actual solution error test, respectively. Note that
% % the initial ===== block must be run regardless to compute the ADI solutions.


% Test 1
idtype = 0;
vtype = 0;
idpar = [2, 3];
tmax = 0.05;
lambda = 0.05;
vpar = [0];
% for levels 6, 7, 8, 9;

% Uncomment block to rerun - dpsis are computed
% =========================================================
% [x, y, t, psi6, psire, psiim, psimod, v] = ...
%     sch_2d_adi(tmax, 6, lambda, idtype, idpar, vtype, vpar);
% [~, ~, ~, psi7, ~, ~, ~, ~] = ...
%     sch_2d_adi(tmax, 7, lambda, idtype, idpar, vtype, vpar);
% [~, ~, ~, psi8, ~, ~, ~, ~] = ...
%     sch_2d_adi(tmax, 8, lambda, idtype, idpar, vtype, vpar);
% [~, ~, ~, psi9, ~, ~, ~, ~] = ...
%     sch_2d_adi(tmax, 9, lambda, idtype, idpar, vtype, vpar);
% 
% % Downsample psi 7,8,9 to level 6
% psi7_ds = psi7(1:2:end, 1:2:end, 1:2:end);
% psi8_ds = psi8(1:4:end, 1:4:end, 1:4:end);
% psi9_ds = psi9(1:8:end, 1:8:end, 1:8:end);
% 
% % Take differences
% dpsi67 = psi7_ds - psi6;
% dpsi78 = psi8_ds - psi7_ds;
% dpsi89 = psi9_ds - psi8_ds;
% =========================================================

% % Manually compute l-2 norms, unsure if MATLAB norm() uses complex modulus
% norm67 = zeros(length(t), 1);
% norm78 = zeros(length(t), 1);
% norm89 = zeros(length(t), 1);
% for time = 1:length(t)
%     sum67 = 0;
%     sum78 = 0;
%     sum89 = 0;
%     for row =1:length(x)
%         for col = 1:length(y)
%             sum67 = sum67 + abs(dpsi67(time, row, col))^2;
%             sum78 = sum78 + abs(dpsi78(time, row, col))^2;
%             sum89 = sum89 + abs(dpsi89(time, row, col))^2;
%         end
%     end
%     norm67(time) = sqrt(sum67 / (length(x) * length(y)));
%     norm78(time) = sqrt(sum78 / (length(x) * length(y)));
%     norm89(time) = sqrt(sum89 / (length(x) * length(y)));
% end
% 
% % Plot scaled differences
% clf;
% figure(1)
% hold on;
% grid on;
% plot(t, norm67,'r-.', t, 4*norm78, 'b-.', t, 16*norm89, 'g-.');
% title("2D Schrodinger Scaled Convergence - Exact Family", 'interpreter', 'latex');
% xlabel("Time", 'interpreter', 'latex');
% ylabel("Solution Difference Norms $$\|\psi^l\|_2$$", 'interpreter', 'latex');
% legend("level 7-6", "level 8-7", "level 9-8", "location", "Northwest");
% axes('Position', [0.68, 0.16, 0.2, 0.2]);
% box on;
% plot(t(45:end), norm67(45:end),'r-.', t(45:end), ...
%     4*norm78(45:end), 'b-.', t(45:end), 16*norm89(45:end), 'g-.', "LineWidth", 1.2);
% hold off;

% ============================================================================

% % Compute exact solutions
% mx = idpar(1);
% my = idpar(2);
% 
% % level 6 exact
% exact6 = zeros(length(t), length(x), length(y));
% for nt = 1:length(t)
%     for i = 1:length(x)
%         for j = 1:length(y)
%             exact6(nt, i, j) = exp(-1i*(mx^2+my^2)*pi^2.*t(nt)) ...
%                 .*sin(mx*pi.*x(i)).*sin(my*pi.*y(j));
%         end
%     end
% end
% 
% % level 7 exact
% nx7 = 2^7 + 1;
% ny7 = 2^7 + 1;
% dx7 = 2^(-7);
% dy7 = dx7;
% dt7 = lambda * dx7;
% nt7 = round(tmax/dt7) + 1;
% 
% x7 = linspace(0.0, 1.0, nx7);
% y7 = linspace(0.0, 1.0, ny7);
% t7 = (0:nt7-1) * dt7;
% 
% exact7 = zeros(length(t7), length(x7), length(y7));
% for nt = 1:length(t7)
%     for i = 1:length(x7)
%         for j = 1:length(y7)
%             exact7(nt, i, j) = exp(-1i*(mx^2+my^2)*pi^2.*t7(nt)) ...
%                 .*sin(mx*pi.*x7(i)).*sin(my*pi.*y7(j));
%         end
%     end
% end
% 
% % level 8 exact
% nx8 = 2^8 + 1;
% ny8 = 2^8 + 1;
% dx8 = 2^(-8);
% dy8 = dx8;
% dt8 = lambda * dx8;
% nt8 = round(tmax/dt8) + 1;
% 
% x8 = linspace(0.0, 1.0, nx8);
% y8 = linspace(0.0, 1.0, ny8);
% t8 = (0:nt8-1) * dt8;
% 
% exact8 = zeros(length(t8), length(x8), length(y8));
% for nt = 1:length(t8)
%     for i = 1:length(x8)
%         for j = 1:length(y8)
%             exact8(nt, i, j) = exp(-1i*(mx^2+my^2)*pi^2.*t8(nt)) ...
%                 .*sin(mx*pi.*x8(i)).*sin(my*pi.*y8(j));
%         end
%     end
% end
% 
% % level 9 exact
% nx9 = 2^9 + 1;
% ny9 = 2^9 + 1;
% dx9 = 2^(-9);
% dy9 = dx9;
% dt9 = lambda * dx9;
% nt9 = round(tmax/dt9) + 1;
% 
% x9 = linspace(0.0, 1.0, nx9);
% y9 = linspace(0.0, 1.0, ny9);
% t9 = (0:nt9-1) * dt9;
% 
% exact9 = zeros(length(t9), length(x9), length(y9));
% for nt = 1:length(t9)
%     for i = 1:length(x9)
%         for j = 1:length(y9)
%             exact9(nt, i, j) = exp(-1i*(mx^2+my^2)*pi^2.*t9(nt)) ...
%                 .*sin(mx*pi.*x9(i)).*sin(my*pi.*y9(j));
%         end
%     end
% end
% 
% % take differences between exact and computed solutions
% dpsi6e = exact6 - psi6;
% dpsi7e = exact7 - psi7;
% dpsi8e = exact8 - psi8;
% dpsi9e = exact9 - psi9;
% 
% % Downsample exact differences
% dpsi7e = dpsi7e(1:2:end, 1:2:end, 1:2:end);
% dpsi8e = dpsi8e(1:4:end, 1:4:end, 1:4:end);
% dpsi9e = dpsi9e(1:8:end, 1:8:end, 1:8:end);
% 
% % Manually compute l-2 norms, unsure if MATLAB norm() uses complex modulus
% norm6e = zeros(length(t), 1);
% norm7e = zeros(length(t), 1);
% norm8e = zeros(length(t), 1);
% norm9e = zeros(length(t), 1);
% for time = 1:length(t)
%     sum6 = 0;
%     sum7 = 0;
%     sum8 = 0;
%     sum9 = 0;
%     for row =1:length(x)
%         for col = 1:length(y)
%             sum6 = sum6 + abs(dpsi6e(time, row, col))^2;
%             sum7 = sum7 + abs(dpsi7e(time, row, col))^2;
%             sum8 = sum8 + abs(dpsi8e(time, row, col))^2;
%             sum9 = sum9 + abs(dpsi9e(time, row, col))^2;
%         end
%     end
%     norm6e(time) = sqrt(sum6 / (length(x) * length(y)));
%     norm7e(time) = sqrt(sum7 / (length(x) * length(y)));
%     norm8e(time) = sqrt(sum8 / (length(x) * length(y)));
%     norm9e(time) = sqrt(sum9 / (length(x) * length(y)));
% end

% % Plot scaled errors
% clf;
% hold on;
% grid on;
% plot(t, norm6e,'r-s', t, 4*norm7e, 'b-s', t, 16*norm8e, 'g-s', t, 64*norm9e, 'm-s');
% title("2D Solution Errors - Exact Family", 'interpreter', 'latex');
% xlabel("Time", 'interpreter', 'latex');
% ylabel("Exact Solution Difference Norms $$\|E(\psi^l)\|_2$$", 'interpreter', 'latex');
% legend("level 6", "level 7", "level 8", "level 9", "location", "Northwest");
% hold off;

