% % Uncomment relevant portion to run

% idtype = 0;
% vtype = 0;
% idpar = [2, 3];
% tmax = 0.05;
% lambda = 0.05;
% vpar = [0];
% 
% [x, y, t, psi, psire, psiim, psimod, v] = ...
%     sch_2d_adi(tmax, 7, lambda, idtype, idpar, vtype, vpar);

% Exact family, no potential
% clf;
% for time = 1:length(t)
%     surf(squeeze(psiim(time, :, :)));
%     xlabel("x");
%     ylabel("y");
%     zlabel("$$\mathrm{Re}(\psi^7)$$", 'interpreter', 'latex');
%     title("Exact family - Zero Potential $$m_x=2, ~ m_y=3$$", 'interpreter', 'latex');
%     txt = ["t = " + num2str(t(time)) + ' s'];
%     text(0, 120, 1, txt, 'Fontsize', 11);
%     drawnow
%     pause();
% end





% ========================================================================
% Boosted Gaussian with double slit
idtype = 1;
vtype = 2;
idpar = [0.5, 1, 0.2, 0.3, 0, 0];
tmax = 0.05;
lambda = 0.05;
vpar = [15, 50, 80, 125, 30000];

[x, y, t, psi, psire, psiim, psimod, v] = ...
    sch_2d_adi(tmax, 7, lambda, idtype, idpar, vtype, vpar);

% === Graphics ===
plotenable = 1;
pausesecs = 0.005;

avienable = 1;
avifilename = 'boosted_slit_im.avi';
aviframerate = 25;


if avienable
    aviobj = VideoWriter(avifilename);
    aviobj.Quality = 100;
    open(aviobj);
end

for time = 1:length(t)
    if plotenable
        clf;
        hold on;
        box on;
        surf(squeeze(psiim(time, :, :)));
        zlim([-0.8, 1.28]);
        view(-26, 21);
        xlabel("x");
        ylabel("y");
        zlabel("$$\mathrm{Im}(\psi^7)$$", 'interpreter', 'latex');
        drawnow;
        if avienable
            if time == 1
                framecount = 3 * aviframerate;
            else
                framecount = 1;
            end
            for iframe = 1 : framecount
                writeVideo(aviobj, getframe(gcf));
            end
        end

        pause(pausesecs);
    end
end
    

if avienable
   close(aviobj);
   fprintf('Created video file: %s\n', avifilename);
end


% ========================================================================
% % Boosted gaussian with barrier across x
% idtype = 1;
% vtype = 1;
% idpar = [10, 10, .5, .5, 1, 1];
% tmax = 0.05;
% lambda = 0.05;
% vpar = [1, 129, 40, 80, 10000];
% 
% [x, y, t, psi, psire, psiim, psimod, v] = ...
%     sch_2d_adi(tmax, 7, lambda, idtype, idpar, vtype, vpar);
% 
% % === Graphics ===
% plotenable = 1;
% pausesecs = 0.005;
% 
% avienable = 1;
% avifilename = 'boosted_barrier_mod.avi';
% aviframerate = 25;
% 
% 
% if avienable
%     aviobj = VideoWriter(avifilename);
%     aviobj.Quality = 100;
%     open(aviobj);
% end
% 
% for time = 1:length(t)
%     if plotenable
%         clf;
%         hold on;
%         box on;
%         surf(squeeze(psimod(time, :, :)));
%         zlim([-0.1*10^(-284), 10*10^(-284)]);
%         view(-26, 28);
%         xlabel("x");
%         ylabel("y");
%         zlabel("$$|\psi^7|$$", 'interpreter', 'latex');
%         drawnow;
%         if avienable
%             if time == 1
%                 framecount = 3 * aviframerate;
%             else
%                 framecount = 1;
%             end
%             for iframe = 1 : framecount
%                 writeVideo(aviobj, getframe(gcf));
%             end
%         end
% 
%         pause(pausesecs);
%     end
% end
%     
% 
% if avienable
%    close(aviobj);
%    fprintf('Created video file: %s\n', avifilename);
% end

