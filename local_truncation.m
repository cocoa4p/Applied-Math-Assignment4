% CURRENT

% Local: the error for a single time step
% log-log fit to estimate order p

%[p,k]: the regressed values for relationship y = k*x^p
%function [p,k] = loglog_fit(x_regression,y_regression,varargin)

function local_truncation()
    tref = 0.492;
    sol = @solution01;
    rate = @rate_func01;

    % Single step function methods

    Embedded=@RK_step_embedded;

    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;...
                    1/5, 0, 0, 0,0,0,0;...
                    3/40, 9/40, 0, 0, 0, 0,0;...
                    44/45, -56/15, 32/9, 0, 0, 0,0;...
                    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
                    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
                    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

    h_step = logspace(-5, 1, 100); % 1e-5 to 1e1 100 points
    
    analytical_diff  = zeros(size(h_step));
    errs_Embedded = zeros(size(h_step));
    

    XA = sol(tref); % finding exact X at the tref X spot

    for r = 1:length(h_step)
        h = h_step(r);

        
        [XB_Embedded, ~] = Embedded(rate, tref, XA, h, DormandPrince);  % one step from exact state
        X_true = sol(tref + h);

        analytical_diff(r) = norm(X_true - XA);
        errs_Embedded(r) = norm(XB_Embedded - X_true);
       

        % errs(i) = norm(XB - X_true);        % norm handles vector/scalar

    end
    % [p,k] = loglog_fit(h, errs(i), 1);
    % plot(errs(i), h);
    
    figure(1);
%     loglog(h_step,analytical_diff,'ko','MarkerFaceColor','k','MarkerSize',1);
%     hold on
    loglog(h_step,errs_Embedded,'ro','MarkerFaceColor','r','MarkerSize',1); hold on;
   

    filter_params = struct();
    filter_params.max_xval = 1;
    filter_params.min_yval = 1e-10; 


%     [p1,k1] = loglog_fit(h_step, analytical_diff, filter_params);
      [p2,k2] = loglog_fit(h_step, errs_Embedded, filter_params)
 
%     loglog(h_step,k1*h_step.^p1,'g','LineWidth',1)
    hold on;
    loglog(h_step,k2*h_step.^p2,'r','LineWidth',1)
%     loglog(h_step,k3*h_step.^p3,'b','LineWidth',1)
% 
%     loglog(h_step,k4*h_step.^p4,'y','LineWidth',1)
%     loglog(h_step,k5*h_step.^p5,'m','LineWidth',1)


    %legend(sprintf('Forward Euler (p = %.4f)', p2));
    xlabel('h step')
    ylabel('Error')
    title('Local Truncation Error to Solution 1')

    
    disp(p2)
    xlabel('')
    ylabel('Error')

    % is y = k*x^p the same as e = O*h^p? 
    % [p,k] = loglog_fit(h_step, errs, 1)
 
end