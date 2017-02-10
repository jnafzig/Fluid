%% This script takes the output of the solver script and generates plot at
% regular time intervals.  Before running, either load solar.mat (for a
% simulation of the solar wind) or run solver.m (for the sod shock tube)
%
% units are arbitrary and don't correspond to a specific physical system.


%% Calculate some numbers

% center of mass
R = sum(diag(x)*densityRec)./sum(densityRec);

% speed of sound
cs = abs((gamma*pressureRec./densityRec).^.5);

% mach number
mach = velocityRec./cs;

% Energy components
kineticRec = 1/2*(velocityRec.*momentumRec);
internalRec = 1/(gamma-1)*(pressureRec);
potentialRec = densityRec.*(phi(1:end-1,ones(1,NumRecTimeSteps))...
                            +phi(2:end,ones(1,NumRecTimeSteps)))/2;

% eliminate zero momentum cases for log plotting
% momentumRec(momentumRec==0) = NaN;
% kineticRec(kineticRec==0) = NaN;

% frame(s) for setting ylim
FrameNumber = 1:NumRecTimeSteps;
% FrameNumber = 1:200;

% step through time
tstep = 1;

% plot just plots 1 through 3
STANDARDPLOT = false;


%# figure
set(gcf, 'Color','white')
Z = peaks; surf(Z);  axis tight
set(gcf, 'nextplot','replacechildren', 'Visible','off');


% set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');
pos(3) = 640; % Select the width of the figure in [cm]
pos(4) = 480; % Select the height of the figure in [cm]
set(gcf, 'Position', pos);

%# preallocate
nFrames = numel(FrameNumber);
% nFrames = 200;
mov(1:nFrames) = struct('cdata',[], 'colormap',[]);


% plot analytic solution for shock tube (STANDARDPLOT only)
RIEMANN = false;

% choose plots
plot1 = (densityRec);
plot1_min = min(min(plot1(:,FrameNumber)))-.1;
plot1_max = max(max(plot1(:,FrameNumber)))+.1;

plot2 = velocityRec;
plot2_min = min(min([plot2(:,FrameNumber),cs]))-.1;
plot2_max = max(max([plot2(:,FrameNumber),cs]))+.1;

plot3 = (pressureRec);
plot3_min = min(min(plot3(:,FrameNumber)))-.1;
plot3_max = max(max(plot3(:,FrameNumber)))+.1;

plot4 = (momentumRec);
plot4_min = min(min(plot4(:,FrameNumber)))-.01;
plot4_max = max(max(plot4(:,FrameNumber)))+.01;

plot5 = (kineticRec);
plot5_min = min(min(plot5(:,FrameNumber)))-.01;
plot5_max = max(max(plot5(:,FrameNumber)))+.01;

plot6 = (internalRec);
plot6_min = min(min(plot6(:,FrameNumber)))-.01;
plot6_max = max(max(plot6(:,FrameNumber)))+.01;


ColOrd = get(gca,'ColorOrder');
% NumRecTimeSteps
for i = 1:tstep:NumRecTimeSteps
    
    if STANDARDPLOT
        
        subplot(1,3,1);
        plot(x,plot1(:,i),'Color',ColOrd(1,:));
        title('density');
        xlim([a,b]);
        ylim([plot1_min,plot1_max]);
        axis square;
        line([a,b],[0,0],...
            'LineStyle','-.',...
            'Color',[.5,.5,.5]);
        
        subplot(1,3,2);
        plot(x,plot2(:,i),'Color',ColOrd(2,:));
        title('velocity');
        xlim([a,b]);
        ylim([plot2_min,plot2_max]);
        axis square;

        subplot(1,3,3);
        plot(x,plot3(:,i),'Color',ColOrd(3,:));
        title('pressure');
        xlim([a,b]);
        ylim([plot3_min,plot3_max]);
        axis square;
        
        if RIEMANN
            
            % calculate exact solution:
            riemann;
            
            I = .5-t(i).*(cL-velocityL);
            
            II = .5-t(i).*(cII-u);
            
            III = .5+t(i).*u;
            
            IV = .5+t(i).*s;
            
            subplot(1,3,1);  % density plot
            line([a,b],[0,0],...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

            line([a,b],[rhoII,rhoII],...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

            line([a,b],[rhoIII,rhoIII],...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

            line([I,I],ylim,...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);
            
            line([II,II],ylim,...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);
            
            line([III,III],ylim,...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

            line([.5+t(i).*(s),.5+t(i).*(s)],ylim,...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

            subplot(1,3,2); % velocity plot
            line([a,b],[0,0],...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

            line([a,b],[u,u],...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

            line([I,I],ylim,...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);
            
            line([IV,IV],ylim,...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

            subplot(1,3,3); % pressure plot
            line([a,b],[0,0],...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

            line([a,b],[p,p],...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

            line([I,I],ylim,...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);
            
            line([IV,IV],ylim,...
                'LineStyle','-.',...
                'Color',[.5,.5,.5]);

        end
        
        
    else
        
        subplot(2,3,1);
        plot(x,plot1(:,i),'Color',ColOrd(1,:));
        title('density');
        xlim([a,b]);
        ylim([plot1_min,plot1_max]);
        axis square;
        line([a,b],[0,0],...
            'LineStyle','-.',...
            'Color',[.5,.5,.5]);
        
%         line([R(i),R(i)],ylim,...
%             'LineStyle','-.',...
%             'Color',[.5,.5,.5]);

        subplot(2,3,2);
        plot(x,[plot2(:,i),cs(:,i)],'Color',ColOrd(2,:));
        title('velocity/sound speed');
        xlim([a,b]);
        ylim([plot2_min,plot2_max]);
        axis square;
        line([a,b],[0,0],...
            'LineStyle','-.',...
            'Color',[.5,.5,.5]);

        subplot(2,3,3);
        plot(x,plot3(:,i),'Color',ColOrd(3,:));
        title('pressure');
        xlim([a,b]);
        ylim([plot3_min,plot3_max]);
        axis square;
        line([a,b],[0,0],...
            'LineStyle','-.',...
            'Color',[.5,.5,.5]);

        subplot(2,3,4);
        plot(x,plot4(:,i),'Color',ColOrd(4,:));
        title('momentum');
        xlim([a,b]);
        ylim([plot4_min,plot4_max]);
        axis square;
        line([a,b],[0,0],...
            'LineStyle','-.',...
            'Color',[.5,.5,.5]);

        subplot(2,3,5);
        plot(x,plot5(:,i),'Color',ColOrd(5,:));
        title('kinetic energy');
        xlim([a,b]);
        ylim([plot5_min,plot5_max]);
        axis square;
        line([a,b],[0,0],...
            'LineStyle','-.',...
            'Color',[.5,.5,.5]);

        subplot(2,3,6);
        plot(x,plot6(:,i),'Color',ColOrd(6,:));
        title('internal energy');
        axis square;
        xlim([a,b]);
        ylim([plot6_min,plot6_max]);
        line([a,b],[0,0],...
            'LineStyle','-.',...
            'Color',[.5,.5,.5]);
    end
    
    fprintf('%5.0f%9.4f\n',i,t(i));
    if any(densityRec(:,i)<0);
        fprintf('negative\n');
    end

    mov(i) = getframe(gcf);
    pause(.05);
end

close(gcf)

%# save as AVI file
% movie2avi(mov, 'vid.avi', 'compression','None', 'fps',30);