clear all

% 1.5 v -> 1.2 Hz
% 2.5 v -> 1.8 Hz
% 3.0 v -> 2.1 Hz


clf 
traj_expt=importdata('Position_3.0 v up 1 ball run 1 15.11.2019 exp 1821.dat');
traj_sim=importdata('Trajectory simulation 2.1 Hz rev1 rspa.dat');
%

subplot(1,3,1)
hold on
x=traj_expt.data(:,1);
y=traj_expt.data(:,2);
L=length(x);
time_duration_expt=L*(1/110); %% 110 frames per second (30 second)


[xc, yc, vxc,vyc,tq]= exp_velo(x,y);
dtq=tq(2)-tq(1);

duration_equal_to_expt=time_duration_expt/dtq;
t_eq=duration_equal_to_expt;
%p=plot(xc,yc,'r', 'LineWidth',5)

color_line(xc,yc,tq,'LineWidth',2)

colorbar
axis square
box on




%axis equal

a=13;


%  gamma_list=[5, 6.3, 6.9];
%     k2_list=[15, 13.2,   13]*a;
%     t2_list1=[.09, 0.04,  0.025];

gamma_list=[7.2];
k2_list=[0.54];
t2_list1=[.013];


[x_sim,y_sim,vx_sim, vy_sim, time_sim,gamma_1,a,k2,t2]= simulation(gamma_list(1), k2_list(1),t2_list1(1));

subplot(1,3,2)
hold on

[pks1, locs1] = findpeaks(xc);
[pks2, locs2] = findpeaks(x_sim);


mm=length(locs1);


%plot(x_sim(1:locs2(mm)),y_sim(1:locs2(mm)))
xsim_plot=x_sim(locs2(mm)+1:2*locs2(mm));
ysim_plot=y_sim(locs2(mm)+1:2*locs2(mm));
t_sim= time_sim(locs2(mm)+1:2*locs2(mm));
color_line(xsim_plot,ysim_plot,t_sim, 'LineWidth',2)

colorbar 
axis square
%axis equal

box on




subplot(1,3,3)
 hold on
 
% Expt 

% vxx=[vxc; vyc];
% [xh,yh]=velo_dist(vxx);
%pp= plot(xh(2:end)/dtq,yh(2:end),'r')

    f{1}= strcat( '1balls_',num2str(3.0,'%.1f'),'v_vx_vy_cameraframe histogram data (mm per sec) truncated.txt');
    data=load(f{1});
    xh=data(:,1);
    yh=data(:,2);
pp=plot(xh(2:end),yh(2:end),'r')
pp.Marker = 'o';
pp.MarkerSize = 15;
pp.MarkerFaceColor='r';
pp.MarkerEdgeColor='k';
pp.LineStyle='none';




 % sim
 [xh,yh]=velo_dist(vx_sim);
  p= plot(xh,yh*4,'k')
  
  p.LineWidth=2;
 
 

set(gca, 'YScale', 'log');
ylim([4e-4 2e-1])
legend_font_size1=18;
axis square
gamma_k_ratio=(gamma_list(1)/k2_list(1))*a;
LL=strcat('$\,\gamma=$',num2str(gamma_list(1),'%1.1f'),...
                          '$,\beta=$',num2str(1/(t2_list1(1)*gamma_list(1)),'%1.1f'), ...
                          '$\gamma$,\, $\gamma \tau_p/m_p$ =', num2str(gamma_k_ratio,'%1.1f'));
                 text(min(xlim), min(ylim)+1e-4 ,LL,'FontSize',legend_font_size1,'Interpreter','latex' );






update_legends()
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf,'Fig3.pdf','Resolution',600)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function[x,y,vx, vy, time,gamma_1,a,k2,t2]= simulation(gamma_list, kfactor,t2_list1)

%t2=1/beta; k=(m_p/tau_p)*a;

dt = 0.001;
t_final = 5000;
n_t = t_final/dt;
t = (0:dt:t_final-dt)';
gamma_1 =gamma_list;% 6.9;
a = 13; %5
%k2 = 50*a;

k2= (a/kfactor)*gamma_list;; %25*gamma_1; %13*a%(gamma_1/2);
w = 5;
fc=0;
ntraj = 1;
variance = zeros(10,1);
tstart = tic;
freq_compared=2.1;

t2_list=[t2_list1];




legend_font_size=22;
legend_font_size1=48;


for ii= 1:length(t2_list)
    
    cnt=1;
    for ff =freq_compared% 1:0.5:10
        
        T = 1/ff;
        %     k = 5*(1-1/f);
        tic
        omega = 2*pi*ff;
        
        if (ff<fc)
            t2 = 1*ff;
            gamma=0;
            k1=a;
            k=0;
            c=0;
            d=1;
        else
            t2 = t2_list(ii);
            gamma=gamma_1;
            k1=a;%a/(10);
            k=k2;
            c=1;
            d=1;
        end
        beta(ii)=1/t2;
        v = zeros(ntraj*n_t,1);
        vx_all = v;
        vy_all = v;
        vx_all_cam = v;
        vy_all_cam = v;
        v_cam = v;
        for j=1:ntraj
            vx = zeros(n_t,1);
            x = zeros(n_t,1);
            vx_cam = vx; vy_cam = vx;
            x_cam = x; y_cam = x;
            x(1) = 0;
            vx(1) = 0;
            vy = vx;
            y = x;
            y(1) = 0;
            vy(1) = 0;
            vp = vx;
            r1 = round(- t2*log(1-rand)/dt);             % random number from exponential distribution with mean t2
            r2 = round(- t2*log(1-rand)/dt);
            while r1==0
                r1 = round(- t2*log(1-rand)/dt);
            end
            %
            dnoise1 = ones(n_t,1);
            dnoise2 = dnoise1;
            %%
            randcount = 1;
            count = 0;
            for i=1:n_t
                if  randcount < r1
                    dnoise1(i) = 0;%0.2*rand(1)+0.8;
                    randcount = randcount+1;
                else
                    dnoise1(i) = 1;
                    count = count+1;
                    if count < w
                        randcount = r1;
                    else
                        r1 = round(- t2*log(1-rand)/dt);
                        if r1 < w
                            r1 = w;
                        end
                        %                     while r1==0
                        %                         r1 = round(- t2*log(1-rand)/dt);
                        %                     end
                        randcount = 1;
                        count = 0;
                    end
                end
            end
            %%
            randcount = 1;
            count = 0;
            for i=1:n_t
                if  randcount < r2
                    dnoise2(i) = 0;%0.2*rand(1)+0.8;
                    randcount = randcount+1;
                else
                    dnoise2(i) = 1;
                    count = count+1;
                    if count < w
                        randcount = r2;
                    else
                        r2 = round(- t2*log(1-rand)/dt);
                        randcount = 1;
                        count = 0;
                    end
                end
            end
            %%
            for i=2:n_t
                vpx(i) = -a*omega*sin(omega*i*dt);
                
                vpy(i) = a*omega*cos(omega*i*dt);
                
                flag_mix=0 ;
                
                if ff<fc
                    vx(i) = vx(i-1)- k1*omega^2*cos(omega*i*dt)*dt;
                    xp(i) =a*omega*sin(omega*i*dt);
                    vy(i) = vy(i-1)- k1*omega^2*sin(omega*i*dt)*dt;
                    xp(i) =a*omega*sin(omega*i*dt);
                    
                else
                    if flag_mix==1
                        q=.7;
                        
                        vx(i) = vx(i-1)+ (1-exp(-((ff-fc)/q)))*(- gamma*(vx(i-1)-vpx(i-1))*(1-dnoise1(i))*dt + 1*(-k*omega*sin(omega*i*dt))*dnoise1(i)*dt)+...
                            (exp(-(ff-fc)/q))*(- k1*omega^2*cos(omega*i*dt)*dt);
                        vy(i) = vy(i-1)+ (1-exp(-((ff-fc)/q)))*(- gamma*(vy(i-1)-vpy(i-1))*(1-dnoise2(i))*dt + 1*(k*omega*cos(omega*i*dt)))*dnoise2(i)*dt +...
                            (exp(-(ff-fc)/q))*(- k1*omega^2*sin(omega*i*dt)*dt);
                        
                    elseif flag_mix==0
                        vx(i) = vx(i-1)- gamma*(vx(i-1)-vpx(i-1))*(1-dnoise1(i))*dt + 1*(-k*omega*sin(omega*i*dt))*dnoise1(i)*dt ;
                        vy(i) = vy(i-1)- gamma*(vy(i-1)-vpy(i-1))*(1-dnoise2(i))*dt + 1*(k*omega*cos(omega*i*dt))*dnoise2(i)*dt   ;
                        
                    end
                    
                end
                x(i) = x(i-1)+vx(i-1)*dt;
                y(i) = y(i-1)+vy(i-1)*dt;
                
                
                
                time(i)=i*dt;
                
                
            end
            
            
            vx_all((j-1)*n_t+1:j*n_t) = vx;
            vy_all((j-1)*n_t+1:j*n_t) = vy;
            
            
            
            
        end
        
        
        
        
        
        
        
        
        
        
        
        toc
    end
end






end


function[xh, yh]= velo_dist(v)

figure

h=histogram(abs(v),100,'Normalization','probability');        % original 'BinWidth'=2
set(gca, 'YScale', 'log');
xh=linspace(h.BinLimits(1)+h.BinWidth/2,h.BinLimits(2)-h.BinWidth/2,h.NumBins)';
yh=h.BinCounts';
yh=yh/sum(yh);
close











end



function[xc, yc, vxc,vyc,tq]= exp_velo(x,y)


dt=1/110;
t=(1:length(x))'*dt;;
L =length(x);

tq = 0:dt/100:max(t);
xc = interp1(t,x,tq);
yc = interp1(t,y,tq);

xc=smoothdata(xc, 'gaussian',2);
yc=smoothdata(yc, 'gaussian',2);
vxc=diff(xc);
vyc=diff(xc);
end


function update_legends()

legend_font_size=18
legend_font_size1=24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,1);
hold on
ax = gca;
ax.FontSize = 14;

box on
title (' Experiment:  $f=2.1 \mathrm{Hz}$ ','FontSize',legend_font_size,'Interpreter','latex')
set(get(gca,'title'),'Position',[max(xlim)-20 max(ylim) 1.00011])


text(min(xlim), 115,'(a)','FontSize',legend_font_size1,'Interpreter','latex' )
xlabel('$x$(mm)','FontSize',legend_font_size1,'Interpreter','latex')
ylabel('$y$(mm)','FontSize',legend_font_size1,'Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2);
hold on
ax = gca;
ax.FontSize = 14;
box on


xlabel('$t\, \mathrm{(s)} \rightarrow $','FontSize',legend_font_size1,'Interpreter','latex')
ylabel('$v\, \mathrm{(mm/s)} $','FontSize',legend_font_size1,'Interpreter','latex')
text(min(xlim), max(ylim)+5,'(b)','FontSize',legend_font_size1,'Interpreter','latex' )
title (' Model:  $f=2.1 \mathrm{Hz}$ ','FontSize',legend_font_size,'Interpreter','latex')
set(get(gca,'title'),'Position',[max(xlim)-20 max(ylim) 1.00011])


xlabel('$x$(mm)','FontSize',legend_font_size1,'Interpreter','latex')
ylabel('$y$(mm)','FontSize',legend_font_size1,'Interpreter','latex')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,3);
legend
hold on
set(gca, 'YScale', 'log');

xlim([0, 150])
ylim([3e-4, .2])

ax = gca;
ax.FontSize = 14;
box on

ylabel('$P(|v_x|)$','FontSize',legend_font_size1,'Interpreter','latex')
xlabel('$|v_x| \mathrm{(mm/s)}   $','FontSize',legend_font_size1,'Interpreter','latex')
text(min(xlim), max(ylim)-0.06,'(c)','FontSize',legend_font_size1,'Interpreter','latex' )
legend('Expt.','Model','FontSize',legend_font_size,'Interpreter','latex')


legend boxoff




end


function h = color_line(x, y, c, varargin)
% color_line plots a 2-D "line" with c-data as color
%
%       h = color_line(x, y, c)
%       by default: 'LineStyle','-' and 'Marker','none'
%
%          or
%       h = color_line(x, y, c, mark) 
%          or
%       h = color_line(x, y, c, 'Property','value'...) 
%             with valid 'Property','value' pairs for a surface object
%
%  in:  x      x-data
%       y      y-data
%       c      3rd dimension for colouring
%       mark   for scatter plots with no connecting line
%
% out:  h   handle of the surface object
% (c) Pekka Kumpulainen 
%     www.tut.fi
h = surface(...
  'XData',[x(:) x(:)],...
  'YData',[y(:) y(:)],...
  'ZData',zeros(length(x(:)),2),...
  'CData',[c(:) c(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','none');
  
if nargin ==4
    switch varargin{1}
        case {'+' 'o' '*' '.' 'x' 'square' 'diamond' 'v' '^' '>' '<' 'pentagram' 'p' 'hexagram' 'h'}
            set(h,'LineStyle','none','Marker',varargin{1})
        otherwise
            error(['Invalid marker: ' varargin{1}])
    end
elseif nargin > 4
    set(h,varargin{:})
end
end

