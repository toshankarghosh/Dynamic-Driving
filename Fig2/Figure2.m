clear all
clf

%%parameters
dt = 0.001;
t_final = 5000;
n_t = t_final/dt;
t = (0:dt:t_final-dt)';
s_factor=[1 1 2 1];
a = 5;
freq_compared_list=[4.5  5.5  6.0    9.5 ]
f_name{1}=strcat('Expt_velo_dist_4p6Hz.dat');
f_name{2}=strcat('Expt_velo_dist_5p6Hz.dat');

f_name{3}=strcat('Expt_velo_dist_6Hz.dat');
f_name{4}=strcat('Expt_velo_dist_9p6Hz.dat');








t2_list=[ .01];
w = 1;
fc=4.4;


ntraj = 1;
variance = zeros(10,1);
tstart = tic;



colour_count=1;
  
    cnt=1;

    for ff = 1:0.1:10

        T = 1/ff;
        tic
        omega = 2*pi*ff;
        
        if (ff<fc)
            t2 = 1*ff;
            gamma=0;
            k1=a;
            k=0;
            c=0;
            d=1;
            k2=0;
        else
            
            fc=4.2;
            delta_f=(ff-fc)/fc; 
            n=1.5;
            gamma= 2.2+ 10.7/delta_f;
            k=gamma * (a/12);
            t2=  0.0042*delta_f^1+0.00064;
            k2=k;
        end
        
        
        GL(cnt)=gamma; % gamma
        KL(cnt)=k;
        T2L(cnt)=k;
        FreqL(cnt)=ff;
        beta(cnt)=1/t2;
        
        
        
        
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
                    dnoise2(i) = 0;
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
                
                
                if ff<fc
                    vx(i) = vx(i-1)- k1*omega^2*cos(omega*i*dt)*dt;
                    xp(i) =a*omega*sin(omega*i*dt);
                    vy(i) = vy(i-1)- k1*omega^2*sin(omega*i*dt)*dt;
                    xp(i) =a*omega*sin(omega*i*dt);
                    
                else
                    
                    
                        vx(i) = vx(i-1)- gamma*(vx(i-1)-vpx(i-1))*(1-dnoise1(i))*dt + 1*(-k*omega*sin(omega*i*dt))*dnoise1(i)*dt ;
                        vy(i) = vy(i-1)- gamma*(vy(i-1)-vpy(i-1))*(1-dnoise2(i))*dt + 1*(k*omega*cos(omega*i*dt))*dnoise2(i)*dt   ;
                        
                    
                    
                end
                x(i) = x(i-1)+vx(i-1)*dt;
                y(i) = y(i-1)+vy(i-1)*dt;
                
                
                
                time(i)=i*dt;
                
                
            end
            
            
            vx_all((j-1)*n_t+1:j*n_t) = vx;
            vy_all((j-1)*n_t+1:j*n_t) = vy;
            
            
            if any(freq_compared_list(:) == ff)==1;
                
                ind=find(freq_compared_list == ff);
                
                f_name{ind}
                velo_dist(vx_all,f_name{ind},ff);
                
            end
            
            
            if ff == 6
                vx_comp6=vx;
                vx_comp_all6=vx_all;
                vpx_comp6=vpx;
                vpy_comp6=vpy;
                x_comp6=x;
                y_comp6=y;
            end
            
        end
        
        
        
        flag=1; % calculate amp/phase from velo data else calculate  from  pos data
        if flag==1
            [ampl,phase_deg, var_amp, var_phase ] = amplitude_signal(vx,vpx, dt,ff);
            phase(cnt,2)=ampl/(2*pi*ff);
            phase(cnt,4)=var_amp/(2*pi*ff);
        else
            [ampl,phase_deg, var_amp, var_phase ] = amplitude_signal(x,xp, dt,ff);
            phase(cnt,2)=ampl;
            phase(cnt,4)=var_amp;
        end
        
        phase(cnt,1)=ff;
        %
        
        phase(cnt,3)=phase_deg;
        phase(cnt,5)=var_phase ;
        cnt=cnt+1;
        
        
        
        
        
        
        
        toc
    end
    %  end
    
    plot_response(phase, fc);
    
    
    
    
    plot_vel_traj(time, vx_comp6,vpx_comp6,k,a , 6);
    
    


update_legends(fc)
 plot_gamma(FreqL,GL,T2L)

 exportgraphics(gcf,'Fig2.png','Resolution',600)
 exportgraphics(gcf,'Fig2.pdf','Resolution',600)
% 




function velo_dist(v,f_name,freq)

%  v=vx;
data_expt=load(f_name);

figure
h=histogram(abs(v),100,'Normalization','probability');        % original 'BinWidth'=2
set(gca, 'YScale', 'log');
xh=linspace(h.BinLimits(1)+h.BinWidth/2,h.BinLimits(2)-h.BinWidth/2,h.NumBins)';
yh=h.BinCounts';
yh=yh/sum(yh);
close


subplot(3,3,7);
hold on

if freq==6
    s_f=2;
else
    s_f=1;
end
    

    p= scatter(data_expt(:,1), data_expt(:,2)/s_f,100, 'filled');
    pp= plot(xh,yh,'k');
    pp.LineWidth=3;








end





function [ampl,phase_deg, var_amp, var_phase ] = amplitude_signal(vc,vp, dt,ff)


L1=length(vc);
% vp=xp;
N=50; % split the time series in to N parts
for i = 1:N
    L=floor(L1/N);
    velo=vc((i-1)*L+1:i*L);
    FFT_v = fft(velo);
    P2 = abs(FFT_v/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    Fs=1/dt;
    freq = Fs*(0:(L/2))/L;
    % semilogx(freq,P1)
    % pause(1)
    [val,idx]=min(abs(freq-ff));
    ampl_1(i)=abs(P1(idx));
    
    
    a1=vp;
    b1=vc;
    %  calculate the phase difference
    phase_rad=acos(dot(a1,b1)/(norm(a1)*norm(b1)));
    phase_deg_1(i)=phase_rad*360/(2*pi);
    
end

ampl=mean(ampl_1);
phase_deg=mean(phase_deg_1);
var_amp=var(ampl_1);
var_phase=var(phase_deg_1);

end







function plot_response(phase,fc)


subplot(3,3,8);
hold on
data_amp=load('Expt_coin_disp_amp_6Hz.dat');

    p=errorbar(data_amp(1:end,1),data_amp(1:end,2),data_amp(1:end,3),'o');
    p.Marker = 'o';
    p.MarkerSize = 10;
    p.MarkerFaceColor='r';



pp=scatter(phase(:,1),phase(:,2),150,'filled')
pp.MarkerEdgeColor='k';
pp.MarkerFaceColor='k';
pp.MarkerFaceAlpha=0.25;








subplot(3,3,9);
hold on
data_phase=load('Expt_coin_disp_phase_6Hz.dat');

    p=errorbar(data_phase(1:end,1),data_phase(1:end,2)+3,data_phase(1:end,3),'o');
    p.Marker = 'o';
    p.MarkerSize = 10;
    p.MarkerFaceColor='r';



pp=scatter(phase(:,1),phase(:,3),150,'filled');
pp.MarkerEdgeColor='k';
pp.MarkerFaceColor='k';
pp.MarkerFaceAlpha=0.25;



% ylabel('$\phi $','FontSize',legend_font_size1,'Interpreter','latex')
% xlabel('$f(Hz) $','FontSize',legend_font_size1,'Interpreter','latex')
%

end



function plot_vel_traj(time, vx,vpx, k2,a,f)


%figure

    subplot(3,3,1:3);
    hold on
    data =load('Expt_velo_traj_6p0Hz.dat');
    t_coin=data(:,1);
    v_coin=data(:,2);
    
    t_plate=data(:,3);
    omega=(2*pi*6);
    v_plate=-(a)*omega*sin(omega*t_coin +0.56);%data(:,4);

    
    hold on
    
    p2=plot(t_plate,v_plate,'r');
    p2.LineWidth = 2;
   
    p1=plot(t_coin,v_coin,'k');
    p1.LineWidth = 2;
    
    




subplot(3,3,4:6);
hold on


    pp2= plot(time,vpx,'r');
    pp2.LineWidth = 2;




    pp1= plot(time,vx,'k');
    pp1.LineWidth = 2;



end


function update_legends(fc)

legend_font_size=24
legend_font_size1=24;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,1:3);
hold on
xlim([65 70])
ylim([-210, 210]);
ax = gca;
ax.FontSize = 20;

box on
title (' Experiment:  $f=6 \, \mathrm{Hz}$ ','FontSize',legend_font_size,'Interpreter','latex')
set(get(gca,'title'),'Position',[max(xlim)-.35 max(ylim) 1.00011])


text(min(xlim), 250,'(a)','FontSize',legend_font_size1,'Interpreter','latex' )
xlabel('$t\, \mathrm{(s)} \rightarrow $','FontSize',legend_font_size1,'Interpreter','latex')
ylabel('$v\, \mathrm{(mm/s)} $','FontSize',legend_font_size1,'Interpreter','latex')
legend('$v_p$','$v_c$','Big Marker','FontSize',legend_font_size,'Interpreter','latex')
legend boxon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,4:6);
hold on
xlim([65 70])
ylim([-210, 210]);
ax = gca;
ax.FontSize = 20;
box on


xlabel('$t\, \mathrm{(s)} \rightarrow $','FontSize',legend_font_size1,'Interpreter','latex')
ylabel('$v\, \mathrm{(mm/s)} $','FontSize',legend_font_size1,'Interpreter','latex')
text(min(xlim), 250,'(b)','FontSize',legend_font_size1,'Interpreter','latex' )
title (' Model:  $f=6\, \mathrm{Hz}$','FontSize',legend_font_size,'Interpreter','latex')
set(get(gca,'title'),'Position',[max(xlim)-.35 max(ylim) 1.00011])

legend('$v_p$','$v_c (M)$','FontSize',legend_font_size,'Interpreter','latex')
legend boxon


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,7);
legend
f=get(gca,'Children');

hold on

set(gca, 'YScale', 'log');

xlim([0, 180])
ylim([1e-5, .05])
%ylim([1e-5, .2])

ax = gca;
ax.FontSize = 20;
box on

ylabel('$P(|v|)$','FontSize',legend_font_size1,'Interpreter','latex')
xlabel('$|v| \mathrm{(mm/s)}   $','FontSize',legend_font_size1,'Interpreter','latex')
text(min(xlim), .1,'(c)','FontSize',legend_font_size1,'Interpreter','latex' )
legend([f(8), f(6), f(4), f(2), f(1)],'4.6 Hz','5.6 Hz', '6.0 Hz', '9.6 Hz','Model','FontSize',20,'Interpreter','latex','Location','southwest')


legend boxoff


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(3,3,8);
hold on
xline(fc,'--k','LineWidth',2);
xlim([1 10])
ylim([0, 5.5]);
ax = gca;
ax.FontSize = 20;

box on
ylabel('$x_{c0} (\mathrm{mm})$','FontSize',legend_font_size1,'Interpreter','latex')
xlabel('$f(Hz) $','FontSize',legend_font_size1,'Interpreter','latex')

text(min(xlim),6.0,'(d)','FontSize',legend_font_size1,'Interpreter','latex' )
legend('Expt.','Model','FontSize',legend_font_size,'Interpreter','latex','Location','southwest')
legend boxoff
ylabel('$x_{c0} (\mathrm{mm})$','FontSize',legend_font_size1,'Interpreter','latex')
xlabel('$f \mathrm{(Hz)}   $','FontSize',legend_font_size1,'Interpreter','latex')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,3,9);
hold on
xline(fc,'--k','LineWidth',2);
xlim([1 10])
ylim([0, 100]);
box on
ax = gca;
ax.FontSize = 20;



ylabel('$\phi\mathrm{(Deg.)} $','FontSize',legend_font_size1,'Interpreter','latex')
xlabel('$f\mathrm{(Hz)} $','FontSize',legend_font_size1,'Interpreter','latex')
text(min(xlim), 110,'(e)','FontSize',legend_font_size1,'Interpreter','latex' )
legend('Expt.','Model','FontSize',legend_font_size,'Interpreter','latex','Location','southeast')
legend boxoff


end


function plot_gamma(FreqL,GL,T2L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(3,3,8);
hold on
p = get(gca, 'Position');
h = axes('Parent', gcf, 'Position', [p(1)+0.13 p(2)+0.12 p(3)-0.135 p(4)-0.14]);
hold(h); 
scatter(h, FreqL,GL,50, 'filled')
box on
ax = gca;
ax.FontSize = 12;
ylabel('$\gamma$','FontSize',18,'Interpreter','latex')
xlabel('$f$','FontSize',18,'Interpreter','latex')





end

