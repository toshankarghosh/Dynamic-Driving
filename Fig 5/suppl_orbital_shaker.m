%% 
clear all
clf

% 1.5 v -> 1.2 Hz
% 2.5 v -> 1.8 Hz
% 3.0 v -> 2.1 Hz

legend_font_size=18
legend_font_size1=24;
legend_font_size2=30;
scale_f=[1,1,1,1.1]

C=["(a)", "(b)", "(c)", "(d)"];
freq_drive=["Exeriment : f=1.2 Hz", "Exeriment : 1.8 Hz", "Exeriment : 2.1 Hz"];
freq=[1.2, 1.8, 2.1];
scale=[5.5 ,4.8 , 4];
drive=[1.5, 2.5 , 3.0];
num_part={"2900","220","40","1"};

color=[0 0.4470 0.7410;
   0.9290, 0.6940, 0.1250;
   0.4660, 0.6740, 0.1880;
   0.8500, 0.3250, 0.0980];

%figure


for j = 1:1



V=drive(j);
f{1} = strcat('2900balls_',num2str(V,'%.1f'),'v_vx_vy_cameraframe histogram data (mm per sec) truncated.txt');
f{2} = strcat( '220balls_',num2str(V,'%.1f'),'v_vx_vy_cameraframe histogram data (mm per sec) truncated.txt');
f{3} =strcat( '40balls_',num2str(V,'%.1f'),'v_vx_vy_cameraframe histogram data (mm per sec) truncated.txt');
f{4}= strcat( '1balls_',num2str(V,'%.1f'),'v_vx_vy_cameraframe histogram data (mm per sec) truncated.txt');

%
%subplot(2,2,j)
hold on


 for i= 1:4
     data=[];
     data=load(f{i});
     xh=data(:,1);
     yh=data(:,2)/scale_f(i);
     pp= plot(xh(1:end-1),yh(1:end-1));
     pp.Marker = 'o';
    pp.MarkerSize = 12;
    pp.MarkerFaceColor=color(i,:);
    pp.MarkerEdgeColor=color(i,:);
    pp.Color=color(i,:);
%     pp.LineStyle=':';
    pp.LineWidth=2;
    set(gca, 'YScale', 'log');
    ylim([4e-5 2e-1])
    %axis square
    ax = gca;
    ax.FontSize =18;
    box on  
    ylabel('$P(|v_x|)$','FontSize',legend_font_size1,'Interpreter','latex')
    xlabel('$|v_x| \mathrm{(mm/s)}   $','FontSize',legend_font_size1,'Interpreter','latex')
    

 end
    % title (freq_drive(j),'FontSize',legend_font_size2,'Interpreter','latex')
    % text(min(xlim), .3,C(j),'FontSize',legend_font_size2,'Interpreter','latex' )
     
    a=13; 
   
    
    
%     gamma_list=[5.2, 6.3, 7.2];
%     k2_list=[0.36, 0.48,   0.54];
%     t2_list1=[0.0750,  0.02, .013];
%     
    
%     gamma_list=[5, 6.3, 6.9];
%     k2_list=[15, 13.2,   15]*a;
%     t2_list1=[.09, 0.04,  0.018];
     
     
     if j<4
            if j==1
            
%              t2_list_n1   = [.07   0.042  0.035 ]; 
%              gamma_list_n1= [ 4.4   4.4   4.4  ];
%              k2_list_n1   = [ 15    15    15 ];
                
            for jj= 1: 1%length(t2_list_n1)    
                
%            [x_sim,y_sim,vx_sim, vy_sim, time_sim]=...
%                simulation(j,freq(j),gamma_list(jj),k2_list(jj),t2_list1(jj));
%             [xh,yh]=velo_dist(vx_sim);
% %             p= plot(xh,yh*scale(j),'Color',color(5-jj,:));
%                           p= plot(xh,yh*scale(j),'k');

             p.LineWidth=3;
            end
%             legend('Expt: N=2900','Expt: N=220','Expt: N=40',' Expt: N=1',...
%                    ' Model: N=1',' Model: N=40',' Model: N=220', ...
%                    'FontSize',legend_font_size,'Interpreter','latex','Location','southwest');
%            
%                gamma_k_ratio=(gamma_list(jj)/k2_list(jj))*a;
                
                
                      
                      
%               LL=strcat(' Model: N=1', ',\gamma=',num2str(gamma_list_n1(jj),'%1.1f'), '$,\beta=$',num2str(1/t2_list_n1(jj),'%1.1f'),...
%                          ',\gamma \tau/m =',num2str( gamma_k_ratio,'%1.1f'));
              legend('Expt: N=2900','Expt: N=220','Expt: N=40',' Expt: N=1',...
                   ' Model: N=1', ...
                   'FontSize',legend_font_size,'Interpreter','latex','Location','northeast');

             legend box off
%                          beta1=(1/(t2_list1(jj)))/gamma_list(jj);

%              LL=strcat('$\,\gamma=$',num2str(gamma_list(jj),'%1.1f'),...
%                           '$,\beta=$',num2str(beta1,'%1.1f'), ...
%                           '$\gamma$,\, $\gamma \tau/m$ =', num2str(gamma_k_ratio,'%1.1f'));
%                  text(min(xlim), min(ylim)+1e-4 ,LL,'FontSize',18,'Interpreter','latex' );

            else
                
                             
                [x_sim,y_sim,vx_sim, vy_sim, time_sim]=...
                simulation(j,freq(j),gamma_list(j),k2_list(j),t2_list1(j));
                [xh,yh]=velo_dist(vx_sim);
                %p= plot(xh,yh*scale(j),'Color',color(4,:));
                                          p= plot(xh,yh*scale(j),'k');

               
                p.LineWidth=3;
                legend('Expt: N=2900','Expt: N=220','Expt: N=40',' Expt: N=1',...
                   ' Model: N=1',...
                   'FontSize',legend_font_size,'Interpreter','latex','Location','northeast');
               
             legend box off
             gamma_k_ratio=(gamma_list(j)/k2_list(j))*a;
             beta1=(1/(t2_list1(j)))/gamma_list(j);
             
             LL=strcat('$\,\gamma=$',num2str(gamma_list(j),'%1.1f'),...
                          '$,\beta=$',num2str(beta1,'%1.1f'), ...
                          '$\gamma$,\, $\gamma \tau/m$ =', num2str(gamma_k_ratio,'%1.1f'));
                 text(min(xlim), min(ylim)+1e-4 ,LL,'FontSize',18,'Interpreter','latex' );
                
            end   
            
            
            
     end
     
     

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%subplot(2,2,4)
hold on
%  data=[];
%      data=load('1balls_1.5v_vx_vy_cameraframe histogram data (mm per sec) truncated.txt');
%      xh=data(:,1);
%      yh=data(:,2);
%      pp= scatter(xh,yh);
%      
     

%t2_list_n1   = [.07   0.042  0.035  0.025]; %
% t2_list_n1= 0.07:0.01:.14;
% 
% t2_list_n1=t2_list_n1/gamma_list(1);



% gamma_list=[5.2, 6.3, 7.2];
%     k2_list=[0.36, 0.48,   0.54];
%     t2_list1=[0.0750,  0.02, .013];
%     


%  gamma_list=[4.8, 6.3, 7.2];
%     k2_list=[0.3, 0.48,   0.54];
%   %  t2_list1=[0.09,  0.02, .013];
%     t2_list_n1=[0.075  0.044  0.039];    


gamma_list=[3.5, 6.3, 7.2];
    k2_list=[0.13, 0.48,   0.54];
  %  t2_list1=[0.09,  0.02, .013];
    t2_list_n1=[0.1  0.044  0.039];  
    
    
    
beta=(1./t2_list_n1)/(gamma_list(1));

    
for jj= 1: length(t2_list_n1)    
                
            [x_sim,y_sim,vx_sim, vy_sim, time_sim]=...
                simulation(j,freq(1),gamma_list(1),k2_list(1),t2_list_n1(jj));

            
            [xh,yh]=velo_dist(vx_sim);
             p= plot(xh,yh*scale(j),'k');%,'Color',color(5-jj,:));
             p.LineWidth=3;

end

set(gca, 'YScale', 'log');
    ylim([4e-5 5e-1])
    %axis square
    ax = gca;
    ax.FontSize =18;
    box on  
    ylabel('$P(|v_x|)$','FontSize',legend_font_size1,'Interpreter','latex')
    xlabel('$|v_x| \mathrm{(mm/s)}   $','FontSize',legend_font_size1,'Interpreter','latex')
        title ('Model','FontSize',legend_font_size2,'Interpreter','latex')

    for iN = 1:length(beta)
        legendCell{iN} = strcat('$\beta$ = ', num2str(beta(iN),'%.1f'),'$\gamma$');
    end
 legend(legendCell,'FontSize',18,'Interpreter','latex','Location','southwest','NumColumns',2)
 legend boxoff
%colororder(parula(length(beta)+1))
     text(min(xlim), 1,"(d)",'FontSize',legend_font_size2,'Interpreter','latex' )
   
    gamma_k_ratio=(gamma_list(1)/k2_list(1))*a;  
 LL=strcat('$\,\gamma=$',num2str(gamma_list(1),'%1.1f'),...
                          '$, \, \gamma \tau/m$ =', num2str(gamma_k_ratio,'%1.1f'), '$, \, f=1.2 \mathrm{Hz}$');
                 text(min(xlim), min(ylim)+.5e-2 ,LL,'FontSize',18,'Interpreter','latex' );
                 
                 
                 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% subplot(2,2,1)
% xlim([0 140])
% 
% subplot(2,2,2)
% xlim([0 200])
% 
% subplot(2,2,3)
% xlim([0 200])
% 
% 
% subplot(2,2,4)
% xlim([0 120])





% exportgraphics(gcf,'Fig5.pdf','Resolution',600)
%  exportgraphics(gcf,'Fig5.png','Resolution',600)
% 



             




function[x,y,vx, vy, time]= simulation(val,freq,gamma,kfactor,t2)

%clear all
%t2=1/beta; k=(m_p/tau_p)*a;
% freq
% gamma
% k2
% t2

dt = 0.001;
t_final = 5000;
n_t = t_final/dt;
t = (0:dt:t_final-dt)';
a = 13; %5
% gamma_list=[5, 6.3, 6.9];
% k2_list=[15, 13.2,   13]*a;
% t2_list1=[.09, 0.04,  0.025];




gamma_1 =gamma; %6.9(for 2.1 Hz);
%k2=k2%(gamma_1/2); (13*a for 2.1 Hz)
t2_list=[t2]; %[ .025] for 2.1 Hz



w = 5;
fc=0;
ntraj = 1;
variance = zeros(10,1);
tstart = tic;
freq_compared=freq





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
            k=gamma * (a/kfactor);
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
