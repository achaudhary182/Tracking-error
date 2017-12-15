% Tracking error
clear all; close all; clc;

dirname=('');
error_files='*.mat';

pixel_size=105.6;
DT=0.3;

variance_of_x_nm=[];
variance_of_y_nm=[];

track_error_file=dir(error_files);

for k=1:numel(track_error_file)
    a=load(track_error_file(k).name);
    
    tracking_erro_x_final=[];
    tracking_erro_y_final=[];
    
    
    for km=1%:numel(a.Molecule),
        %jp=1:numel(a.Molecule(km).Results(:,1));
        framek=a.Molecule(km).Results(:,1);   %Lists the frames in a trajectory
        timek=a.Molecule(km).Results(:,2);    %Lists the instantaneous times
        xk=a.Molecule(km).Results(:,3);   % Lists the x-position
        yk=a.Molecule(km).Results(:,4);   % Lists the y-position
        distk=a.Molecule(km).Results(:,5);    % Distance (calculated from x and y)
        fwhmk=a.Molecule(km).Results(:,6);    % Full Width Half Maximum
        intensityk=a.Molecule(km).Results(:,7);   % Intensity
        errork=a.Molecule(km).Results(:,8);   % Position error 
        
        %%%% X-position calculation 
        mu_x=mean(xk);  % mean of x-position
        x_upd=xk-mu_x;  % x-axis
        x_std=std(xk);  % standard deviation of x-position
        x_var=var(xk);   % variance in x-position
        
        %plot histogram
        figure(1)
        bins = -10:1:10;
        nx = hist(x_upd,bins);
        bar(bins, nx./sum(nx));
        title('x-position histogram');
        set(gca,'fontsize',18);
        
        % Gaussian fit
        p1_x = -.5 * ((xk - mu_x)/x_std) .^ 2;
        p2_x = (x_std * sqrt(2*pi));
        f_x = exp(p1_x) ./ p2_x; 
       
        hold on,
        figure(2)
        plot(x_upd,f_x,'ko','MarkerSize',4);
        variance_of_x_nm=[variance_of_x_nm,x_var];
        %title('tracking error x-position');
        tracking_error_x=2*x_std;
        tracking_erro_x_final=[tracking_erro_x_final,tracking_error_x];
        mean(tracking_erro_x_final)
        title(sprintf('Tracking error x (nm) =%f',tracking_error_x));
        %a=polyfit(x_upd,f_x,2);
        grid on
        set(gca,'fontsize',18);
        
        %%%% Y-position calculation (Gaussian fit)
        mu_y=mean(yk);  % mean of y-position
        y_upd=yk-mu_y;  % y-axis
        y_std=std(yk);  % standard deviation of y-position
        y_var=var(yk);   % variance in y-position
                
        %plot 
        figure(3)
        %bins = -20:2:20;
        ny = hist(y_upd,bins);
        bar(bins, ny./sum(ny));
        title('y-position histogram');
        set(gca,'fontsize',18);
        
        % Gaussian fit
        p1_y = -.5 * ((yk - mu_y)/y_std) .^ 2;  % -1/2.^((y-mu)./(sigma))
        p2_y = (y_std * sqrt(2*pi));    % (sigma*sqrt(2*pi))
        f_y = exp(p1_y) ./ p2_y; % exp(-1/2.^((y-mu)./(sigma)))./((sigma*sqrt(2*pi)))
        figure(4), 
        plot(y_upd,f_y,'.');
        variance_of_y_nm=[variance_of_y_nm,y_var];
        tracking_error_y=2*y_std;
        tracking_erro_y_final=[tracking_erro_y_final,tracking_error_y];
        mean(tracking_erro_y_final)
        title(sprintf('Tracking error y (nm) =%f',tracking_error_y));
        xlabel('nm');ylabel('Probability');
        % legend('x-position tracking error','y-position tracking error');
        % a=polyfit(y_upd,f_y,2);
        grid on
        set(gca,'fontsize',18);
        
        % Gaussian fit of x and y combined
        r=sqrt((xk.^2)+(yk.^2));
        mu_r=mean(r);
        r_upd=r-mu_r;
        r_std=std(r);
        p1_r = -.5 * ((r - mu_r)/r_std) .^ 2;  % -1/2.^((y-mu)./(sigma))
        p2_r = (r_std * sqrt(2*pi));    % (sigma*sqrt(2*pi))
        f_r = exp(p1_r) ./ p2_r; % exp(-1/2.^((y-mu)./(sigma)))./((sigma*sqrt(2*pi)))
        figure(7)
        plot(r_upd,f_r,'.');
        tracking_error_r=2*r_std;
        title(sprintf('Tracking error r (nm) =%f',tracking_error_r));
        xlabel('nm');ylabel('Probability');
        grid on
        set(gca,'fontsize',18);
        
        figure(6)
        plot(xk,yk);
        xlabel('x(nm)'),ylabel('y(nm)');set(gca,'fontsize',18);
    end
   
end