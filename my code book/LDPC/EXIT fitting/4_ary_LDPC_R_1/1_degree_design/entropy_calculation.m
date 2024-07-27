clear all;clc;

step_size=1e-5; % relate to the finess
m=4;            % interval size

interval=[-m/2, m/2];
disp('interval I is'); disp(interval);

point_set=[-m/2,-m/4, 0, m/4];
disp('point set is'); disp(point_set);

interval_1=[point_set(1):step_size:point_set(2)-step_size];    % @ region [-m/2,-m/4]
interval_2=[point_set(2):step_size:point_set(3)-step_size];    % @ region [-m/4, 0]
interval_3=[point_set(3):step_size:point_set(4)-step_size];    % @ region [0, m/4]
interval_4=[point_set(4):step_size:interval(2)];               % @ region [m/4, m/2]

%%% Calculate h(y|c0, c1)

    t=2.18;
    
    pdf_1=@(z) exp(-t.*z.^2)./(exp(-t.*(z-point_set(1)).^2)  +exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z+m-point_set(4)).^2));   % @ region [-m/2,-m/4]
    pdf_2=@(z) exp(-t.*z.^2)./(exp(-t.*(z-point_set(1)).^2)  +exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2));     % @ region [-m/4, 0]
    pdf_3=@(z) exp(-t.*z.^2)./(exp(-t.*(z-m-point_set(1)).^2)+exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2));     % @ region [0, m/4]
    pdf_4=@(z) exp(-t.*z.^2)./(exp(-t.*(z-m-point_set(1)).^2)+exp(-t.*(z-m-point_set(2)).^2)+exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2));     % @ region [m/4, m/2]
    
    norm=sum(pdf_1(interval_1).*step_size)+sum(pdf_2(interval_2).*step_size)+sum(pdf_3(interval_3).*step_size)+sum(pdf_4(interval_4).*step_size);
    disp('now norm is');    disp(norm);
    
    entropy_1=-sum(pdf_1(interval_1).*log2(pdf_1(interval_1)).*step_size);
    entropy_2=-sum(pdf_2(interval_2).*log2(pdf_2(interval_2)).*step_size);
    entropy_3=-sum(pdf_3(interval_3).*log2(pdf_3(interval_3)).*step_size);
    entropy_4=-sum(pdf_4(interval_4).*log2(pdf_4(interval_4)).*step_size);
    
    entropy=entropy_1+entropy_2+entropy_3+entropy_4;
    var=sum((interval_1.^2) .* pdf_1(interval_1) .* step_size)+sum( (interval_2.^2) .* pdf_2(interval_2) .* step_size)+...
        sum((interval_3.^2) .* pdf_3(interval_3) .* step_size)+sum( (interval_4.^2) .* pdf_4(interval_4) .* step_size);

% Get rate R from above simulation
R=log2(m)-entropy;
Gau_var=2.^(2.*entropy)./(2*pi*exp(1));    % equivalent Gaussian variance
loss_dB=10*log10(var./Gau_var);

% h(z)
disp('h(y|c0, c1) =');
disp(entropy);

%%% Calculate h(y|c0)

% @ region [-m/2,-m/4]
pdf_1=@(z) 0.5*exp(-t.*z.^2)./(exp(-t.*(z-point_set(1)).^2)  +exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z+m-point_set(4)).^2))+...
           0.5*exp(-t.*(z-point_set(2)).^2)./(exp(-t.*(z-point_set(1)).^2)  +exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z+m-point_set(4)).^2));
% @ region [-m/4, 0]
pdf_2=@(z) 0.5*exp(-t.*z.^2)./(exp(-t.*(z-point_set(1)).^2)  +exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2))+...
           0.5*exp(-t.*(z-point_set(2)).^2)./(exp(-t.*(z-point_set(1)).^2)  +exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2));
% @ region [0, m/4]
pdf_3=@(z) 0.5*exp(-t.*z.^2)./(exp(-t.*(z-m-point_set(1)).^2)+exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2))+...
           0.5*exp(-t.*(z-point_set(2)).^2)./(exp(-t.*(z-m-point_set(1)).^2)+exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2));
% @ region [m/4, m/2]
pdf_4=@(z) 0.5*exp(-t.*z.^2)./(exp(-t.*(z-m-point_set(1)).^2)+exp(-t.*(z-m-point_set(2)).^2)+exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2))+...
           0.5*exp(-t.*(z-m-point_set(2)).^2)./(exp(-t.*(z-m-point_set(1)).^2)+exp(-t.*(z-m-point_set(2)).^2)+exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2));

norm=sum(pdf_1(interval_1).*step_size)+sum(pdf_2(interval_2).*step_size)+sum(pdf_3(interval_3).*step_size)+sum(pdf_4(interval_4).*step_size);
    disp('now norm is');    disp(norm);
    
    entropy_1=-sum(pdf_1(interval_1).*log2(pdf_1(interval_1)).*step_size);
    entropy_2=-sum(pdf_2(interval_2).*log2(pdf_2(interval_2)).*step_size);
    entropy_3=-sum(pdf_3(interval_3).*log2(pdf_3(interval_3)).*step_size);
    entropy_4=-sum(pdf_4(interval_4).*log2(pdf_4(interval_4)).*step_size);
    
    entropy_new=entropy_1+entropy_2+entropy_3+entropy_4;
    
% h(z)
disp('h(y|c1) =');
disp(entropy_new);    


I_uc_i=2-entropy_new;  % I(c0_j;y_j) = I(c1_j;y_j) = H(y)-H(y|c_0) = H(y)-H(y|c_1) = log_2(m)-H(y|c_1)
I_uc_ni=entropy_new-entropy; % I(c0_j;y_j|c1_j) = I(c1_j;y_j|c0_j) = H(y|c_1)-H(y|c_1,c_0)
disp('I_uc_i=');
disp(I_uc_i);
disp('I_uc_ni=');
disp(I_uc_ni);
disp('sum rate =');
disp(I_uc_i+I_uc_ni);