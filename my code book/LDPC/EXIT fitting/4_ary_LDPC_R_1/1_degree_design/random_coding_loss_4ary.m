clear all;clc;

step_size=1e-5; % relate to the finess

m=4;    % interval size
interval=[-m/2, m/2];
disp('interval I is'); disp(interval);

point_set=[-m/2,-m/4, 0, m/4];
disp('point set is'); disp(point_set);

interval_1=[point_set(1):step_size:point_set(2)-step_size];    % @ region [-m/2,-m/4]
interval_2=[point_set(2):step_size:point_set(3)-step_size];    % @ region [-m/4, 0]
interval_3=[point_set(3):step_size:point_set(4)-step_size];    % @ region [0, m/4]
interval_4=[point_set(4):step_size:interval(2)];               % @ region [m/4, m/2]

t_test=[2.1:0.01:2.2];%best is near 2.005
entropy=zeros(1, numel(t_test));

for i=1:1:numel(t_test),
    t=t_test(i);
    
    pdf_1=@(z) exp(-t.*z.^2)./(exp(-t.*(z-point_set(1)).^2)  +exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z+m-point_set(4)).^2));   % @ region [-m/2,-m/4]
    pdf_2=@(z) exp(-t.*z.^2)./(exp(-t.*(z-point_set(1)).^2)  +exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2));     % @ region [-m/4, 0]
    pdf_3=@(z) exp(-t.*z.^2)./(exp(-t.*(z-m-point_set(1)).^2)+exp(-t.*(z-point_set(2)).^2)  +exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2));   % @ region [0, m/4]
    pdf_4=@(z) exp(-t.*z.^2)./(exp(-t.*(z-m-point_set(1)).^2)+exp(-t.*(z-m-point_set(2)).^2)+exp(-t.*(z-point_set(3)).^2)+exp(-t.*(z-point_set(4)).^2)); % @ region [m/4, m/2]
    
    norm=sum(pdf_1(interval_1).*step_size)+sum(pdf_2(interval_2).*step_size)+sum(pdf_3(interval_3).*step_size)+sum(pdf_4(interval_4).*step_size);
    disp('now norm is');    disp(norm);
    
    entropy_1=-sum(pdf_1(interval_1).*log2(pdf_1(interval_1)).*step_size);
    entropy_2=-sum(pdf_2(interval_2).*log2(pdf_2(interval_2)).*step_size);
    entropy_3=-sum(pdf_3(interval_3).*log2(pdf_3(interval_3)).*step_size);
    entropy_4=-sum(pdf_4(interval_4).*log2(pdf_4(interval_4)).*step_size);
    
    entropy(i)=entropy_1+entropy_2+entropy_3+entropy_4;
    var(i)=sum((interval_1.^2) .* pdf_1(interval_1) .* step_size)+sum( (interval_2.^2) .* pdf_2(interval_2) .* step_size)+...
           sum((interval_3.^2) .* pdf_3(interval_3) .* step_size)+sum( (interval_4.^2) .* pdf_4(interval_4) .* step_size);
    
end

% Get rate R from above simulation
R=log2(m)-entropy;
Gau_var=2.^(2.*entropy)./(2*pi*exp(1));    % equivalent Gaussian variance

loss_dB=10*log10(var./Gau_var);
plot(R, loss_dB);

%save four_ary_result;