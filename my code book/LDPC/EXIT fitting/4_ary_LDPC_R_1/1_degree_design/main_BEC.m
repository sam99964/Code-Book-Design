clear all;clc; % Ref: Brink_TCOM_04 

R=1;

Iuc_i =0.5128;
Iuc_ni=0.5189;

step=0.01;
delta=0.001;
n=10^5; % block length-->[total VND=2*n]

% ��Ӽ����� input
Ia_CND=[0:step:1];

% �M�w CND �� distribution
Dci=[8];
aci=[1];
bci=aci.*Dci./(sum(aci.*Dci));

% [ CND �� output ] = [ VND �� input ]
Ie_CND=zeros(1, numel(Ia_CND));
for i=1:1:numel(Dci),
    degree=Dci(i);
    tmp=bci(i)* CND_formula(Ia_CND,degree);
    Ie_CND=Ie_CND+tmp;
end

%figure(1);plot(Ie_CND, Ia_CND);

% [ VND �� output ] = [ CND �� input ]
Dvi=[2:1:20];   % �i�H��tune�A�U���Q�n�� linprog �ӸѰ��D�A�u�����ܼƬO bvi

% �y�z lb �M ub, lb <= x <= ub
    lb = zeros(numel(Dvi),1);
    ub =  ones(numel(Dvi),1);

% �y�z Aeq �M beq, Aeq*x=beq
    Aeq=zeros(2,numel(Dvi));
    beq=zeros(2, 1);
    % (1) sum(bvi)=1
    Aeq(1,:)=ones(1, numel(Dvi));
    beq(1) = 1;
    % (2) Total edge ������
    Aeq(2,:)=1./Dvi;
    beq(2,:)=(sum((1./Dci).*bci) /(2-R) )*2;

% �y�z A �M b, Ax<=b �A�n�� Icb_new >= �쥻�� Icb
% �]���o�̪����󦡬O <= �A�ҥH�n�����G���W�t��

    Ie_DET=(1-Ie_CND)*Iuc_i+Ie_CND*Iuc_ni;

    A=zeros( numel(Ia_CND), numel(Dvi) );
    for i=1:1:numel(Dvi),
        degree=Dvi(i);
        tmp=VND_formula( Ie_CND, Ie_DET, degree );
        A(:,i)=tmp';
    end
    A=-A;
    
    b=zeros( numel(Ia_CND), 1);
    b(end)=1;   % �̫�@�I��̳�=1
    b( 1:(end-1), 1)= Ia_CND( 1: (end-1))+delta;
    b=-b;
    
% �y�z f�A�n�� f'*x �V�p�V�n
% �Ʊ�EXIT chart�U�������n�V�p�V�n
    f=zeros( numel(Dvi), 1);
    for i=1:1:numel(Dvi),
        degree=Dvi(i);
        tmp=VND_formula( Ie_CND, Ie_DET, degree );
        interval=Ie_CND(2:end)-Ie_CND(1:(end-1));
        area=sum(tmp(1:end-1).*interval);
        f(i)=area;
    end
    
% �i�H apply linprog
x = linprog(f,A,b,Aeq,beq,lb,ub);
bvi=x';
% check
disp(sum(bvi));

% �e�X����
Ie_VND_new=zeros(1, numel(Ie_CND));
for i=1:1:numel(Dvi),
    degree=Dvi(i);
    temp=VND_formula( Ie_CND, Ie_DET, degree );
    Ie_VND_new=Ie_VND_new+bvi(i)*temp;
end

x_axis=Ie_CND;
y_axis=[Ia_CND; Ie_VND_new];
plot(x_axis, y_axis);

dv_node= (2*n)*(bvi./Dvi)./sum(bvi./Dvi);
dc_node= n*(2-R)*(bci./Dci)./sum(bci./Dci);
total_edge=sum(dc_node.*aci.*Dci);

lala=floor(dv_node);


