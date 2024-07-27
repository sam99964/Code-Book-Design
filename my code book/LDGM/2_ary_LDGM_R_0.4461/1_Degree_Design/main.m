clear all;clc;

beta=1.1;
degree=[2];

for i=1:1:25,
    d_pre=degree(i);
    d_new=ceil(d_pre*beta);
    degree=[degree, d_new];
end

R=0.4461;
db=10;
finess=1e-4;

Aeq(1,:)=ones(1,numel(degree)); 
beq(1,:)=1;
% --> summation edge fraction = 1

d=degree;    % use above empirical result
Aeq(2,:)=ones(1,numel(degree))./d;
beq(2,:)=1/R/db;
% --> node number with code rate constraint

lb=zeros(numel(degree), 1);
ub= ones(numel(degree), 1);

x=[finess:finess:1-finess];

for i=1:1:numel(degree),
    power=degree(i);
    x_temp=ones(1,numel(x));
    for j=1:1:power-1,
        x_temp=x_temp.*x;
    end
    x_d_m1(i,:)=x_temp;
end


for i=1:1:numel(degree),
    power=degree(i);
    x_temp=ones(1,numel(x));
    for j=1:1:power-2,
        x_temp=x_temp.*x;
    end
    x_d_m2(i,:)=x_temp;
end

d_m1=d-1;

f=@(edge) max( (edge')*x_d_m1+ (db-1) * (1-x) .*( (d_m1.*edge') * x_d_m2 ) );
x0=zeros(numel(d), 1);

% options = optimset('Algorithm','interior-point');
% options = optimset('Algorithm','active-set');
options = optimset('Algorithm','sqp');
[edge_frac, f_val, exit_flag]=fmincon(f, x0, [], [], Aeq, beq, lb, ub,[],options);

disp('Ic_thr');
disp(1/f_val);

node_frac=edge_frac./(degree')*R*db;

n=116906;
node_num=n*node_frac;
node_num = node_num';
table = zeros(3,26);
table(1,:) = degree(1,:);
table(2,:) = node_num(:);
table(3,:) = floor(node_num(:));
a = sum(table(3,:))
total = table(1,:)*table(2,:)'
s = [0,0,0,0,0,1,2,1,0,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,1]
final_node = s+table(3,:);
total = sum(final_node)
t = degree * final_node'
rate = t/10/total
% Paper's data comparison:



















