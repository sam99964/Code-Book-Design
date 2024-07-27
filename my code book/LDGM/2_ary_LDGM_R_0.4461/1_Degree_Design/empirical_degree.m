clear all;clc;


beta=1.1;
degree=[2];

for i=1:1:40,
    d_pre=degree(i);
    d_new=ceil(d_pre*beta);
    degree=[degree, d_new];
end