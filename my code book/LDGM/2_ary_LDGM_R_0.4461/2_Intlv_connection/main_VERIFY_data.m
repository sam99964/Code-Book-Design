clear all;clc;

Dci=[1:1:11, 13:2:23, 24, 26, 27, 30, 33, 34, 37, 41, 42, 46, 49, 51, 57, 59, 63];
 
dc_node=[342, 54331, 18339, 7779, 3000, 2745, 2491, 2119, 917, 974, 946, 786, 653, ...
         577, 539, 518, 1, 475, 1, 415, 365, 338, 1, 343, 320, 1, 271, 1, 211, 136, 1, 64];

Dvi=[10];
dv_node=44610;

total_CND_edge=sum(Dci.*dc_node);
total_VND_edge=sum(Dvi.*dv_node);

disp('total CND nodes');
disp(sum(dc_node));

disp('total VND nodes');
disp(sum(dv_node));

disp('total CND edge');
disp(total_CND_edge);

disp('total VND edge');
disp(total_VND_edge);