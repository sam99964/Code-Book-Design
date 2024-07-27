clear all;clc;


Dci=[2,3,4,5,6,7,8,9,10,11,13,15,17,19,24,27,30,33,57,63,70];

dc_node=[65292, 19361, 9011, 5793, 2386, 3095, 2228, 790, 1012, 866, 1468, 1546, 559, 29, 622, 923, 735, 2, 350, 800, 38];

Dvi=[10];
dv_node=52154;

% interleaver connection principle
% (1) no parallel edges
% (2) S-random requirement
% (2) no two CND have the same check equations

CND_num=sum(dc_node);
VND_num=sum(dv_node);
Dc=fliplr(Dci);
Dv=Dvi;
Frac_CND_num=fliplr(dc_node);
Frac_VND_num=dv_node;
Total_edge_CND=sum(Frac_CND_num.*Dc);
Total_edge_VND=sum(Frac_VND_num.*Dv);

% Test total_edge_number
diff=Total_edge_CND-Total_edge_VND;
disp('Difference between Total_CND_edge and Total_VND_edge');
disp(diff);

% CND node Initializaiton
CND_node_num=zeros(2, Total_edge_CND, 'int32');
edge_index=1;
i=1; 
remain=Frac_CND_num(i);
for node_num=1:1:CND_num,
    if (remain==0)
        %disp(edge_index);
        i=i+1;
        remain=Frac_CND_num(i);
    end
    CND_node_num(1,edge_index:edge_index+Dc(i)-1)=node_num;
    CND_node_num(2,edge_index:edge_index+Dc(i)-1)=Dc(i);
    edge_index=edge_index+Dc(i);
    remain=remain-1;
end

% VND node Initializaiton
VND_node_num=zeros(2, Total_edge_VND, 'int32');
edge_index=1;
i=1; 
remain=Frac_VND_num(i);
for node_num=1:1:VND_num,
    if (remain==0)
        disp(edge_index);
        i=i+1;
        remain=Frac_VND_num(i);
    end
    VND_node_num(1,edge_index:edge_index+Dv(i)-1)=node_num;
    VND_node_num(2,edge_index:edge_index+Dv(i)-1)=Dv(i);
    edge_index=edge_index+Dv(i);
    remain=remain-1;
end

% Deinterleaver Initialization
Deintlv=zeros(2, Total_edge_CND, 'int32');
% 1st row index X, value Y : edge 從 CND 的第Xth edge, 連到 VND的第Yth edge
% 2nd row index X, value Y : VND 第Yth edge, 有沒有被連接, 0 代表沒有被連接

% Deinterleaver Generation !!!
S_CND=max(Dc);   % Adaptively, 隨著CND的degree 做調整
S_VND=20;        % Note VND's degree is 10 

CND_index=1;        % initial
edge_index=1;       % initial
pre_degree=Dc(1);   % initial
CND_equ_table=zeros(1,70);   % initial
CND_equ_node =[];   % initial
Error_Event=[];     % initial 

while ( CND_index <= CND_num )
    
    % 如果degree改變，要把talbe清空
    now_degree=CND_node_num(2, edge_index); 
    if ( now_degree ~= pre_degree )
        CND_equ_table=zeros(1,now_degree);
        if (25 <now_degree)
            S_CND=now_degree;
            %disp('now S_CND is');
            %disp(S_CND);
        else
            S_CND=20;
            %disp('now S_CND is');
            %disp(S_CND);
        end
    end
    
    fail=0;
    search=1; % 執行
    
    while (search==1)
        search=0;   % 一旦遇到不符合的情形，馬上把search=1反過來-->重新找過
        
        pool=find( Deintlv(2,:) == 0 ); % =0 代表沒有被連接, pool 代表還可以選的連線
        
        % Purpose: 選出 now_degree 個 candidate
        if (length(pool) > 400), % 可以一次選擇所有CND需要的node, 重複機率低, !!!!!!!!!!!!!!!!! can be tuned
            candidate=pool( ceil(rand([1,now_degree]).*length(pool)));
        else % 用 randperm 去找
            position=randperm(length(pool));
            candidate=pool(position(1:now_degree));
        end
               
        % 暫時先假設 candidate 是 OK 的
        Deintlv(1, edge_index:1:edge_index+now_degree-1 ) = candidate;     % 從 CND 第 edge_index 連接到第 candidate VND edge
        Deintlv(2, candidate ) = 1;                                        % 從 VND 第 candidate edge 已經被連接
        
        % 找出所對應到第幾個node,而非edge
        CND_equ_node=VND_node_num(1, candidate);    % int32 form
        CND_equ_node=sort( CND_equ_node );          % int32 form
               
        % 檢驗 parallel edges, design principle (1)
        parallel_check= length(CND_equ_node)-length(unique(CND_equ_node));
               
        if (parallel_check > 0), % 代表 parallel edge 發生, 重新找一次
            disp(CND_index);
            disp(' parallel edge error ');
            Deintlv(1, edge_index:1:edge_index+now_degree-1) = 0;
            Deintlv(2, candidate ) = 0;
            search=1;    % search again
            fail=fail+1;
            error=[CND_index, 1];   % violate design principle 1
            Error_Event=[Error_Event; error];
            save ('Error_LDGM_wei.mat', 'Error_Event');
        else
            % 檢驗 same check node equation, design principle (3)
            Equ_check=ismember(CND_equ_node, CND_equ_table, 'rows');
            if (Equ_check == 1),% 代表有重複的 CND equation 發生, 重新找一次
                disp(CND_index);
                disp(' same CND equation error ');
                Deintlv(1, edge_index:1:edge_index+now_degree-1) = 0;
                Deintlv(2, candidate ) = 0;
                search=1;   % search again
                fail=fail+1;
                error=[CND_index, 3];   % violate design principle 3
                Error_Event=[Error_Event; error];
                save ('Error_LDGM_wei.mat', 'Error_Event');
            else % 目前通過 (1) 與 (3) 的測試
                 % 檢驗 S-random, design principle (2) 
                if (edge_index > max(Dc)),
                    pre_edge=[];
                    now_edge=[];
                    S_range=edge_index-S_CND+1;
                    for a=0:1:(now_degree-1),
                        pre_edge=[pre_edge,(edge_index-1+a):-1:(S_range+a)];
                        now_edge=[now_edge, repmat(edge_index+a, 1, S_CND-1)  ];
                    end
                    
                    if( ismember(1, abs(Deintlv(1,pre_edge) - Deintlv(1, now_edge)) <=S_VND ) ),
                        disp(CND_index);
                        disp(' S-random fail');
                        search=1; % 要重新再找一次
                        Deintlv(1, edge_index:1:edge_index+now_degree-1) = 0;
                        Deintlv(2, candidate ) = 0;
                        fail=fail+1;
                        error=[CND_index, 2];   % violate design principle 2
                        Error_Event=[Error_Event; error];
                        save ('Error_LDGM_wei.mat', 'Error_Event');
                    end
                else
                    S_range=1;
                    for now_edge=edge_index:1:(edge_index+now_degree-1),
                        pre_edge=now_edge-1:-1:S_range;
                        if( ismember(1, abs(Deintlv(1,pre_edge) - Deintlv(1, now_edge)) <=S_VND )),
                            disp(CND_index);
                            disp(' S-random fail');
                            search=1; % 要重新再找一次
                            Deintlv(1, edge_index:1:edge_index+now_degree-1) = 0;
                            Deintlv(2, candidate ) = 0;
                            fail=fail+1;
                            error=[CND_index, 2];   % violate design principle 2
                            Error_Event=[Error_Event; error];
                            save ('Error_LDGM_wei.mat', 'Error_Event');
                            disp("fail")
                        end
                    end
                end
            end            
        end
        
        if (fail>500)
            save LDGM_BPSK_wei;
            error('stop'); 
        end
    end
    
    CND_equ_table=[CND_equ_table; CND_equ_node];    % !!! 
    CND_equ_node=[];
    CND_index=CND_index+1;
    pre_degree=now_degree;
    edge_index=edge_index+now_degree;
    
end

save LDGM_BPSK_wei;