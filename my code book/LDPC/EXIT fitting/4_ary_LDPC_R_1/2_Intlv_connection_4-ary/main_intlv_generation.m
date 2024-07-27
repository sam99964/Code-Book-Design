clear all;clc;
% R=0.68, LDPC
% block length is 10^5
% interleaver connection principle
% (1) no parallel edges
% (2) S-random requirement
% (2) no two CND have the same check equations


Dvi=[2, 1, 2, 1, 2, 1, 2, 1, 2, 3, 7, 6, 7, 6, 7, 6, 7, 6, 7, 8, 11, 22, 23];
dv_node=[34736, 1, 35121, 1, 20896, 1, 17987, 1, 2648, 73495, 1769, 1, 9417, 1, 11418, 1, 11363, 1, 2903, 2583, 2, 4543, 4923];
         
Dci=[8];
dc_node=[116905];

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
if (diff ~=0)
    disp('Difference between Total_CND_edge and Total_VND_edge');
    disp(diff);
    error('stop'); 
end

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
        %disp(edge_index);
        i=i+1;
        remain=Frac_VND_num(i);
    end
    VND_node_num(1,edge_index:edge_index+Dv(i)-1)=node_num;
    VND_node_num(2,edge_index:edge_index+Dv(i)-1)=Dv(i);
    edge_index=edge_index+Dv(i);
    remain=remain-1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fileID = fopen('C:\Users\sam01\Desktop\python encode\VND_node_num.txt','w');
% fprintf(fileID,'%d\n',VND_node_num(1,:)-1);
% fclose(fileID);
% fileID = fopen('C:\Users\sam01\Desktop\python encode\CND_node_num.txt','w');
% fprintf(fileID,'%d\n',CND_node_num(1,:)-1);
% fclose(fileID);




%%%%%%%%%%%%%%%%%%%%%%%%%

% Deinterleaver Initialization
Deintlv=zeros(2, Total_edge_CND, 'int32');
% 1st row index X, value Y : edge �q CND ����Xth edge, �s�� VND����Yth edge
% 2nd row index X, value Y : VND ��Yth edge, ���S���Q�s��, 0 �N���S���Q�s��

% Deinterleaver Generation !!!
S_VND=46;           % Note max VND degree is 23
S_CND=32;           % Note max CND degree is 8

CND_index=1;        % initial
edge_index=1;       % initial
pre_degree=Dc(1);   % initial
CND_equ_table=zeros(1,8);   % initial
CND_equ_node =[];   % initial
Error_Event=[];     % initial 

while ( CND_index <= CND_num )
    
    % �p�Gdegree���ܡA�n��talbe�M��
    now_degree=CND_node_num(2, edge_index); 
    if ( now_degree ~= pre_degree )
        CND_equ_table=[];
    end
    
    fail=0;
    search=1; % ����
    
    while (search==1)
        search=0;   % �@���J�줣�ŦX�����ΡA���W��search=1�ϹL��-->���s��L
        
        pool=find( Deintlv(2,:) == 0 ); % =0 �N���S���Q�s��, pool �N���٥i�H�諸�s�u
        
        % Purpose: ��X now_degree �� candidate
        if (length(pool) > 400), % �i�H�@����ܩҦ�CND�ݭn��node, ���ƾ��v�C, !!!!!!!!!!!!!!!!! can be tuned
            candidate=pool( ceil(rand([1,now_degree]).*length(pool)));
        else % �� randperm �h��
            position=randperm(length(pool));
            candidate=pool(position(1:now_degree));
        end
               
        % �Ȯɥ����] candidate �O OK ��
        Deintlv(1, edge_index:1:edge_index+now_degree-1 ) = candidate;     % �q CND �� edge_index �s����� candidate VND edge
        Deintlv(2, candidate ) = 1;                                        % �q VND �� candidate edge �w�g�Q�s��
        
        % ��X�ҹ�����ĴX��node,�ӫDedge
        CND_equ_node=VND_node_num(1, candidate);    % int32 form
        CND_equ_node=sort( CND_equ_node );          % int32 form
               
        % ���� parallel edges, design principle (1)
        parallel_check= length(CND_equ_node)-length(unique(CND_equ_node));
               
        if (parallel_check > 0), % �N�� parallel edge �o��, ���s��@��
            disp(CND_index); disp(' parallel edge error ');
            Deintlv(1, edge_index:1:edge_index+now_degree-1) = 0;
            Deintlv(2, candidate ) = 0;
            search=1;    % search again
            fail=fail+1;
            error=[CND_index, 1];   % violate design principle 1
            Error_Event=[Error_Event; error];
            %save ('Error_LDGM_simple.mat', 'Error_Event');
        else
            % ���� same check node equation, design principle (3)
            Equ_check=ismember(CND_equ_node, CND_equ_table, 'rows');
            if (Equ_check == 1),% �N�������ƪ� CND equation �o��, ���s��@��
                disp(CND_index);
                disp(' same CND equation error ');
                Deintlv(1, edge_index:1:edge_index+now_degree-1) = 0;
                Deintlv(2, candidate ) = 0;
                search=1;   % search again
                fail=fail+1;
                error=[CND_index, 3];   % violate design principle 3
                Error_Event=[Error_Event; error];
                %save ('Error_LDGM_simple.mat', 'Error_Event');
            else % �ثe�q�L (1) �P (3) ������
                 % ���� S-random, design principle (2) 
                 if (edge_index > S_CND),
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
                        search=1; % �n���s�A��@��
                        Deintlv(1, edge_index:1:edge_index+now_degree-1) = 0;
                        Deintlv(2, candidate ) = 0;
                        fail=fail+1;
                        error=[CND_index, 2];   % violate design principle 2
                        Error_Event=[Error_Event; error];
                        %ave ('Error_LDGM_simple.mat', 'Error_Event');
                    end
                else
                    S_range=1;
                    for now_edge=edge_index:1:(edge_index+now_degree-1),
                        pre_edge=now_edge-1:-1:S_range;
                        if( ismember(1, abs(Deintlv(1,pre_edge) - Deintlv(1, now_edge)) <=S_VND )),
                            disp(CND_index);
                            disp(' S-random fail');
                            search=1; % �n���s�A��@��
                            Deintlv(1, edge_index:1:edge_index+now_degree-1) = 0;
                            Deintlv(2, candidate ) = 0;
                            fail=fail+1;
                            error=[CND_index, 2];   % violate design principle 2
                            Error_Event=[Error_Event; error];
                            %save ('Error_LDGM_simple.mat', 'Error_Event');
                        end
                    end
                end
            end            
        end
        
        if (fail>1000)
            save LDPC_BPSK_fail;
            error('stop'); 
        end
    end
    
    CND_equ_table=[CND_equ_table; CND_equ_node];    % !!! 
    CND_equ_node=[];
    CND_index=CND_index+1;
    pre_degree=now_degree;
    edge_index=edge_index+now_degree;
    
end

save LDPC_BPSK_simple;