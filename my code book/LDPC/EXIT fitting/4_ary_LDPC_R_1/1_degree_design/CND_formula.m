% Ref: Brink_TCOM_04 
function [ Ie_CND ] = CND_formula( Ia_CND, dc )
            
    if(dc==1)
        disp ('Do not use degree 1 CND');
        error('stop');
    else
        % Target function is as below: formula (9) in ref
        % I=1-J_fun( sqrt(dc-1) * inv_J(1-Ia_CND) );
        
        %--- inv_J part
        inv_J_1=1-Ia_CND;
        
        I_asterisk=0.3646;
        a_sigma_1=1.09542;  b_sigma_1=0.214217;  c_sigma_1=2.33727;
        a_sigma_2=0.706692; b_sigma_2=0.386013;  c_sigma_2=-1.75017;
        
        g1=find(inv_J_1<=I_asterisk);
        std_var_g1= a_sigma_1.*inv_J_1(g1).*inv_J_1(g1)+b_sigma_1.*inv_J_1(g1)+c_sigma_1.*sqrt(inv_J_1(g1));
                
        g2=find(inv_J_1>I_asterisk & inv_J_1<1);
        std_var_g2= -1*a_sigma_2.*log(b_sigma_2.*(1-inv_J_1(g2)))-c_sigma_2*inv_J_1(g2);
                
        g3=find(inv_J_1>=1);
        std_var_g3=10*ones(1,length(g3));
        
        % order:  g3, g2, g1
        std_var=[std_var_g3, std_var_g2, std_var_g1];        
               
        std_var_F=sqrt(dc-1) * std_var;
        
        %--- J function part
        stdvar_asterisk=1.6363;
        a_J_1=-0.0421061;  b_J_1=0.209252;  c_J_1=-0.00640081;
        a_J_2=0.00181491; b_J_2=-0.142675;  c_J_2=-0.0822054; d_J_2=0.0549608;
        
        g1=find(std_var_F <= stdvar_asterisk);
        I_g1=a_J_1*std_var_F(g1).^3+b_J_1*std_var_F(g1).^2+c_J_1*std_var_F(g1);
    
        g2=find(std_var_F > stdvar_asterisk & std_var_F <10);
        I_g2= 1- exp(a_J_2*std_var_F(g2).^3+b_J_2*std_var_F(g2).^2+c_J_2*std_var_F(g2)+d_J_2);
    
        g3=find(std_var_F >= 10);
        I_g3=ones(1, length(g3));
        % order:  g3, g2, g1
        Ie_CND=[I_g3, I_g2, I_g1];
        
        Ie_CND=1-Ie_CND;
    end
end
