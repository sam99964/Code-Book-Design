% Ref: Brink_TCOM_04 
function [ Ie_VND ] = VND_formula( Ia_VND, Ie_DET, dv )

    % Target function is as below: formula (22) in ref
    % Ie_VND=J_fun( sqrt( (dv-1) * inv_J(Ia_VND).^2 + inv_J(Ie_DET).^2 ) );
    
    %--- inv_J part
    % inv_J(Ia_VND)
    I_asterisk=0.3646;
    a_sigma_1=1.09542;  b_sigma_1=0.214217;  c_sigma_1=2.33727;
    a_sigma_2=0.706692; b_sigma_2=0.386013;  c_sigma_2=-1.75017;
    
    g1=find(Ia_VND<=I_asterisk);
    std_var_g1= a_sigma_1.*Ia_VND(g1).*Ia_VND(g1)+b_sigma_1.*Ia_VND(g1)+c_sigma_1.*sqrt(Ia_VND(g1));
        
    g2=find(Ia_VND>I_asterisk & Ia_VND<1);
    std_var_g2= -1*a_sigma_2.*log(b_sigma_2.*(1-Ia_VND(g2)))-c_sigma_2*Ia_VND(g2);
        
    g3=find(Ia_VND==1);
    std_var_g3=10*ones(1,length(g3));
    
    std_var=[std_var_g1, std_var_g2, std_var_g3];
    std_var_1=(dv-1)*(std_var.^2);
    %disp('size of std_var_1');    disp(size(std_var_1));
    % inv_J(Ie_DET)
    g1=find(Ie_DET<=I_asterisk);
    std_var_g1= a_sigma_1.*Ie_DET(g1).*Ie_DET(g1)+b_sigma_1.*Ie_DET(g1)+c_sigma_1.*sqrt(Ie_DET(g1));
        
    g2=find(Ie_DET>I_asterisk & Ie_DET<1);
    std_var_g2= -1*a_sigma_2.*log(b_sigma_2.*(1-Ie_DET(g2)))-c_sigma_2*Ie_DET(g2);
        
    g3=find(Ie_DET==1);
    std_var_g3=10*ones(1,length(g3));
    
    std_var=[std_var_g1, std_var_g2, std_var_g3];
    std_var_2=(std_var.^2);
    %disp('size of std_var_2');    disp(size(std_var_2));
    
    std_var_2=std_var_2.*ones(1, numel(std_var_1));
    %disp('size of std_var_2');    disp(size(std_var_2));
    std_var=(std_var_1+std_var_2).^0.5;
    
    %--- J function part
    stdvar_asterisk=1.6363;
	a_J_1=-0.0421061;  b_J_1=0.209252;  c_J_1=-0.00640081;
    a_J_2=0.00181491; b_J_2=-0.142675;  c_J_2=-0.0822054; d_J_2=0.0549608;
    
    g1=find(std_var <= stdvar_asterisk);
    I_g1=a_J_1*std_var(g1).^3+b_J_1*std_var(g1).^2+c_J_1*std_var(g1);
    
    g2=find(std_var > stdvar_asterisk & std_var <10);
    I_g2= 1- exp(a_J_2*std_var(g2).^3+b_J_2*std_var(g2).^2+c_J_2*std_var(g2)+d_J_2);
    
    g3=find(std_var >= 10);
    I_g3=ones(1, length(g3));
    
    Ie_VND=[I_g1, I_g2, I_g3];
    
end

