close all;
clear all;
clc;
Z_series=[7.2864+70.329i 5.5861+34.418i 5.5849+29.195i 5.56+35.531i 5.5331+30.664i;
   5.5861+34.418i 7.2864+70.329i 5.5861+34.418i 5.56+35.531i 5.56+35.531i;
   5.5849+29.195i 5.5861+34.418i 7.2864+70.329i 5.5331+30.664i 5.56+35.531i
   5.56+35.531i 5.56+35.531i 5.5331+30.664i 149.28+91.839i 5.5337+34.476i;
   5.5331+30.664i 5.56+35.531i 5.56+35.531i 5.5337+34.476i 149.28+91.839i];

Z_4x4=[Z_series(1,1)-(Z_series(1,5)*Z_series(5,1)/Z_series(5,5)) Z_series(1,2)-(Z_series(1,5)*Z_series(5,2)/Z_series(5,5)) Z_series(1,3)-(Z_series(1,5)*Z_series(5,3)/Z_series(5,5)) Z_series(1,4)-(Z_series(1,5)*Z_series(5,4)/Z_series(5,5));
    Z_series(2,1)-(Z_series(2,5)*Z_series(5,1)/Z_series(5,5)) Z_series(2,2)-(Z_series(2,5)*Z_series(5,2)/Z_series(5,5)) Z_series(2,3)-(Z_series(2,5)*Z_series(5,3)/Z_series(5,5)) Z_series(2,4)-(Z_series(2,5)*Z_series(5,4)/Z_series(5,5));
    Z_series(3,1)-(Z_series(3,5)*Z_series(5,1)/Z_series(5,5)) Z_series(3,2)-(Z_series(3,5)*Z_series(5,2)/Z_series(5,5)) Z_series(3,3)-(Z_series(3,5)*Z_series(5,3)/Z_series(5,5)) Z_series(3,4)-(Z_series(3,5)*Z_series(5,4)/Z_series(5,5));
    Z_series(4,1)-(Z_series(4,5)*Z_series(5,1)/Z_series(5,5)) Z_series(4,2)-(Z_series(4,5)*Z_series(5,2)/Z_series(5,5)) Z_series(4,3)-(Z_series(4,5)*Z_series(5,3)/Z_series(5,5)) Z_series(4,4)-(Z_series(4,5)*Z_series(5,4)/Z_series(5,5))];

Z_3x3 =[Z_4x4(1,1)-(Z_4x4(1,4)*Z_4x4(4,1)/Z_4x4(4,4)) Z_4x4(1,2)-(Z_4x4(1,4)*Z_4x4(4,2)/Z_4x4(4,4)) Z_4x4(1,3)-(Z_4x4(1,4)*Z_4x4(4,3)/Z_4x4(4,4));
        Z_4x4(2,1)-(Z_4x4(2,4)*Z_4x4(4,1)/Z_4x4(4,4)) Z_4x4(2,2)-(Z_4x4(2,4)*Z_4x4(4,2)/Z_4x4(4,4)) Z_4x4(2,3)-(Z_4x4(2,4)*Z_4x4(4,3)/Z_4x4(4,4));
        Z_4x4(3,1)-(Z_4x4(3,4)*Z_4x4(4,1)/Z_4x4(4,4)) Z_4x4(3,2)-(Z_4x4(3,4)*Z_4x4(4,2)/Z_4x4(4,4)) Z_4x4(3,3)-(Z_4x4(3,4)*Z_4x4(4,3)/Z_4x4(4,4))]
    
    Z_Z=[(Z_3x3(1,1)+Z_3x3(2,2)+Z_3x3(3,3))/3]
    Z_M=[(Z_3x3(1,2)+Z_3x3(1,3)+Z_3x3(2,3))/3]
    Z_seq = [(Z_3x3(1,1)+Z_3x3(2,2)+Z_3x3(3,3))/3+2*((Z_3x3(1,2)+Z_3x3(1,3)+Z_3x3(2,3))/3) 0 0;     % O sequence
             0 (Z_3x3(1,1)+Z_3x3(2,2)+Z_3x3(3,3))/3-((Z_3x3(1,2)+Z_3x3(1,3)+Z_3x3(2,3))/3) 0;       % positive sequence
             0 0 (Z_3x3(1,1)+Z_3x3(2,2)+Z_3x3(3,3))/3-((Z_3x3(1,2)+Z_3x3(1,3)+Z_3x3(2,3))/3)]      % negative sequence
     
     Z_p_o = [(Z_3x3(1,1)+Z_3x3(2,2)+Z_3x3(3,3))/3-((Z_3x3(1,2)+Z_3x3(1,3)+Z_3x3(2,3))/3)]
     
     Z_n_o = [(Z_3x3(1,1)+Z_3x3(2,2)+Z_3x3(3,3))/3-((Z_3x3(1,2)+Z_3x3(1,3)+Z_3x3(2,3))/3)]
      
     
     Z_o_o = [(Z_3x3(1,1)+Z_3x3(2,2)+Z_3x3(3,3))/3+2*((Z_3x3(1,2)+Z_3x3(1,3)+Z_3x3(2,3))/3)]
              
          
    V_base = 230e3;
    S_base = 100e6;

     
     Zb = (V_base^2)/S_base
     Z_pu = Z_seq/Zb
     
     Z_p = Z_p_o./Zb  % calculating p.u. sequence matrix
     Z_n = Z_n_o./Zb
     Z_o = Z_o_o./Zb
     
     %% power flow
     
     Y_p = inv(Z_p)
     Y_n = inv(Z_n)
     Y_o = inv(Z_o)
     
     Y_seq  = inv(Z_seq)
     Z_t = 0.1i;
     Y_t = inv(Z_t);
     Y = [Y_t+Y_p(1,1) -Y_p(1,1);-Y_p(1,1) Y_p(1,1)]
     G=real(Y(1,2))
     B=imag(Y(1,2))
     V1_1 = 1; V2_1 = 1; th1_1 = 0; th2_1 = 0; G21 = G; G22 = -G; B21 =B; B22 = -B;
 
     P_calc_1 = V2_1*V1_1*((G21*cos(th2_1-th1_1))+(B21*sin(th2_1-th1_1))) + V2_1*V2_1*((G22*cos(th2_1-th2_1))+(B22*sin(th2_1-th2_1)))
     Q_calc_1 = V2_1*V1_1*((G21*sin(th2_1-th1_1))-(B21*cos(th2_1-th1_1))) + V2_1*V2_1*((G22*sin(th2_1-th2_1))-(B22*cos(th2_1-th2_1)))
     
     P_gen_1 = 0; Q_gen_1 = 0; P_load_1 = 1.058; Q_load_1 = 0;
     
     P2_net_1 = P_gen_1 - P_load_1
     Q2_net_1 = Q_gen_1 - Q_load_1
     P2_delta_1 = P2_net_1 - P_calc_1
     Q2_delta_1 = Q2_net_1 - Q_calc_1
     
     J111=-Q_calc_1-V2_1^2*B22;
     J112=P_calc_1+V2_1^2*G22;
     J121=P_calc_1-V2_1^2*G22;
     J122=Q_calc_1-V2_1^2*B22;
     J1=[J111 J112;J121 J122]
     a= inv(J1)*[P2_delta_1;Q2_delta_1]
     th2_delta_1=a(1)
     V2_delta_1=a(2)
     th2_2 = th2_1 + th2_delta_1
     V2_2 = V2_1 + V2_delta_1
    
     %2nd Iteration
     P_calc_2 = V2_2*V1_1*((G21*cos(th2_2-th1_1))+(B21*sin(th2_2-th1_1))) + V2_2*V2_2*((G22*cos(th2_2-th2_2))+(B22*sin(th2_2-th2_2)))
     
     c=V2_2*V1_1*((G21*sin(th2_2-th1_1))-(B21*cos(th2_2-th1_1)))
     d=V2_2*V2_2*((G22*sin(th2_2-th2_2))-(B22*cos(th2_2-th2_2)))
     Q_calc_2 = c+d
     
    
     P2_delta_2 = P2_net_1 - P_calc_2
     Q2_delta_2 = Q2_net_1 - Q_calc_2
     
     J211=-Q_calc_2-V2_2^2*B22;
     J212=P_calc_2+V2_2^2*G22;
     J221=P_calc_2-V2_2^2*G22;
     J222=Q_calc_2-V2_2^2*B22;
     J=[J211 J212;J221 J222]
     b=inv(J)
     a2= b*[P2_delta_2;Q2_delta_2]
     th2_delta_2=a2(1)
     V2_delta_2=a2(2)
     th2_3=th2_2 + th2_delta_2
     V2_3=V2_2 + V2_delta_2
     
     %3rd Iteration
     P_calc_3 = V2_3*V1_1*((G21*cos(th2_3-th1_1))+(B21*sin(th2_3-th1_1))) + V2_3*V2_3*((G22*cos(th2_3-th2_3))+(B22*sin(th2_3-th2_3)))
     
     Q_calc_3 = V2_3*V1_1*((G21*sin(th2_3-th1_1))-(B21*cos(th2_3-th1_1)))+ V2_3*V2_3*((G22*sin(th2_3-th2_3))-(B22*cos(th2_3-th2_3)))
     
     
    
     P2_delta_3 = P2_net_1 - P_calc_3
     Q2_delta_3 = Q2_net_1 - Q_calc_3
     
     J311=-Q_calc_3-V2_3^2*B22;
     J312=P_calc_3+V2_3^2*G22;
     J321=P_calc_3-V2_3^2*G22;
     J322=Q_calc_3-V2_3^2*B22;
     J3=[J311 J312;J321 J322]
     b2=inv(J3)
     a3= b2*[P2_delta_3;Q2_delta_3]
     th2_delta_3=a3(1)
     V2_delta_3=a3(2)
     th2_4=th2_3 + th2_delta_3
     V2_4=V2_3 + V2_delta_3
     
     %4th Iteration
     
     %4th Iteration
     P_calc_4 = V2_4*V1_1*((G21*cos(th2_4-th1_1))+(B21*sin(th2_4-th1_1))) + V2_4*V2_4*((G22*cos(th2_4-th2_4))+(B22*sin(th2_4-th2_4)))
     
     Q_calc_4 = V2_4*V1_1*((G21*sin(th2_4-th1_1))-(B21*cos(th2_4-th1_1)))+ V2_4*V2_4*((G22*sin(th2_4-th2_4))-(B22*cos(th2_4-th2_4)))
         
    
     P2_delta_4 = P2_net_1 - P_calc_4
     Q2_delta_4 = Q2_net_1 - Q_calc_4
     
     J411=-Q_calc_4-V2_4^2*B22;
     J412=P_calc_4+V2_4^2*G22;
     J421=P_calc_4-V2_4^2*G22;
     J422=Q_calc_4-V2_4^2*B22;
     J4=[J411 J412;J421 J422]
     b3=inv(J4)
     a4= b3*[P2_delta_4;Q2_delta_4]
     th2_delta_4=a4(1)
     V2_delta_4=a4(2)
     th2_5=th2_4 + th2_delta_4
     V2_5=V2_4 + V2_delta_4
     
     
      %% Fault Current Calculation
     
     Ts = [1 1 1; 1 -0.5-0.866i -0.5+0.866i; 1 -0.5+0.866i -0.5-0.866i];
     Ib = 251;
     % 3-Phase to Ground Fault: Only Positive Sequence Equivalent Circuit
     % is required.
     Zg = 0; %+ve seq values of gen, transmission line, transformer and fault impedance.
     Z_tl_1 = Z_pu(2,2); % Z_series in pu.
     Z_t = 0.1i;
     Zf = 0.01/Zb; %Zf in ohm per 100km
     E = 1; %pu
     Z1 = Zg + Z_tl_1 + Z_t;
     If_3ph_1 = (E/(Z1+Zf));
     If_3ph_2 = 0;
     If_3ph_0 = 0;
     If_3ph_seq = [If_3ph_0; If_3ph_1;If_3ph_2]
     If_3ph = (Ts*If_3ph_seq*Ib)/1000
     
     % 1-phase-to-ground fault
     Z_t1_0 = Z_pu(1,1);
     Z0=Z_t1_0+Z_t;
     Z_t1_2 = Z_pu(3,3);
     Z2 = Z_t1_2+Z_t;
     If_1ph_seq=1/(Z0+Z1+Z2+3*Zf)*[E;E;E]
     If_1ph = (Ts*If_1ph_seq*Ib)/1000
     
     % 2-phase fault
     If_2ph_seq=1/(Z1+Z2+Zf)*[0;E;-E]
     If_2ph = (Ts*If_2ph_seq*Ib)/1000
     
     % 2-phase-to-ground fault
     Z_th_2phg=(Z1+Zf)*(Z2+Zf)+(Z1+Zf)*(Z0+Zf)+(Z2+Zf)*(Z0+Zf);
     If_2phg_0=(-(Z2+Zf)*E)/Z_th_2phg;
     If_2phg_1=(((Z2+Zf)+(Z0+Zf))*E)/Z_th_2phg;
     If_2phg_2=(-(Z0+Zf)*E)/Z_th_2phg;
     If_2phg_seq=[If_2phg_0;If_2phg_1;If_2phg_2]
     If_2phg=(Ts*If_2phg_seq*Ib)/1000