function dmdt = LLG_Solver_DiffEq(t,m)
global alpha hd...
       max_hfx  max_hfy  max_hfz...
       b_f f_ramp a_f pwf...
       max_hx_clk_p  max_hy_clk_p  max_hz_clk_p...
       max_hx_clk_d  max_hy_clk_d  max_hz_clk_d...
       b_clk clk_ramp a_clk pw_clk...              
       max_hx_bias_p  max_hy_bias_p  max_hz_bias_p...
       max_hx_bias_d  max_hy_bias_d  max_hz_bias_d...
       b_bias bias_ramp a_bias pw_bias clk_full f_full bias_full;
%  persistent cc
%  if isempty(cc)
%      cc=1;
%  end
%  cc=cc+1;
mx = m(1);
my = m(2);
mz = m(3);
%%%%%%%%%          Pulse shaping the extenrnal torques        %%%%%%%%%%%%%
%%%%%%%%%                    External Field                       %%%%%%%%%
%--------------------------------------------------------------------------
%%%%% Before pulse
if t <= b_f                                           
hfx = 0;
hfy = 0;
hfz = 0;
%%%%% During turn-ON
elseif f_ramp~=0 && t > b_f && t <= (b_f+f_ramp)    
hfx = (max_hfx/f_ramp) * (t-b_f);
hfy = (max_hfy/f_ramp) * (t-b_f);
hfz = (max_hfz/f_ramp) * (t-b_f);
%%%%% Pulse width
elseif  t > (b_f+f_ramp) && t <= (b_f+f_ramp+pwf)     
hfx = max_hfx;
hfy = max_hfy;
hfz = max_hfz;
%%%%% During turn-OFF
elseif f_ramp~=0 && t > (b_f+f_ramp+pwf) && t<= (b_f+2*f_ramp+pwf)
hfx = (max_hfx/f_ramp) * (f_ramp - (t-(b_f+f_ramp+pwf)) ) * f_full;
hfy = (max_hfy/f_ramp) * (f_ramp - (t-(b_f+f_ramp+pwf)) ) * f_full;
hfz = (max_hfz/f_ramp) * (f_ramp - (t-(b_f+f_ramp+pwf)) ) * f_full;
%%%%% After turn-OFF
elseif  t > (b_f+2*f_ramp+pwf) && t<= (b_f+2*f_ramp+pwf+a_f)
hfx = 0;
hfy = 0;
hfz = 0;
end
%save('pulse','hfz','-append')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                 Clock Spin-Torque                       %%%%%%%%%
%--------------------------------------------------------------------------
%%%%% Before pulse
if t <= b_clk                                           
hx_clk_p = 0;
hy_clk_p = 0;
hz_clk_p = 0;
%%%%%%%%%%%%%
hx_clk_d = 0;
hy_clk_d = 0;
hz_clk_d = 0;
%%%%% During turn-ON
elseif  clk_ramp~=0 && t > b_clk && t <= (b_clk+clk_ramp)                  
hx_clk_p = (max_hx_clk_p/clk_ramp) * (t-b_clk);
hy_clk_p = (max_hy_clk_p/clk_ramp) * (t-b_clk);
hz_clk_p = (max_hz_clk_p/clk_ramp) * (t-b_clk);
%%%%%%%%%%%%
hx_clk_d = (max_hx_clk_d/clk_ramp) * (t-b_clk);
hy_clk_d = (max_hy_clk_d/clk_ramp) * (t-b_clk);
hz_clk_d = (max_hz_clk_d/clk_ramp) * (t-b_clk);
%%%%% Pulse width
elseif  t > (b_clk+clk_ramp) && t <= (b_clk+clk_ramp+pw_clk)  
hx_clk_p = max_hx_clk_p;
hy_clk_p = max_hy_clk_p;
hz_clk_p = max_hz_clk_p;
%%%%%%%%%%%%
hx_clk_d = max_hx_clk_d;
hy_clk_d = max_hy_clk_d;
hz_clk_d = max_hz_clk_d;
%%%%% During turn-OFF
elseif clk_ramp~=0 && t > (b_clk+clk_ramp+pw_clk) && t<= (b_clk+2*clk_ramp+pw_clk)
hx_clk_p = (max_hx_clk_p/clk_ramp) * (clk_ramp - (t-(b_clk+clk_ramp+pw_clk)) ) * clk_full;
hy_clk_p = (max_hy_clk_p/clk_ramp) * (clk_ramp - (t-(b_clk+clk_ramp+pw_clk)) ) * clk_full;
hz_clk_p = (max_hz_clk_p/clk_ramp) * (clk_ramp - (t-(b_clk+clk_ramp+pw_clk)) ) * clk_full;
%%%%%%%%%%%%
hx_clk_d = (max_hx_clk_d/clk_ramp) * (clk_ramp - (t-(b_clk+clk_ramp+pw_clk)) ) * clk_full;
hy_clk_d = (max_hy_clk_d/clk_ramp) * (clk_ramp - (t-(b_clk+clk_ramp+pw_clk)) ) * clk_full;
hz_clk_d = (max_hz_clk_d/clk_ramp) * (clk_ramp - (t-(b_clk+clk_ramp+pw_clk)) ) * clk_full;
%%%%% After turn-OFF
elseif  t > (b_clk+2*clk_ramp+pw_clk) && t<= (b_clk+2*clk_ramp+pw_clk+a_clk)
hx_clk_p = 0;
hy_clk_p = 0;
hz_clk_p = 0;
%%%%%%%%%%%%
hx_clk_d = 0;
hy_clk_d = 0;
hz_clk_d = 0;
end
%hz_clk_d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                 Bias Spin-Torque                        %%%%%%%%%
%--------------------------------------------------------------------------
%%%%% Before pulse
if t <= b_bias                                           
hx_bias_p = 0;
hy_bias_p = 0;
hz_bias_p = 0;
%%%%%%%%%%%%%%%
hx_bias_d = 0;
hy_bias_d = 0;
hz_bias_d = 0;
%%%%% During turn-ON
elseif bias_ramp~=0 && t > b_bias && t <= (b_bias+bias_ramp)                  
hx_bias_p = (max_hx_bias_p/bias_ramp) * (t-b_bias);
hy_bias_p = (max_hy_bias_p/bias_ramp) * (t-b_bias);
hz_bias_p = (max_hz_bias_p/bias_ramp) * (t-b_bias);
%%%%%%%%%%%%%%
hx_bias_d = (max_hx_bias_d/bias_ramp) * (t-b_bias);
hy_bias_d = (max_hy_bias_d/bias_ramp) * (t-b_bias);
hz_bias_d = (max_hz_bias_d/bias_ramp) * (t-b_bias);
%%%%% Pulse width
elseif  t > (b_bias+bias_ramp) && t <= (b_bias+bias_ramp+pw_bias)  
hx_bias_p = max_hx_bias_p;
hy_bias_p = max_hy_bias_p;
hz_bias_p = max_hz_bias_p;
%%%%%%%%%%%%%%
hx_bias_d = max_hx_bias_d;
hy_bias_d = max_hy_bias_d;
hz_bias_d = max_hz_bias_d;
%%%%% During turn-OFF
elseif bias_ramp~=0 && t > (b_bias+bias_ramp+pw_bias) && t<= (b_bias+2*bias_ramp+pw_bias)
hx_bias_p = (max_hx_bias_p/bias_ramp) * (bias_ramp - (t-(b_bias+bias_ramp+pw_bias)) ) * bias_full;
hy_bias_p = (max_hy_bias_p/bias_ramp) * (bias_ramp - (t-(b_bias+bias_ramp+pw_bias)) ) * bias_full;
hz_bias_p = (max_hz_bias_p/bias_ramp) * (bias_ramp - (t-(b_bias+bias_ramp+pw_bias)) ) * bias_full;
%%%%%%%%%%%%%%
hx_bias_d = (max_hx_bias_d/bias_ramp) * (bias_ramp - (t-(b_bias+bias_ramp+pw_bias)) ) * bias_full;
hy_bias_d = (max_hy_bias_d/bias_ramp) * (bias_ramp - (t-(b_bias+bias_ramp+pw_bias)) ) * bias_full;
hz_bias_d = (max_hz_bias_d/bias_ramp) * (bias_ramp - (t-(b_bias+bias_ramp+pw_bias)) ) * bias_full;
%%%%% After turn-OFF
elseif  t > (b_bias+2*bias_ramp+pw_bias) && t<= (b_bias+2*bias_ramp+pw_bias+a_bias)
hx_bias_p = 0;
hy_bias_p = 0;
hz_bias_p = 0;
%%%%%%%%%%%%%
hx_bias_d = 0;
hy_bias_d = 0;
hz_bias_d = 0;
end
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Internal Fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% uniaxial field
hux = 0;
huy = 0;
huz = mz;
%%%%%%%%%%%%%%%%%%%%
%%%%   Demagnetzing field  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdx = 0;             % This is the in-plane x-axis demagnetzing field
hdy = - hd*my;     % This is the OutOfPlane y demagnetzing field
hdz = 0;             % This is the in-plane z-axis demagnetzing field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%             Precession                               %%%%%%%%%%%
Phx = hux + hdx + hfx + hx_clk_p + hx_bias_p;
Phy = huy + hdy + hfy + hy_clk_p + hy_bias_p;
Phz = huz + hdz + hfz + hz_clk_p + hz_bias_p;
%%%%%%%%%%%
Pdmxdt = 1 *( my*Phz - mz*Phy);
Pdmydt = 1 *( mz*Phx - mx*Phz);
Pdmzdt = 1 *( mx*Phy - my*Phx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%             Damping                                  %%%%%%%%%%
Dhx = hux + hdx + hfx + hx_clk_d + hx_bias_d; 
Dhy = huy + hdy + hfy + hy_clk_d + hy_bias_d; 
Dhz = huz + hdz + hfz + hz_clk_d + hz_bias_d; 
mdh = mx*Dhx + my*Dhy + mz*Dhz;
mdm = mx*mx + my*my + mz*mz;
%fx=[Phx Dhx];
%%%%%%%%%%%
Ddmxdt = mdh*mx - Dhx*mdm;
Ddmydt = mdh*my - Dhy*mdm;
Ddmzdt = mdh*mz - Dhz*mdm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dmxdt = -(1/(2*alpha)) * Pdmxdt - (1/2) * Ddmxdt;
dmydt = -(1/(2*alpha)) * Pdmydt - (1/2) * Ddmydt;
dmzdt = -(1/(2*alpha)) * Pdmzdt - (1/2) * Ddmzdt;
%%%%%%%%%%%%%%%%%%%%%%%%  Solution   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dmdt  = [dmxdt; dmydt; dmzdt];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%