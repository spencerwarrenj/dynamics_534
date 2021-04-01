% Project.m
%% task 1
% get all data in SI units

%mass
lb_to_kg = 0.453592;
m_cm = 9730*lb_to_kg;
m_sm = 9690*lb_to_kg;
m_prop = 37295*lb_to_kg;
tot_mass = m_cm+m_sm+m_prop;

%com
in_to_m = 0.0254;
com_cm = [1043.1;-.1;7.8].*in_to_m;
com_sm = [908.2;0.7;-0.6].*in_to_m;
com_prop = [905.9;5.6;-2.4].*in_to_m;

%Inertia
slug_ft2_to_kg_m2 = 1.35581795;
I_cm = [4474,0,0;0,3919,0;0,0,3684].*slug_ft2_to_kg_m2;
I_sm = [6222,0,0;0,10321,0;0,0,10136].*slug_ft2_to_kg_m2;
I_prop = [19162,0,0;0,19872,0;0,0,26398].*slug_ft2_to_kg_m2;

%% Task A.a General Case
%a)find the center of mass as a vector relative to the A frame
r_csmp = (com_cm*m_cm+com_sm*m_sm+com_prop*m_prop)/tot_mass

%b)find total intertia matrix of csm with prop about B at the com

%find inertia for cm
d_cm = com_cm-r_csmp;
I_cm_o = I_cm+m_cm*(d_cm'*d_cm*eye(3)-d_cm*d_cm');

%find inertia for sm
d_sm = com_sm-r_csmp;
I_sm_o = I_sm+m_sm*(d_sm'*d_sm*eye(3)-d_sm*d_sm');

%find inertia for prop
d_prop = com_prop-r_csmp;
I_prop_o = I_prop+m_prop*(d_prop'*d_prop*eye(3)-d_prop*d_prop');

%find total inertia
I_csmp_o = I_cm_o+I_sm_o+I_prop_o