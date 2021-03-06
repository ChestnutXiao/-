%% 高网作业第三小腿
% Date:2020/11/11
% Author:Cai Xiao
clear all

%% 读取数据
mpc = loadcase('case14');
[YBUS, YF, YT] = makeYbus(mpc);
Full_Ybus=full(YBUS);
Full_Ybus=Full_Ybus(14:-1:1,14:-1:1);
complexmatrixplot(Full_Ybus,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])

%节点    1  2  3  4  5  6  7  8  9  10  11 12 13  14 
I_input=[0 1.1 0 0.9 0 1.3 0 0.8 0 1.05 0  0  1.2 0]';

I_input=I_input(14:-1:1);
%% 求出边界节点上的等值支路和等值注入电流，探究是否使用 WARD 等值的求解差异  
V_1=Full_Ybus\I_input;

Y_EE=Full_Ybus(1:5,1:5);
Y_BB=Full_Ybus(6:9,6:9);
Y_EB=Full_Ybus(1:5,6:9);
Y_BE=Full_Ybus(6:9,1:5);
widetilde_Y_BB=Y_BB-Y_BE*inv(Y_EE)*Y_EB;
Y_2=[widetilde_Y_BB,Full_Ybus(6:9,10:14);Full_Ybus(10:14,6:9),Full_Ybus(10:14,10:14)];
complexmatrixplot(Y_2,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])

I_E=I_input(1:5);
I_B=I_input(6:9);
widetilde_I_B=I_B-Y_BE*inv(Y_EE)*I_E;

I_2=[widetilde_I_B;I_input(10:14)];
V_2=Y_2\I_2;

%% 求逆
I_input=I_input(14:-1:1);
V_1=V_1(end:-1:1);

Y=makeYbus(mpc);
Z=inv(Y);

%% 求出以节点13,14对地为端口的戴维南等值参数，包括等值阻抗和等值戴维南电动势
Z_eq=Z(13,13)+Z(14,14)-2*Z(13,14);
V_eq=V_1(13)-V_1(14);

%% 当支路（12,13）开断后，对（2）中求出的戴维南等值参数进行修正  
M_alpha=zeros(14,1);
M_alpha(12)=1;
M_alpha(13)=-1;

M_L=zeros(14,1);
M_L(13)=1;
M_L(14)=-1;
y_alpha=1/(mpc.branch(19,3)+mpc.branch(19,4)*1i);
Y_alphaalpha=1/(-y_alpha+M_alpha'*Z*M_alpha);
Z_Lalpha=M_L'*Z*M_alpha;
V_alpha=(Z*M_alpha)'*I_input;

Z_eq_2=Z_eq-Z_Lalpha*Y_alphaalpha*Z_Lalpha';
V_eq_2=V_eq-Z_Lalpha*Y_alphaalpha*V_alpha;