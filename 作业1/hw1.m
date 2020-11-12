% 电力系统第一次作业
% Date:2020/10/9

%% 数据准备
mpc = loadcase('case14');
branch14=mpc.branch;
bus14=mpc.bus;
baseMVA=mpc.baseMVA;

%% （1）不考虑变压器变比，生成 IEEE 14 节点系统的节点不定导纳矩阵Y_0
Y0=zeros(15,15);
for j=1:size(branch14,1)
    %     Y0(branch14(j,1),branch14(j,2))=branch14(j,3)+branch14(j,4)*1i;
    y=1/(branch14(j,3)+branch14(j,4)*1i);
    M=zeros(15,1);
    M(branch14(j,1))=1;
    M(branch14(j,2))=-1;
    Y0=Y0+M*y*M';% 支路的贡献
    if branch14(j,5)>0
        y2=(0.5*branch14(j,5)*1i);
        M1=zeros(15,1);
        M1(branch14(j,1))=1;
        M1(15)=-1;
        Y0=Y0+M1*y2*M1';% 线路接地导纳1
        M2=zeros(15,1);
        M2(branch14(j,2))=1;
        M2(15)=-1;
        Y0=Y0+M2*y2*M2';% 线路接地导纳2
    end
end
for j=1:size(bus14,1)
    y=bus14(j,6)*1i/baseMVA;
    M=zeros(15,1);
    M(j)=1;
    M(15)=-1;   
    Y0=Y0+M*y*M';% 节点接地导纳
end
% Y0=Y0+(diag(bus14(:,5))+diag(bus14(:,6))*1i)/baseMVA;% 节点接地导纳
complexmatrixplot(Y0,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])
Y=Y0(1:end-1,1:end-1);
complexmatrixplot(Y,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])

save('data_1','Y')
%% （2）考虑变压器变比， 生成 IEEE 14 节点系统对应导纳矩阵Y
Y_widetilde=zeros(14,14);
for j=1:size(branch14,1)
%     Y0(branch14(j,1),branch14(j,2))=branch14(j,3)+branch14(j,4)*1i;
    y=1/(branch14(j,3)+branch14(j,4)*1i);
    M=zeros(14,1);
    if branch14(j,9)>0
        M(branch14(j,1))=1/branch14(j,9);
        M(branch14(j,2))=-1;
    else 
        M(branch14(j,1))=1;
        M(branch14(j,2))=-1;
    end
    Y_widetilde=Y_widetilde+M*y*M';
    Y_widetilde(branch14(j,1),branch14(j,1))=Y_widetilde(branch14(j,1),branch14(j,1))+branch14(j,5)*1i/2;
    Y_widetilde(branch14(j,2),branch14(j,2))=Y_widetilde(branch14(j,2),branch14(j,2))+branch14(j,5)*1i/2;
end
Y_widetilde=Y_widetilde+(diag(bus14(:,5))+diag(bus14(:,6))*1i)/baseMVA;
complexmatrixplot(Y_widetilde,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])
%% （3）对比 matpower 中的 makeYbus 程序...
[YBUS, YF, YT] = makeYbus(mpc);
Full_Ybus=full(YBUS);

% 可视化部分

complexmatrixplot(Full_Ybus,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])

complexmatrixplot(Full_Ybus-Y_widetilde,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])
%% （4）采用修正的形式求解变压器变比变化后的导纳矩阵
Y_widetilde2=Y;
for j=1:size(branch14,1)
%     Y0(branch14(j,1),branch14(j,2))=branch14(j,3)+branch14(j,4)*1i;
    if branch14(j,9)>0
        y=1/(branch14(j,3)+branch14(j,4)*1i);
        M=zeros(14,1);
        M2=zeros(14,1);
        
        M(branch14(j,1))=1/branch14(j,9);
        M(branch14(j,2))=-1;
        M2(branch14(j,1))=1;
        M2(branch14(j,2))=-1;
        Y_widetilde2=Y_widetilde2+M*y*M'-M2*y*M2';
    end
end
complexmatrixplot(Y_widetilde2,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])
complexmatrixplot(Y_widetilde2-Y_widetilde,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])

