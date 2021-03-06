%% 第二小题
% Date:2020/10/24
% Author:Cai Xiao
clear all
%% 读取历史数据
load('data_1.mat')
n=size(Y,1);
complexmatrixplot(Y,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])
Y_copy=Y;
%% 2.1.1 基于不考虑变压器变比的 IEEE 14 节点系统导纳矩阵Y ， 实现 LDU 分解
% LDU分解部分(不考虑稀疏）
ZeroJudge2=0;
ZeroJudge3=0;
Judge2=0;
Judge3=0;
for p=1:n
    for j=p+1:n
        ZeroJudge2=ZeroJudge2+1;% 检验非零元数目
        if Y(p,j)~=0
            Judge2=Judge2+1;% 检验非零元数目
            Y(p,j)=Y(p,j)/Y(p,p);
            for i=p+1:n
                ZeroJudge3=ZeroJudge3+1;% 检验非零元数目
                if Y(i,p)~=0
                    Judge3=Judge3+1;% 检验非零元数目
                    Y(i,j)=Y(i,j)-Y(i,p)*Y(p,j);
                end
            end
        end
    end
    for k=p+1:n
        if Y(k,p)~=0
            Y(k,p)=Y(k,p)/Y(p,p);
        end
    end
end
complexmatrixplot(Y,'ColorBar','On');% 可视化
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])
%% 检验注入元
zhuru=zeros(n,n);
for i=1:n
    for j=1:n
        if Y_copy(i,j)==0&&Y(i,j)~=0
            zhuru(i,j)=1;
        end
    end
end
complexmatrixplot(zhuru,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])

%% 基于 LDU 分解，采用稀疏矩阵技术，求解该导纳矩阵的逆（即阻抗矩阵Z ）


[U,JU,IU,L,IL,JL,D]=cxLDU(Y);
N=14;
Z=eye(N);
Z(N,N)=1/D(N);
for i=N-1:-1:1
    for j=N:-1:i+1
        for k=IU(i):IU(i+1)-1
            if JU(k)<j
                Z(i,j)=Z(i,j)-U(k)*Z(JU(k),j); 
            else
                Z(i,j)=Z(i,j)-U(k)*Z(j,JU(k)); 
            end
        end
    end
    Z(i,i)=1/D(i);
    for k=IU(i):IU(i+1)-1
        Z(i,i)=Z(i,i)-U(k)*Z(i,JU(k)); 
    end
end
complexmatrixplot(Z,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])

complexmatrixplot(inv(Y_copy),'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])
%%
% 面向支路
y_1=diag([1/(0.016+0.058*1i) 0.081*1i 0.081*1i]);
M_1=zeros(14,3);
M_1(5,1)=1;
M_1(6,1)=-1;
M_1(5,2)=1;
M_1(6,3)=1;
deltaY_1=M_1*y_1*M_1';% 支路的贡献
complexmatrixplot(deltaY_1,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])
% 面向节点
y_2=[0.081*1i+1/(0.016+0.058*1i) -1/(0.016+0.058*1i);-1/(0.016+0.058*1i) 0.081*1i+1/(0.016+0.058*1i)];
M_2=zeros(14,2);
M_2(5,1)=1;
M_2(6,2)=1;
deltaY_2=M_2*y_2*M_2';% 支路的贡献
% complexmatrixplot(deltaY_2,'ColorBar','On');
% colormap(flipud(hot)*0.6+[0.4 0.4 0.4])

complexmatrixplot(Y_copy,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])

Y_add=Y_copy+deltaY_1;

complexmatrixplot(Y_add,'ColorBar','On');
colormap(flipud(hot)*0.6+[0.4 0.4 0.4])