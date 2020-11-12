function [U,JU,IU,L,IL,JL,D]=cxLDU(A)
% A为矩阵的输入，要求为方阵
% [U,JU,IU,L,IL,JL,D]为稀疏格式的数据输出
[m,n]=size(A);
if m~=n
    disp('矩阵不为方阵！');
end

[row,col]=find(A);
U=[];
JU=[];
IU=[1];
L=zeros(1,1);
IL=zeros(1,1);
JL=zeros(1,1);
%% 计算下三角部分
for i=1:m
    JU_temp=sort(col(row==i));
    JU_temp=JU_temp(JU_temp>i);
    JU=[JU,JU_temp'];
    IU=[IU,size(JU,2)+1];
    Utemp=A(i,JU_temp);
    U=[U,Utemp];
end

%% 计算上三角部分
IL=JU;
JL=IU;
L=U;


%% 对角线部分
D=diag(A);