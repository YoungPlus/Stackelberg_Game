%论文复现
%通过KKT条件理论以及二进制展开法、大Ｍ法等线性化手段将双层均衡模型转化为混合整数线性规划
%段声志,陈皓勇,郑晓东,黄剑平,邓盛盛.碳市场背景下发电商竞价策略及电力市场均衡分析[J].电测与仪表,2022,59(05)
clear;
clc;
%% 给定参数
G=[1,1,3,4,5]; %机组所在节点
D=[0,2,3,4,0]; %负荷所在节点
P_Gmax=[40;170;200;520;600]; %机组中标电量上限
P_Gmax=[P_Gmax,P_Gmax,P_Gmax];
P_Gmax=P_Gmax/3;
P_Gmin=0; %机组中标电量下限
P_Dmax=[0;0;300;300;400]; %用户中标电量上限
P_Dmax=[P_Dmax,P_Dmax,P_Dmax];
P_Dmax=P_Dmax/3;
P_Dmin=0; %用户中标电量下限
lamda_G=[280,310,360;300,320,342;305,337,380;290,315,340;260,308,343]; %机组除碳成本外的边际发电成本
phi=[0.88,0.64,0.85,0.81,0.8]; %机组碳排放强度
lamda_D=[0,0,0;600,500,448;550,450,440;580,460,430;0,0,0]; %用户报价
delta_t=1; %考虑时间段为1h
Q_co2=[900,788]; %碳配额
%线路电纳
B=[ 0,          1/0.0281,   0,          1/0.0304,   1/0.0064;
    1/0.0281,   0,          1/0.0108,   0,          0;
    0,          1/0.0108,   0,          1/0.0297,   0;
    1/0.0304,   0,          1/0.0297,   0,          1/0.0297;
    1/0.0064,   0,          0,          1/0.0297,   0];
%线路功率传输极限
P_Lmax=inf(5,5);P_Lmax(1,2)=400;P_Lmax(2,1)=400;P_Lmax(4,5)=240;P_Lmax(5,4)=240;
%% 决策变量
alpha=sdpvar(5,3); %机组报价
P_G=sdpvar(5,3); %机组中标电量
P_D=sdpvar(5,3); %用户中标电量
delta=sdpvar(5,1); %节点电压相角
%上层决策变量
%下层决策变量

%% 上层模型
%目标函数
Obj=sum(alpha.*P_G*delta_t,'all')-sum(lamda_D.*P_D*delta_t,'all');
%约束条件
Cons=[P_Gmin<=P_G<=P_Gmax,P_Dmin<=P_D<=P_Dmax]; %Cons7-8
for i=1:5
    g=find(G==i);d=find(D==i);
    sum_3=0;
    for j=1:5
        sum_3=sum_3+B(i,j)*(delta(i)-delta(j));
    end
    Cons=[Cons,sum(P_D(d,:),'all')-sum(P_G(g,:),'all')+sum_3==0]; %Cons9
end
for i=1:5
    for j=1:5
        Cons=[Cons,B(i,j)*(delta(i)-delta(j))<=P_Lmax(i,j)]; %Cons10
    end
end
Cons=[Cons,delta(1)==0]; %Cons11

%Cons for MC
for i=1:5
    for j=1:3
        Cons=[Cons,alpha(i,j)==lamda_G(i,j)];
    end
end

%% 下层模型
%目标函数
%约束条件

%% 模型求解
ops=sdpsettings('solver','cplex');
optimize(Cons,Obj,ops);