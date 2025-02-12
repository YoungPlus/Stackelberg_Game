%用yalmip的kkt命令
clear;
clc;
%% 参数
price_day_ahead=[0.35;0.33;0.3;0.33;0.36;0.4;0.44;0.46;0.52;0.58;0.66;0.75;0.81;0.76;0.8;0.83;0.81;0.75;0.64;0.55;0.53;0.47;0.40;0.37]; %日前市场电价
price_b=1.2*price_day_ahead; %实时市场购入电价
price_s=1.2*price_day_ahead; %实时市场卖出电价
lb=0.8*price_day_ahead; %充电电价Ce下界
ub=1.2*price_day_ahead; %充电电价Ce上界
%电动汽车可充电时段
T_1=[1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1];
T_2=[1;1;1;1;1;1;1;1;0;0;0;0;1;1;1;0;0;0;0;1;1;1;1;1];
T_3=[0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0];
index1=find(T_1==0);index2=find(T_2==0);index3=find(T_3==0);
%% 定义变量
%上层决策变量
Ce=sdpvar(24,1); %电价
Pb=sdpvar(24,1); %日前购电
Pb_day=sdpvar(24,1); %实时购电
Ps_day=sdpvar(24,1); %实时售电
Pch=sdpvar(24,1); %储能充电
Pdis=sdpvar(24,1); %储能放电
z=binvar(24,1); %购售电状态
u=binvar(24,1); %储能状态
%下层决策变量
Pc1=sdpvar(24,1); %一类车充电功率
Pc2=sdpvar(24,1); %二类车充电功率
Pc3=sdpvar(24,1); %三类车充电功率
%状态变量
S=sdpvar(24,1); %储荷容量
for t=2:24
    S(t)=S(t-1)+0.9*Pch(t)-Pdis(t)/0.9; %储能充放电效率eta=0.9
end
%% 内层模型
CI=[sum(Pc1)==50*(0.9*24-9.6),sum(Pc2)==20*(0.9*24-9.6),sum(Pc3)==10*(0.9*24-9.6),Pc1>=0,Pc2>=0,Pc3>=0,Pc1<=50*3,Pc2<=20*3,Pc3<=10*3,Pc1(index1)==0,Pc2(index2)==0,Pc3(index3)==0]; %电量需求约束,Cons(10)-(12)
OI=sum(Ce.*(Pc1+Pc2+Pc3)); %内层目标函数
ops=sdpsettings('solver','gurobi','kkt.dualbounds',0);
[K,details] = kkt(CI,OI,Ce,ops); %内层模型最优性条件。建立KKT系统，其中Ce为参量
%% 外层模型
CO=[lb<=Ce<=ub,mean(Ce)==0.5,Pb>=0,0<=Pb_day<=1000*z,0<=Ps_day<=Pdis.*(1-z),0<=Pch<=1000*u,0<=Pdis<=1000*(1-u)]; %边界约束
CO=[CO,Pc1+Pc2+Pc3+Pch-Pdis==Pb+Pb_day-Ps_day]; %能量平衡
CO=[CO,sum(0.9*Pch-Pdis/0.9)==0,S(24)==2500,S>=0,S<=5000]; %SOC约束
OO=-(details.b'*details.dual+details.f'*details.dualeq)+sum(price_s.*Ps_day-price_day_ahead.*Pb-price_b.*Pb_day); %外层目标函数
%% 模型求解
optimize([K,CI,CO,boundingbox([CI,CO]),details.dual<=1],-OO)
