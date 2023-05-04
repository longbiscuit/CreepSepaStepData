%% 将分级加载蠕变数据整理为分别加载
clc;
clear all;
close all;
global iterGbestval;
global BestSol;
global GFES; %计数，看误差函数调用几次
global NumRDP;
BestSol.Position=[];
BestSol.Cost=[];

%% 需要自己修改的
q=[50000,100000,200000,300000];%每个阶段所施加的等效应力q
SepaTimeEnd=[23.76666667,49.39444444,74.60833333,96.35833333];%每个阶段结束时间点
SepaData=readmatrix("SepaData.txt");%输入分级加载txt文档的路径和名称
%% 需要自己修改的

fig1=figure('Name','Raw creep data with two column');
plot(SepaData(:,1),SepaData(:,2),'ko');hold on;
xlabel('{\it t}')
ylabel('{\it \epsilon}_{a}')
set(gca,'FontName','Times New Roman','FontSize',12)
[AllGroupTwoColMeasStartZero,AllGroupTwoColMeas,GroupEndPoint,sepaPointNum,EndIndexInSepaData]=GetSepaPointNum(SepaData,SepaTimeEnd);%每级多少个数据点

fig2=figure('Name','translate');
for ii=1:numel(SepaTimeEnd)
    plot(AllGroupTwoColMeasStartZero{ii}(:,1),AllGroupTwoColMeasStartZero{ii}(:,2),'bo');hold on;
end
xlabel('{\it t}')
ylabel('{\it \epsilon}_{a}')
set(gca,'FontName','Times New Roman','FontSize',12)

% 根据每段数据分别拟合参数
nVar=3;
VarSize=[1 nVar];
VarMin=[1e-27,   0.00,  -0.99];%
VarMax=[1.80,   3.00,  0.00];%
NumRDP=[ 27,      27,     27];%小数点后右侧有几位小数
fhd= str2func('TimeHardening');

pop_size=3^(4+floor(log2(sqrt(nVar))));%2^(4+floor(log2(sqrt(nVar)))); % Do not change this value
group_size=3; % Do not change this value
fes_max=100000*nVar %  can change this value
iter_max=ceil(fes_max/pop_size) % can change this value
group_num=10*nVar; % group_num=20+round(nVar/10); % can change this value

% Unused Settings
jingdu=0 % precision
% STADE algorithm was used for optimization
for ii=1:numel(SepaTimeEnd)

    [gbest,gbestval,FES,suc,suc_fes]= DE32_func(jingdu,fhd,nVar,...
    group_size,group_num,iter_max,fes_max,VarMin,VarMax,NumRDP,...
    AllGroupTwoColMeasStartZero{ii},q(ii));
    sepaPara(:,ii)=gbest
    if ii==1
        %AllGroupTwoColMeas{ii-1}(end,:) 平移前起始点
        tadd=SepaData(:,1);%取第一组参数
        %tadd=[AllGroupTwoColMeasStartZero{ii}(:,1);AllGroupTwoColMeasStartZero{ii+1}(:,1)+GroupEndPoint{ii}(1,1)];
        AllGroupTwoColModelStartZero{ii}=GetTimeHardening_t_ea(sepaPara(:,ii),q(ii),tadd);
    else
        tadd=SepaData(EndIndexInSepaData(ii-1):end,1)-GroupEndPoint{ii-1}(1,1);%第ii阶段及之后所有所有数据减去ii-1段结尾那行数据，相当于移到0点
        %tadd=[AllGroupTwoColMeasStartZero{ii}(:,1)];
        AllGroupTwoColModelStartZero{ii}=GetTimeHardening_t_ea(sepaPara(:,ii),q(ii),tadd);
    end
    
    if ii==1
        AllGroupTwoColModel{ii}=AllGroupTwoColModelStartZero{ii}+AllGroupTwoColMeas{ii}(1,:);%平移回去，当ii=1是起始点为0相当于不平移
    else
        AllGroupTwoColModel{ii}=AllGroupTwoColModelStartZero{ii}+AllGroupTwoColMeas{ii-1}(end,:);%正常等级加载t对应的预测值
    end
    
    figure(fig2)
    plot(AllGroupTwoColModelStartZero{ii}(:,1),AllGroupTwoColModelStartZero{ii}(:,2),'b.');hold on
    
    figure(fig1)%绘制根据全部或部分（使用全部数据时趋势可能不太好，所以有时用部分数据反而更好）试验数据拟合的理论值
    %下面绘图可以绘制全部，也可以仅绘制预测部分
    %plot(AllGroupTwoColModel{ii}(:,1),AllGroupTwoColModel{ii}(:,2),'k.');hold on;
     plot(AllGroupTwoColModel{ii}(sepaPointNum(ii):end,1),AllGroupTwoColModel{ii}(sepaPointNum(ii):end,2),'k.-');hold on;
   
end

% 开始相减
translate_ea{1}=AllGroupTwoColMeas{1};%最终转换的分别加载数据
writematrix(translate_ea{1,1},[num2str(q(1)) '.dat'])
for ii=2:numel(SepaTimeEnd)
        for jj=1:sepaPointNum(ii)
            minus_ea{ii}(jj,1)=AllGroupTwoColModel{ii}(jj,1);%预测线上x，也是下一组实测x（根据下一组实测x获得，所以不用插值，因为x是对应的）
            minus_ea{ii}(jj,2)=AllGroupTwoColMeas{ii}(jj,2)-AllGroupTwoColModel{ii-1}(sepaPointNum(ii-1)+jj-1,2);%实测y-预测延长线y
            translate_ea{ii}(jj,1)=AllGroupTwoColMeas{ii}(jj,1)-AllGroupTwoColMeas{ii}(1,1);%GroupEndPoint{ii-1}(1,1);%实测x-上一阶段结束点x进行平移
            translate_ea{ii}(jj,2)=interp1(translate_ea{ii-1}(:,1),translate_ea{ii-1}(:,2),minus_ea{ii}(jj,1)-GroupEndPoint{ii-1}(1,1),'linear','extrap')...
            +minus_ea{ii}(jj,2);
            
        end
    figure(fig1)
    plot(translate_ea{ii}(:,1),translate_ea{ii}(:,2),'r.');hold on;
    writematrix(translate_ea{1,ii},[num2str(q(ii)) '.dat'])
end
savefig(fig1,'fig1.fig')
savefig(fig2,'fig2.fig')
% writecell(translate_ea)
% type 'translate_ea.txt'





%% 1-根据分级加载数据找到每级加载所含试验点数,并返回每级加载数据
% 输入：
% SepaData 原始数据
% SepaTimeEnd 每级数据结束的时间点
%
% 输出：
% AllGroupTwoColMeas 将各级数据分开，每级一组
% AllGroupTwoColMeasStartZero 将各级数据分开，每级一组，且起始点移动到原点
% GroupEndPoint 每组最后一行数据点
% sepaPointNum 每级数据有多少个数据点
function [AllGroupTwoColMeasStartZero,AllGroupTwoColMeas,GroupEndPoint,sepaPointNum,EndIndexInSepaData]=GetSepaPointNum(SepaData,SepaTimeEnd)
sepaPointNum=zeros(numel(SepaTimeEnd),1);
for ii=1:numel(SepaTimeEnd)
    mm(ii)=find(SepaData(:,1)==SepaTimeEnd(ii));
    if ii==1
        AllGroupTwoColMeas{ii}=SepaData(1:mm(ii),:);
        AllGroupTwoColMeasStartZero{ii}=AllGroupTwoColMeas{ii}-AllGroupTwoColMeas{ii}(1,:)% 将数分级加载数据起始点归零
        
    else
        AllGroupTwoColMeas{ii}=SepaData(mm(ii-1):mm(ii),:);% 取第2级试验数据的第2行第1列：AllGroupTwoColMeas{2}(2,1)
        AllGroupTwoColMeasStartZero{ii}=AllGroupTwoColMeas{ii}-AllGroupTwoColMeas{ii-1}(end,:)% 将数分级加载数据起始点归零，方法是本阶段数据减去上阶段最后一行数据
    end
    GroupEndPoint{ii}=AllGroupTwoColMeas{ii}(end,:);
end
sepaPointNum(1)=mm(1);
sepaPointNum(2:end)=diff(mm)+1;

EndIndexInSepaData=mm;

end

%% 2-根据参数、q和t获得对应的ea
function [t_ea_model]=GetTimeHardening_t_ea(para,q,t)
t_ea_model(:,1)=t;
A=para(1);n=para(2);m=para(3);
t_ea_model(:,2)=A/(m+1)*q^n*t.^(m+1);
end


%% fitness TimeHardening
function z = TimeHardening( x,varargin )
global BestSol;
global GFES %计数，看误差函数调用几次
global NumRDP;
[~,n]=size(x);
z=zeros(1,n);
t_ea_Measure=varargin{1,1};% https://blog.csdn.net/weixin_43859329/article/details/103919801
q=varargin{1,2};

for jj=1:n%每列一组参数
    para1=x( : , jj )';
    %保留小数位数
    for ii=1:numel(para1)
        para1(ii)=round(para1(ii),NumRDP(ii));
    end
    
    if ~isempty(BestSol.Cost)% 如果BestSol.Cost不为空
        if all(para1==BestSol.Position)
            z(jj)=BestSol.Cost;
        else
        t_ea_Model=GetTimeHardening_t_ea (para1 ,q,t_ea_Measure(:,1));
%         z(jj)=sum((t_ea_Measure(:,2)-t_ea_Model(:,2)).^2);%使用全部数据
         [rowth,~]=size(t_ea_Measure);
        z(jj)=sum((t_ea_Measure(round(rowth/1.5):end,2)-t_ea_Model(round(rowth/1.5):end,2)).^2);
        GFES=GFES+1;
        end
    else% BestSol.Cost为空，说明是第一代计算
        t_ea_Model=GetTimeHardening_t_ea (para1 ,q,t_ea_Measure(:,1));
        z(jj)=sum((t_ea_Measure(:,2)-t_ea_Model(:,2)).^2);
        GFES=GFES+1;
    end
    
    if isinf(z(jj))%有时参数不协调（比如kappa>lambda）会出现z无穷大
        z(jj)=1e300;
    end
    if isnan(z(jj))
        z(jj)=1e300;
    end
end
end



% nVar=9;
% VarSize=[1 nVar];
% VarMin=[0.8,   0.05,  0.002,    0.02,    0.17,    0.0,   0.17,   0.000,  0.0];%
% VarMax=[1.8,   0.45,  0.100,    0.40,    3.00,     0.0,   1.50,   0.99,  10.0];%
% NumRDP=[ 3,      3,     3,       3,      3,      3,       3,     3,     3];%小数点后右侧有几位小数
% fhd= str2func('UserDefineCallMex');
% 
% pop_size=30;%2^(4+floor(log2(sqrt(D)))); % Do not change this value
% group_size=3; % Do not change this value
% fes_max=20000*nVar %  can change this value
% iter_max=ceil(fes_max/pop_size) % can change this value
% group_num=10*nVar; % group_num=20+round(nVar/10); % can change this value
% 
% % Unused Settings
% jingdu=1e-8 % precision
% para2=[1,0,0,ErrorType];%目标函数是否加约束，是否输出理论值,（0读入csuhReadInList.txt，1读入csuhReadInListZong.txt）
% % STADE algorithm was used for optimization
% [gbest,gbestval,FES,suc,suc_fes]= DE32_func(jingdu,fhd,nVar,...
%     group_size,group_num,iter_max,fes_max,VarMin,VarMax,NumRDP,para2);




%%%%%%%%%%%%%%%%%差分进化算法求函数极值%%%%%%%%%%%%%%%%%
function [gbest,gbestval,FES,suc,suc_fes]=DE32_func(jingdu,fhd,Dimension,group_size,group_num,Max_Gen,Max_FES,Xmin,Xmax,NumRDP,varargin)
% profile on
global iterGbestval;
global BestSol;
global GFES;
GFES=0;
rand('state',sum(100*clock));

FES=0;
suc=0;
suc_fes=0;
counter=0;
% SETTING
randKind='UR'% UR:Uniform Random; SB:SOBOL random
useLocalSTA=1% 1:is use sta local search; 0:NOT use sta local search
useDynamicParameter=1% 1:is use Dynamic Parameter F and CR; 0:NOT use
useDNP=1;%是否用更多的个体寻找更好的初始解
%%%%%%%%%%%%%%%%%%%%%%%%%初始化%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;                            %清除所有变量
% close all;                            %清图
% clc;                                  %清屏
% NP=20;                                %个体数目
% D=2;                                  %变量的维数
% G=100;                                %最大进化代数
% F=0.5;                                %变异算子
% CR=0.1;                               %交叉算子
% Xs=4;                                 %上限
% Xx=-4;                                %下限

NP=group_num;                                %个体数目
D=Dimension;                                  %变量的维数
G=Max_Gen;                                %最大进化代数
F0=0.5;                                %变异算子F=[0,2]  F=0.5;
% CR=0.1;                               %交叉算子CR=[0,1] CR=0.1;
Xs=Xmax;                                 %上限
Xx=Xmin;                                %下限

%%%%%%%%%%%%%%%%%%%%%%%%%赋初值%%%%%%%%%%%%%%%%%%%%%%%%
if useDNP==1
    tempNP=NP*D*10;
else
    tempNP=NP
end

x=zeros(D,NP);                        %初始种群
v=zeros(D,NP);                        %变异种群
u=zeros(D,NP);                        %选择种群
xini=zeros(D,tempNP);


% randKind='SB';% UR:Uniform Random; SB:SOBOL random
switch randKind
    case 'UR'
        xini=rand(D,tempNP).*(Xs-Xx)'+Xx';              %赋初值,每列是一组参数
    case 'SB'
        p = sobolset(D);
        xini=p(1:tempNP,:)'.*(Xs-Xx)'+Xx';
    otherwise
        xini=rand(D,tempNP).*(Xs-Xx)'+Xx';              %赋初值,,每列是一组参数
end

%%%%%%%%%%%%%%%%%%%%计算目标函数%%%%%%%%%%%%%%%%%%%%%%%
for m=1:tempNP
    %保留小数位数
    for num=1:D
        xini(num,m)=round(xini(num,m),NumRDP(num));
    end
    
    Ob(m)=feval(fhd,xini(:,m),varargin{:});
end
[B,I]=sort(Ob)

x=xini(:,I(1:NP));
Ob=Ob(I(1:NP))

FES=FES+NP;
trace(1)=min(Ob);
%%%%%%%%%%%%%%%%%%%%%%%差分进化循环%%%%%%%%%%%%%%%%%%%%%
% Omega = [1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8];% sta
Omegai=mean(Xs-Xx)/10;
for ii=1 :max(NumRDP)+1
    Omegai=Omegai/5;
    Omega(ii)=Omegai;
end

for gen=1:G
    % Dynamic mutation operator
    if useDynamicParameter==1
        Lambda=exp(1-((G)/(G+1-gen)));
        F=F0*(2^Lambda);
        CR=0.5*(1+rand);
    else
        F=0.5;                                %变异算子F=[0,2]  F=0.5;
        CR=0.1;                               %交叉算子CR=[0,1] CR=0.1;
    end
    %%%%%%%%%%%%%%%%%%%%%%变异操作%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%r1,r2,r3和m互不相同%%%%%%%%%%%%%%%
    for m=1:NP
        r1=randi([1,NP],1,1);%[1,NP]之间的一个1行1列的随机整数
        while (r1==m)
            r1=randi([1,NP],1,1);
        end
        r2=randi([1,NP],1,1);
        while (r2==m)|(r2==r1)
            r2=randi([1,NP],1,1);
        end
        r3=randi([1,NP],1,1);
        while (r3==m)|(r3==r1)|(r3==r2)
            r3=randi([1,NP],1,1);
        end
        v(:,m)=x(:,r1)+F*(x(:,r2)-x(:,r3));
    end
    %%%%%%%%%%%%%%%%%%%%%%交叉操作%%%%%%%%%%%%%%%%%%%%%%%
    r=randi([1,D],1,1);%[1,D]之间的一个1行1列的随机整数
    for n=1:D%遍历每行
        cr=rand(1);
        if (cr<=CR)|(n==r)
            u(n,:)=v(n,:);
        else
            u(n,:)=x(n,:);
        end
    end
    %%%%%%%%%%%%%%%%%%%边界条件的处理%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%边界吸收%%%%%%%%%%%%%%%%%%%%%%%%%
    for n=1:D
        for m=1:NP
            if u(n,m)<Xx(n)
                if rand<0.2
                    u(n,m)=Xx(n);
                else
                    u(n,m)=Xx(n)+rand*(Xs(n)-Xx(n));
                end
            end
            if u(n,m)>Xs(n)
                if rand<0.2
                    u(n,m)=Xs(n);
                else
                    u(n,m)=Xx(n)+rand*(Xs(n)-Xx(n));
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%选择操作%%%%%%%%%%%%%%%%%%%%%%%
    for m=1:NP
        %保留小数位数
        for num=1:D
            u(num,m)=round(u(num,m),NumRDP(num));
        end
        Ob1(m)=feval(fhd,u(:,m),varargin{:});
    end
    
    FES=FES+NP;
    
    for m=1:NP
        if Ob1(m)<Ob(m)
            x(:,m)=u(:,m);
            Ob(m)=Ob1(m);%有了这一句，就不用再对所有的x求一遍误差函数了
        end
    end
    
    %     for m=1:NP
    %         Ob(m)=feval(fhd,x(:,m),varargin{:});
    %     end
    
    FES=FES+NP;
    
    [minOb,minObIndex]=min(Ob);
    
    % sta Local search
    %     useLocalSTA=1;% 1:is use sta local search; 0:NOT use sta local search
    if (useLocalSTA==1)
        SE=30;
        Range = [Xx;Xs];%
        oldbestcost=Ob(minObIndex);
        [stax,staxcost] = rotate_w(fhd,x(:,minObIndex)',SE,Range,Omega,varargin{:});
        if oldbestcost>=staxcost
            x(:,minObIndex)  =stax';
            Ob(minObIndex)  =staxcost;
        else
            disp('error?');
        end
    end
    
    
    if gen>1
        if norm(trace(gen-1)-trace(gen)) < jingdu % can be changed
            counter = counter + 1;
            if counter > 0.02*G % can be changed
                %给新个体
                [B,I]=sort(Ob);
                Ob=Ob(I);
                x=x(:,I);
                
                xini=rand(D,tempNP).*(Xs-Xx)'+Xx';
                for m=1:tempNP
                    Obxini(m)=feval(fhd,xini(:,m),varargin{:});
                end
                [B,I]=sort(Obxini);
                Ob((round(NP/2)+1):NP)=Obxini(1:round(NP/2));
                x(:,(round(NP/2)+1):NP)=xini(:,1:round(NP/2));
                counter=0;
                %                 break;
            end
        else
            counter = 0;
        end
    end
    
    %     if FES>Max_FES
    %         break;
    %     end
    if GFES>Max_FES
        break;
    end
    
    [minOb,minObIndex]=min(Ob);
    
    trace(gen+1)=min(Ob);
    % output result
    iter=gen;
    gbest=x(:,minObIndex)';
    gbestval=min(Ob);
    iterGbestval=[iterGbestval ;[iter,gbestval]];
    BestSol.Position=gbest;
    BestSol.Cost=gbestval;
    %     if iter==1
    disp(['iter=' num2str(iter) ', position=' num2str(gbest) ', fitness=' num2str(gbestval)]);
    %     else
    %         if gbestval<=iterGbestval(iter-1,end)
    %             disp(['iter=' num2str(iter) ', position=' num2str(gbest) ', fitness=' num2str(gbestval)]);
    %         end
    %     end
    
    
end
[SortOb,Index]=sort(Ob);
x=x(:,Index);
X=x(:,1);                              %最优变量
Y=min(Ob);                             %最优值

gbest=X;
gbestval=Y;
% profile viewer
disp(['iter=' num2str(iter) ', gbest_position=' num2str(gbest') ', gbest_fitness=' num2str(gbestval)]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%画图%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(trace);
% xlabel('迭代次数')
% ylabel('目标函数值')
% title('DE目标函数曲线')



%-----------------------rotate_w------------------------
function [Best,fBest] = rotate_w(funfcn,Best,SE,Range,Omega,varargin)
[Best,alpha] = update_alpha(funfcn,Best,SE,Range,Omega,varargin{:});
for i = 1:10
    [Best,fBest] = rotate(funfcn,Best,SE,Range,alpha,varargin{:});
end
end

%-----------------------update_alpha------------------------
function [Best,alpha] = update_alpha(funfcn,Best,SE,Range,Omega,varargin)

Pop_Lb = repmat(Range(1,:),SE,1);
Pop_Ub = repmat(Range(2,:),SE,1);
m  = length(Omega);

%fBest  = feval(funfcn,Best,varargin{:});
[Best,fBest] = fitness(funfcn,Best,varargin{:});
alpha = 1;
Best0 = Best;
for i = 1:m
    State = op_rotate(Best0,SE,Omega(i)); %
    %Apply  for State > Pop_Ub or State < Pop_Lb
    changeRows = State > Pop_Ub;
    State(find(changeRows)) = Pop_Ub(find(changeRows));
    changeRows = State < Pop_Lb;
    State(find(changeRows)) = Pop_Lb(find(changeRows));
    %Apply  for State > Pop_Ub or State < Pop_Lb
    [tempBest,tempfBest] = fitness(funfcn,State,varargin{:});
    if tempfBest < fBest
        Best = tempBest;
        fBest = tempfBest;
        alpha = Omega(i);
    end
end
end

%%-----------------------rotate------------------------
function [Best,fBest,flag] = rotate(funfcn,oldBest,SE,Range,alpha,varargin)

beta = 1;
Pop_Lb = repmat(Range(1,:),SE,1);
Pop_Ub = repmat(Range(2,:),SE,1);

Best = oldBest;

%fBest  = feval(funfcn,Best,varargin{:});
[Best,fBest] = fitness(funfcn,Best,varargin{:});

State = op_rotate(Best,SE,alpha); %
%Apply  for State > Pop_Ub or State < Pop_Lb
changeRows = State > Pop_Ub;
State(find(changeRows)) = Pop_Ub(find(changeRows));
changeRows = State < Pop_Lb;
State(find(changeRows)) = Pop_Lb(find(changeRows));
%Apply  for State > Pop_Ub or State < Pop_Lb
[newBest,fGBest] = fitness(funfcn,State,varargin{:});
if fGBest < fBest
    fBest = fGBest;
    Best = newBest;
    flag = 1;
else
    flag = 0;
end

if flag ==1
    State = op_translate(oldBest,Best,SE,beta);%
    %Apply  for State > Pop_Ub or State < Pop_Lb
    changeRows = State > Pop_Ub;
    State(find(changeRows)) = Pop_Ub(find(changeRows));
    changeRows = State < Pop_Lb;
    State(find(changeRows)) = Pop_Lb(find(changeRows));
    %Apply  for State > Pop_Ub or State < Pop_Lb
    [newBest,fGBest] = fitness(funfcn,State,varargin{:});
    if fGBest < fBest
        fBest = fGBest;
        Best = newBest;
    end
end
end

function y = op_rotate(Best,SE,alpha)
% 屿?¨ 旋转
% n = length(Best);
% R1 = 2*rand(SE,n)-1;
% R2 = 2*rand(SE,1)-1;
% y = repmat(Best,SE,1) + alpha*repmat(R2,1,n).*R1./repmat(sqrt(sum(R1.*R1,2)),1,n);

% %???????????ú????·??????????Щ???????÷?ē? ??
n = length(Best);
State=zeros(SE,n);
d=1/n; % d=1/n ?ù?·????d=1 ???·??
for ii=1:SE
    Rr1= 2*rand(1,n)-1;
    %     Rr2= 2*rand(1,n)-1;
    %     Rr2= 2*rand-1;% ?
    Rr2= rand^(d);%
    State(ii,:)=Best+alpha*Rr2.*Rr1/norm(Rr1);
end
y=State;

% %分解
% n = length(Best);
% R1 = 2*rand(SE,n)-1;%SE行n列的矩阵
% R2 = 2*rand(SE,1)-1;
% A=repmat(Best,SE,1);%形状和R1相同
% B=repmat(R2,1,n);%形状和R1也相吿
% C=sqrt(sum(R1.*R1,2));%R1和R1每个元素之积，形状和R1相同，然后对行元素求和?????形状和R2相同＿sum(A,2) is a column vector containing the sum of each row.
% y= A+ alpha*B.*R1./repmat(C,1,n);%

% %改写上面的超球体代码为循环形式，用于改编为有些不支持矩阵操作的语訿
% n = length(Best);
% State=zeros(SE,n);
% R1=zeros(n);
% for ii=1:SE
%     F2=0.0;
%     for dd=1:n
%         R1(dd)=2*rand-1;
%         F2=F2+ R1(dd)^2;
%     end
%     F2=sqrt(F2)+eps;
%     C=F2;
%     R2=2*rand-1;
%     for dd=1:n
%         State(ii,dd)=Best(dd)+alpha*R2*R1(dd)/C;%注意这儿丿??是一绿
%     end
% end
% y=State;

% %论文中公弿超立方体
% n = length(Best);
% State=zeros(SE,n);
% F2=0.0;
% for dd=1:n
%     F2=F2+ Best(dd)^2;
% end
% F2=sqrt(F2)+eps;
% for ii=1:SE
%     for dd=1:n
%         R1=2*rand-1;
%         State(ii,dd)=Best(dd)+alpha*R1*Best(dd)/(n*F2);%注意这儿丿??是一绿
%     end
% end
% y=State;

end

%-----------------------op_translate_zbl------------------------
function y = op_translate(oldBest,newBest,SE,beta)
% 屿?¨搜索 丿??线
%原始代码
n = length(oldBest);
y = repmat(newBest',1,SE) + beta/(norm(newBest-oldBest)+eps)*reshape(kron(rand(SE,1),(newBest-oldBest)'),n,SE);
y = y';

% A=[1 2;3 4]
% B=[4 3 ;2 1]
% C=kron(A,B)
%分解
% n = length(oldBest);
% A=repmat(newBest',1,SE);%n*SE
% B=norm(newBest-oldBest)+eps;
% C=kron(rand(SE,1),(newBest-oldBest)');
% y = A + beta/B*reshape(C,n,SE);


% %改写上面为循环形式，用于改编为有些不支持矩阵操作的语訿
% n = length(oldBest);
% State=zeros(SE,n);
% F2=0.0;
% for dd=1:n
%     F2=F2+ (newBest(dd)- oldBest(dd))^2;
% end
% F2=sqrt(F2)+eps;
% B=F2;
% for ii=1:SE
%     R1=rand;
%     for dd=1:n
%         State(ii,dd)=newBest(dd)+beta/B*R1*(newBest(dd)-oldBest(dd));%注意这儿丿??是一绿
%     end
% end
% y=State;

end


function [Best,fBest] = fitness(funfcn,State,varargin)
% calculate fitness

SE = size(State,1);
fState = zeros(SE,1);
% for i = 1:SE
%     fState(i) = feval(funfcn,State(i,:),varargin{:});
% end
fState = feval(funfcn,State',varargin{:});
[fGBest, g] = min(fState);
fBest = fGBest;
Best = State(g,:);
end






