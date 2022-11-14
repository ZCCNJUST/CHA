tic
clear;%  清楚内存中所有的或指定的变量
CHA1=0;

cishu=1;
step=75/256;
%K=randn(1024,1024)+1i*randn(1024,1024);%  randn(m,n)生成正态分布的随机二维数组256*256，用于生成一个随机复数，表示随机散射介质，即传输矩阵K
for xunhun=1:cishu%  清楚内存中所有的或指定的变量
n=16;
vs=1.5;
K=randn(n^2,n^2)+1i*randn(n^2,n^2);
K_shift=randn(n^2,n^2)+1i*randn(n^2,n^2);
A=zeros(256,1);%  生成全零数组1024*1的列矩阵
A(:,1)=exp(1i*0);%  将二维矩阵的第一列赋值，实际上表示一个16*16的方阵，是入射光源

% load K_shift;
% K_shift=randn(256,64)+1i*randn(256,64);
 K_new=K;
for i=1:step*256
    K_new(:,i)=K_shift(:,i);
end

%参考光强的计算
I0=K*A;
I00=abs(I0).^2;
Irefvalue=mean(I00);

%初始化：生成一个1024阶的哈达玛矩阵
% H=hadamard(256);
% H=walsh(256);
H=cchdm(256);

%设置一个初始最佳相位图,初始相位均为0
F0=zeros(256,1);
for ii=1:15
for i=1:250
%利用四步移相机制计算每一阶矩阵的相位图
f1=F0+(0*180/2)*((H(:,i)+1)/2);
f2=F0+(1*180/2)*((H(:,i)+1)/2);
f3=F0+(2*180/2)*((H(:,i)+1)/2);
f4=F0+(3*180/2)*((H(:,i)+1)/2);

%假设聚焦点121
f1=f1/180*pi;
f2=f2/180*pi;
f3=f3/180*pi;
f4=f4/180*pi;
I1=abs(K(121,:)*exp(1i*f1)).^2;
I2=abs(K(121,:)*exp(1i*f2)).^2;
I3=abs(K(121,:)*exp(1i*f3)).^2;
I4=abs(K(121,:)*exp(1i*f4)).^2;

toangle=(I1-I3)+1i*(I2-I4);
f_opt=angle(toangle)/pi*180;%幅角主值，在-pi到pi区间
F1=F0+f_opt*((H(:,i)+1)/2);
for j=1:256
    if(F1(j,1)>180)
        F1(j,1)=F1(j,1)-360;
    end
    if(F1(j,1)<-180)
        F1(j,1)=360+F1(j,1);
    end
end
F0=F1;
toPaint=abs(K*exp(1i*(F0/180*pi))).^2+vs*Irefvalue*abs(randn(n^2,1));
Av_paint(i+250*(ii-1),1)=toPaint(121,1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%初始化：生成一个1024阶的哈达玛矩阵
% H=hadamard(256);
%设置一个初始最佳相位图,初始相位均为0
%F0=zeros(256,1);%要延续扰动前的最佳相位对吧
for ii=1:5
for i=1:250
%利用四步移相机制计算每一阶矩阵的相位图
f1=F0+(0*180/2)*((H(:,i)+1)/2);
f2=F0+(1*180/2)*((H(:,i)+1)/2);
f3=F0+(2*180/2)*((H(:,i)+1)/2);
f4=F0+(3*180/2)*((H(:,i)+1)/2);

%假设聚焦点121
f1=f1/180*pi;
f2=f2/180*pi;
f3=f3/180*pi;
f4=f4/180*pi;
I1=abs(K_new(121,:)*exp(1i*f1)).^2;
I2=abs(K_new(121,:)*exp(1i*f2)).^2;
I3=abs(K_new(121,:)*exp(1i*f3)).^2;
I4=abs(K_new(121,:)*exp(1i*f4)).^2;

toangle=(I1-I3)+1i*(I2-I4);
f_opt=angle(toangle)/pi*180;%幅角主值，在-pi到pi区间
F1=F0+f_opt*((H(:,i)+1)/2);
for j=1:256
    if(F1(j,1)>180)
        F1(j,1)=F1(j,1)-360;
    end
    if(F1(j,1)<-180)
        F1(j,1)=360+F1(j,1);
    end
end
F0=F1;
toPaint=abs(K_new*exp(1i*(F0/180*pi))).^2;
Av_paint(3750+i+250*(ii-1),1)=toPaint(121,1);
end
end

% Iref=abs(K*A).^2;
% figure(1)
% imshow(reshape(Iref,32,32),[]);

F0=F0/180*pi;
Iout=abs(K_new*exp(1i*F0)).^2;
% figure(2)
% imshow(reshape(Iout,16,16),[]);

Av=Iout(121,1)/Irefvalue;
Av_paint=Av_paint/Irefvalue;
for i=1:5000
Av_paint_new((i-1)*4+1:i*4,1)=Av_paint(i,1);
end
CHA1=CHA1+Av_paint_new;
end
%%
CHA1=CHA1/cishu;
%%
figure()
x2=1:20000; %确定图表范围
plot(x2,CHA1,'LineWidth',1)
title('Nosie:1.5<I_0>','FontWeight','bold')
xlabel('Measurements','FontWeight','bold')
ylabel('Enhancement','FontWeight','bold')
set(gca, 'FontSize',10) % 设置坐标轴字体是 8
axis([0 20000 0 250])
legend('CHA');

toc
disp(['运行时间: ',num2str(toc)]);