%% 生成信号


function X=Signal_simulation(IF,PW,PRT,SNR,Noise_Power,Fs,Time,n,Siganl_Tag,Bandwidth)%,Array_Tag,kelmy,kelmz,theta,fe,gamma,eita)
%IF:中频(MHZ),IF数值需高于100
%PW:脉宽(us)
%PRT:脉冲重复周期(us)
%SNR:信噪比(dB)
%Noise_Power:噪声功率(W)
%Fs:采样频率(HZ)
%Time:采样时长(ms)
%n:快拍数(个)
%Signal_Tag:信号标签(1.单载频信号，2.Chirp信号，3.NS+Chirp信号)
%Bandwidth:带宽(MHz)
%Array_Tag:阵列标签(1.L阵，2.圆阵，3.矩阵)
%kelmy:Y轴，若阵型为圆阵则为总阵元数(个)
%,kelmz:Z轴阵元数，若为圆阵则为圆阵半径
%theta,fe:方位角，俯仰角(度)
%gamma,eita:极化辅助角，极化相位差(度)
%-----------------------------------生成信号-----------------------------------%


twpi = 2*pi;
rad = pi/180;  %角度与弧度的转换
deg = 180/pi;  %弧度与角度的转换

PW=PW*1e-6;
PRT=PRT*1e-6;
Time=Time*1e-3;
Signal_Power=Noise_Power*10^(SNR/10);%信号功率
Signal_Amplitude=sqrt(2*Signal_Power);%信号幅度
Center_frequency=IF*1e6;%中频
Ts=1/Fs;%采样间隔
Sample_frequency=Fs;%采样频率
Sample_number_PW=round(PW/Ts);%脉宽内的采样点数
Sample_number_zeros=round((PRT-PW)/(Ts));%补零点数
Sample_number_PRT=Sample_number_zeros+Sample_number_PW;%一个脉冲重复周期中的采样点数
Sample_index=[0:Sample_number_PW-1];%采样下标
t_base=linspace(0,PW,Sample_number_PW);%一个PW中的采样点数
t_base0=linspace(0,PRT-PW,Sample_number_zeros);%一个PRT中所有补零的点
Signal_zeros=0*t_base0;%将需要补零的点化为0

Signal_Sample_number=round(Time*Fs);
Pulse_number=ceil(Time/PRT);
Signal=[];
% for i_s=1:Signal_number
%单载频信号
if Siganl_Tag==1
    Array_frequency_n=Center_frequency+rand*1e9;%使中心频率有所偏移以防止相干
    for i=1:Pulse_number%第i_s个信号的第i个PW
        Signal_n_PW{1,i}=exp(1i*2*pi*Array_frequency_n/Sample_frequency*Sample_index);
    end
    for i=1:Pulse_number%第i_s个信号的第i个PRT
        Signal_n_PRT{1,i}=[Signal_n_PW{1,i},Signal_zeros];
    end
    Signal_n=zeros(1,Pulse_number*Sample_number_PRT);
    for i=1:Pulse_number%第i_s个信号的Pulse_number个PRT
        Signal_n(1,((i-1)*Sample_number_PRT+1):(i*Sample_number_PRT))=Signal_n_PRT{1,i};
    end
    Signal=[Signal;Signal_n];
    %disp('选择的信号为单载频信号');
    
%Chirp信号
elseif Siganl_Tag==2
    Bandwidth=Bandwidth*1e6;%带宽
    K=Bandwidth/PW;%斜率
    
    Array_frequency_n=Center_frequency+rand*1e7;
    for i=1:Pulse_number%第i_s个信号中第i个PW
        Signal_n_PW{1,i}=exp(1i*2*pi*(Array_frequency_n*t_base+0.5*K*t_base.^2));
    end
    for i=1:Pulse_number%第i_s个信号中第i个PRT
        Signal_n_PRT{1,i}=[Signal_n_PW{1,i},Signal_zeros];  
    end
    for i=1:Pulse_number%第i_s个信号的Pulse_number个PRT写成时域矢量
        Signal_n(1,((i-1)*Sample_number_PRT+1):(i*Sample_number_PRT))=Signal_n_PRT{1,i};
    end
    Signal=[Signal;Signal_n];
    %disp('选择的信号为Chirp信号');
else
    disp('输入Signal_Tag错误');
end
% end

Signal=Signal_Amplitude*Signal;
Signal=Signal(1:Signal_Sample_number);
X=Signal(1:n);

%% -----------------------------------画信号时域与频域图-----------------------------------%
% Sample=real(Signal(1:Sample_number_PW));
% %Sample=imag(Signal(1:Sample_number_PW));
% %FFT_Signal=abs(fft(Signal(1:2*Sample_number_PW-1)));
% FFT_Signal=abs(fft(Signal(1:Sample_number_PW)));
% figure(1);
% plot(Sample);%实部时域图
% xlabel('时间轴 t /s');
% ylabel('信号S');
% title('信号 (时域)');
% figure(2);
% plot(0:Sample_frequency/Sample_number_PW:(Sample_frequency-Sample_frequency/Sample_number_PW),FFT_Signal);%频域图
% xlabel('频率轴 f /Hz');
% ylabel('能量 P');
% title('信号 (频域)');



% %-----------------------------------生成阵列-----------------------------------%
% c=3e8;%光速
% lamda=c/Center_frequency;%波长
% dd=lamda/2;%阵元间距
% A=[];%阵列流形
% theta=theta*rad;%方位角
% fe=fe*rad;%俯仰角
% gamma=gamma*rad;%极化辅助角
% eita=eita*rad;%极化相位差   
% n=length(Signal);%信号长度
% 
% 
% %L阵
% if Array_Tag==1
%     kelm=kelmy+kelmz;%阵元总数
%     dy=0:dd:(kelmy-1)*dd;%Y轴阵元分布
%     dz=dd:dd:kelmz*dd;%Z轴阵元分布
% 
%     for i=1:Signal_number
%         y=exp(-1i*twpi/lamda*dy.'*sin(theta(1,i))*cos(fe(1,i)));
%         z=exp(-1i*twpi/lamda*dz.'*sin(fe(1,i)));
%         As=[y;z];%空域导向矢量
%         polar=[cos(theta(1,i))*cos(gamma(1,i))+sin(fe(1,i))*sin(theta(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i));-cos(fe(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i))];%极化域导向矢量
%         Steering_vector=kron(As,polar);%做克罗内克积
%         A=[A,Steering_vector];
%     end
%  Noise= sqrt(Noise_Power)*( sqrt(0.5)*randn(2*kelm,n) +sqrt(0.5)*1i*randn(2*kelm,n));
% disp('选择的阵型为L阵');
% %圆阵
% elseif Array_Tag==2
%    for i=1:Signal_number
%        as=[];%空域导向矢量
%        R=kelmz;%圆阵半径
%        for i_k=1:kelmy
%             m=i_k-1;
%             a_kelm=exp(-1i*R*twpi*(cos(twpi*m/kelm)*cos(fe(i))*sin(theta(i))+sin(twpi*m/kelm)*sin(fe(i))));
%             as=[as;a_kelm];
%        end
%         polar=[cos(theta(1,i))*cos(gamma(1,i))+sin(fe(1,i))*sin(theta(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i));-cos(fe(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i))];%x方向极化-角度域导向矢量
%         Steering_vector=kron(as,polar);
%         A=[A,Steering_vector];
%    end
%  Noise= sqrt(Noise_Power)*( sqrt(0.5)*randn(2*kelmy,n) +sqrt(0.5)*1i*randn(2*kelmy,n));
% disp('选择的阵型为圆阵');
% %矩阵
% elseif Array_Tag==3
%     kelm=kelmy*kelmz;%阵元总数             
%     dy=0:dd:(kelmy-1)*dd;%y轴阵元坐标序列
%     dz=0:dd:(kelmz-1)*dd;%z轴
%     for i=1:Signal_number
%     as_y=exp(-1i*twpi*dy.'*sin(theta(1,i))*cos(fe(1,i)));%y轴导向矢量
%     as_z=exp(-1i*twpi*dz.'*sin(fe(1,i)));%z轴导向矢量
%     as=[as_y;as_z];%面阵空域信号导向矢量
%     polar=[cos(theta(1,i))*cos(gamma(1,i))+sin(fe(1,i))*sin(theta(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i));-cos(fe(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i))];%极化-角度域导向矢量
%     Steering_vector=kron(as,polar);%面阵导向矢量
%     A=[A,Steering_vector];
%     end
%      Noise= sqrt(Noise_Power)*( sqrt(0.5)*randn(2*kelm,n) +sqrt(0.5)*1i*randn(2*kelm,n));
% disp('选择的阵型为矩阵');
% else disp('输入Signal_Tag错误');
% end
% 
% 
% %-----------------------------------生成阵列接受信号-----------------------------------%
% X=A*Signal+Noise;
% end
