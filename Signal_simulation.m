%% �����ź�


function X=Signal_simulation(IF,PW,PRT,SNR,Noise_Power,Fs,Time,n,Siganl_Tag,Bandwidth)%,Array_Tag,kelmy,kelmz,theta,fe,gamma,eita)
%IF:��Ƶ(MHZ),IF��ֵ�����100
%PW:����(us)
%PRT:�����ظ�����(us)
%SNR:�����(dB)
%Noise_Power:��������(W)
%Fs:����Ƶ��(HZ)
%Time:����ʱ��(ms)
%n:������(��)
%Signal_Tag:�źű�ǩ(1.����Ƶ�źţ�2.Chirp�źţ�3.NS+Chirp�ź�)
%Bandwidth:����(MHz)
%Array_Tag:���б�ǩ(1.L��2.Բ��3.����)
%kelmy:Y�ᣬ������ΪԲ����Ϊ����Ԫ��(��)
%,kelmz:Z����Ԫ������ΪԲ����ΪԲ��뾶
%theta,fe:��λ�ǣ�������(��)
%gamma,eita:���������ǣ�������λ��(��)
%-----------------------------------�����ź�-----------------------------------%


twpi = 2*pi;
rad = pi/180;  %�Ƕ��뻡�ȵ�ת��
deg = 180/pi;  %������Ƕȵ�ת��

PW=PW*1e-6;
PRT=PRT*1e-6;
Time=Time*1e-3;
Signal_Power=Noise_Power*10^(SNR/10);%�źŹ���
Signal_Amplitude=sqrt(2*Signal_Power);%�źŷ���
Center_frequency=IF*1e6;%��Ƶ
Ts=1/Fs;%�������
Sample_frequency=Fs;%����Ƶ��
Sample_number_PW=round(PW/Ts);%�����ڵĲ�������
Sample_number_zeros=round((PRT-PW)/(Ts));%�������
Sample_number_PRT=Sample_number_zeros+Sample_number_PW;%һ�������ظ������еĲ�������
Sample_index=[0:Sample_number_PW-1];%�����±�
t_base=linspace(0,PW,Sample_number_PW);%һ��PW�еĲ�������
t_base0=linspace(0,PRT-PW,Sample_number_zeros);%һ��PRT�����в���ĵ�
Signal_zeros=0*t_base0;%����Ҫ����ĵ㻯Ϊ0

Signal_Sample_number=round(Time*Fs);
Pulse_number=ceil(Time/PRT);
Signal=[];
% for i_s=1:Signal_number
%����Ƶ�ź�
if Siganl_Tag==1
    Array_frequency_n=Center_frequency+rand*1e9;%ʹ����Ƶ������ƫ���Է�ֹ���
    for i=1:Pulse_number%��i_s���źŵĵ�i��PW
        Signal_n_PW{1,i}=exp(1i*2*pi*Array_frequency_n/Sample_frequency*Sample_index);
    end
    for i=1:Pulse_number%��i_s���źŵĵ�i��PRT
        Signal_n_PRT{1,i}=[Signal_n_PW{1,i},Signal_zeros];
    end
    Signal_n=zeros(1,Pulse_number*Sample_number_PRT);
    for i=1:Pulse_number%��i_s���źŵ�Pulse_number��PRT
        Signal_n(1,((i-1)*Sample_number_PRT+1):(i*Sample_number_PRT))=Signal_n_PRT{1,i};
    end
    Signal=[Signal;Signal_n];
    %disp('ѡ����ź�Ϊ����Ƶ�ź�');
    
%Chirp�ź�
elseif Siganl_Tag==2
    Bandwidth=Bandwidth*1e6;%����
    K=Bandwidth/PW;%б��
    
    Array_frequency_n=Center_frequency+rand*1e7;
    for i=1:Pulse_number%��i_s���ź��е�i��PW
        Signal_n_PW{1,i}=exp(1i*2*pi*(Array_frequency_n*t_base+0.5*K*t_base.^2));
    end
    for i=1:Pulse_number%��i_s���ź��е�i��PRT
        Signal_n_PRT{1,i}=[Signal_n_PW{1,i},Signal_zeros];  
    end
    for i=1:Pulse_number%��i_s���źŵ�Pulse_number��PRTд��ʱ��ʸ��
        Signal_n(1,((i-1)*Sample_number_PRT+1):(i*Sample_number_PRT))=Signal_n_PRT{1,i};
    end
    Signal=[Signal;Signal_n];
    %disp('ѡ����ź�ΪChirp�ź�');
else
    disp('����Signal_Tag����');
end
% end

Signal=Signal_Amplitude*Signal;
Signal=Signal(1:Signal_Sample_number);
X=Signal(1:n);

%% -----------------------------------���ź�ʱ����Ƶ��ͼ-----------------------------------%
% Sample=real(Signal(1:Sample_number_PW));
% %Sample=imag(Signal(1:Sample_number_PW));
% %FFT_Signal=abs(fft(Signal(1:2*Sample_number_PW-1)));
% FFT_Signal=abs(fft(Signal(1:Sample_number_PW)));
% figure(1);
% plot(Sample);%ʵ��ʱ��ͼ
% xlabel('ʱ���� t /s');
% ylabel('�ź�S');
% title('�ź� (ʱ��)');
% figure(2);
% plot(0:Sample_frequency/Sample_number_PW:(Sample_frequency-Sample_frequency/Sample_number_PW),FFT_Signal);%Ƶ��ͼ
% xlabel('Ƶ���� f /Hz');
% ylabel('���� P');
% title('�ź� (Ƶ��)');



% %-----------------------------------��������-----------------------------------%
% c=3e8;%����
% lamda=c/Center_frequency;%����
% dd=lamda/2;%��Ԫ���
% A=[];%��������
% theta=theta*rad;%��λ��
% fe=fe*rad;%������
% gamma=gamma*rad;%����������
% eita=eita*rad;%������λ��   
% n=length(Signal);%�źų���
% 
% 
% %L��
% if Array_Tag==1
%     kelm=kelmy+kelmz;%��Ԫ����
%     dy=0:dd:(kelmy-1)*dd;%Y����Ԫ�ֲ�
%     dz=dd:dd:kelmz*dd;%Z����Ԫ�ֲ�
% 
%     for i=1:Signal_number
%         y=exp(-1i*twpi/lamda*dy.'*sin(theta(1,i))*cos(fe(1,i)));
%         z=exp(-1i*twpi/lamda*dz.'*sin(fe(1,i)));
%         As=[y;z];%������ʸ��
%         polar=[cos(theta(1,i))*cos(gamma(1,i))+sin(fe(1,i))*sin(theta(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i));-cos(fe(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i))];%��������ʸ��
%         Steering_vector=kron(As,polar);%�������ڿ˻�
%         A=[A,Steering_vector];
%     end
%  Noise= sqrt(Noise_Power)*( sqrt(0.5)*randn(2*kelm,n) +sqrt(0.5)*1i*randn(2*kelm,n));
% disp('ѡ�������ΪL��');
% %Բ��
% elseif Array_Tag==2
%    for i=1:Signal_number
%        as=[];%������ʸ��
%        R=kelmz;%Բ��뾶
%        for i_k=1:kelmy
%             m=i_k-1;
%             a_kelm=exp(-1i*R*twpi*(cos(twpi*m/kelm)*cos(fe(i))*sin(theta(i))+sin(twpi*m/kelm)*sin(fe(i))));
%             as=[as;a_kelm];
%        end
%         polar=[cos(theta(1,i))*cos(gamma(1,i))+sin(fe(1,i))*sin(theta(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i));-cos(fe(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i))];%x���򼫻�-�Ƕ�����ʸ��
%         Steering_vector=kron(as,polar);
%         A=[A,Steering_vector];
%    end
%  Noise= sqrt(Noise_Power)*( sqrt(0.5)*randn(2*kelmy,n) +sqrt(0.5)*1i*randn(2*kelmy,n));
% disp('ѡ�������ΪԲ��');
% %����
% elseif Array_Tag==3
%     kelm=kelmy*kelmz;%��Ԫ����             
%     dy=0:dd:(kelmy-1)*dd;%y����Ԫ��������
%     dz=0:dd:(kelmz-1)*dd;%z��
%     for i=1:Signal_number
%     as_y=exp(-1i*twpi*dy.'*sin(theta(1,i))*cos(fe(1,i)));%y�ᵼ��ʸ��
%     as_z=exp(-1i*twpi*dz.'*sin(fe(1,i)));%z�ᵼ��ʸ��
%     as=[as_y;as_z];%��������źŵ���ʸ��
%     polar=[cos(theta(1,i))*cos(gamma(1,i))+sin(fe(1,i))*sin(theta(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i));-cos(fe(1,i))*sin(gamma(1,i))*exp(1i*eita(1,i))];%����-�Ƕ�����ʸ��
%     Steering_vector=kron(as,polar);%������ʸ��
%     A=[A,Steering_vector];
%     end
%      Noise= sqrt(Noise_Power)*( sqrt(0.5)*randn(2*kelm,n) +sqrt(0.5)*1i*randn(2*kelm,n));
% disp('ѡ�������Ϊ����');
% else disp('����Signal_Tag����');
% end
% 
% 
% %-----------------------------------�������н����ź�-----------------------------------%
% X=A*Signal+Noise;
% end
