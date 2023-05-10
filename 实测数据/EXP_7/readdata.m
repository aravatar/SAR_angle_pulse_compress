% close all;
clear all
clc

fid = fopen('DATA_1.BIN','rb');
[adc_data]=fread(fid);
fclose(fid);

% fid = fopen('TEST/exp_4/DATA_1.BIN','rb');
% [adc_data_16]=fread(fid,'uint16');
% fclose(fid);

%%

azi_num = length(adc_data)/(2080*2+32);
adc_data = adc_data(1:(2080*2+32)*azi_num);
% data_d = zeros(length(adc_data)/2);
% data = zeros(length(adc_data)/2);
[adc_data] = reshape(adc_data,(2080*2+32),length(adc_data)/(2080*2+32));
Gps_row_data = adc_data(4161:4192,:);

adc_data = adc_data(1:4160,:);
[adc_data] = reshape(adc_data,1,size(adc_data,1)*size(adc_data,2));

[m,n] = size(Gps_row_data);
Gps_data = zeros(8, n);

for i = 1:4
   Gps_data(i,:) =  Gps_row_data((i-1)*4+4, :)*2^24 + Gps_row_data((i-1)*4+3, :)*2^16 + Gps_row_data((i-1)*4+2, :)*2^8 + Gps_row_data((i-1)*4+1, :);
   if Gps_data(i,:) > 2^31
       Gps_data(i,:) = Gps_data(i,:) - 2^32;
   end
end

for i = 5:8
   Gps_data(i,:) =  Gps_row_data((i-5)*2+2 +16, :)*2^8 +  Gps_row_data((i-5)*2+1 +16, :);
   if i<8
       if Gps_data(i,:) > 2^15
           Gps_data(i,:) = Gps_data(i,:) - 2^16;
       end
   end
end

for i=1:(length(adc_data)/2)
    data_d(i) = adc_data(2*i-1)*256 + adc_data(2*i);
end

for i=1:length(data_d) 
    if(data_d(i)>32767)
        data(i) = data_d(i) - 65536;
    else
        data(i) = data_d(i);
    end
end
FrmDtLen = 2048;
FrmHdLen = 32;
FrmLen = FrmDtLen + FrmHdLen;
[data] = reshape(data,FrmLen,length(data)/FrmLen);
data = data(FrmHdLen+1:FrmLen,:);
[data] = reshape(data,1,FrmDtLen*size(data,2));
%%%%实部虚部分离
data_real(length(data)/2)=0;
data_image(length(data)/2)=0;
for i=1:(length(data)/4)
    data_real(2*i-1)=data(4*i-3);
    data_real(2*i)=data(4*i-2);
    
    data_image(2*i-1)=data(4*i-1);
    data_image(2*i)=data(4*i);
end
%%%AOA
% lambda=3.90e-3;
frames=length(data)/FrmDtLen;
channel=1;
ch_num=1;
N=1024;
range=zeros(frames,N);
range_fft=zeros(frames,N);
range_x=0.05:0.05:51.2;
for i=1:frames
    range(i,:)=data_real(ch_num*N*(i-1)+channel:ch_num*N*(i-1)+channel+N-1)+1i*data_image(ch_num*N*(i-1)+channel:ch_num*N*(i-1)+channel+N-1);
    win = hamming(N);
    range(i,:) = range(i,:).*win.';%距离向加窗
    range_fft(i,:) = fft(range(i,:).*(range_x).^4);
end

figure;
imagesc(abs(range_fft));

%% 参数
C=3e8;%光速
fc=77e9;%雷达载频
lambda=C/fc;%波长
V=1;%雷达移动速度
% La = 0.886*2*V/3e9;%合成孔径长度
% Ta=La/V;%合成孔径时间
H=4;%高度
Y=6;%雷达作用距离
Xmin=0;%目标区域方位向范围
Xmax=10;%[Xmin,Xmax]
PRF=2000;%脉冲重复频率
Nslow=PRF*(Xmax-Xmin)/V;%慢时间采样数
Br=3e9;%距离向带宽
Tr=52e-6;%发射信号时宽
K=Br/Tr;%调频斜率
pr=C/(2*K*Tr);%%距离分辨率（方位~）
La=2*pr;
Ta=La/V;%合成孔径时间
Fr=5e6;%距离向采样频率
rmin=sqrt(H^2+Y^2);

% bmax=floor(PRF/3);%多普勒最大带宽
% Ka = 2*V^2/(rmin*lambda);   %多普勒调频率
% Ta=floor(bmax/Ka);
% La=Ta*V;

% %% 距离向匹配滤波
% rmax=sqrt(H^2+Y^2+(La/2)^2);
% t=linspace(2*rmin/C-Tr/2,2*rmax/C+Tr/2,N);%线性调频信号范围
% for i = 1:frames
%     range_fft(i,:) = range_fft(i,:).*(exp(-1j*2*pi*fc*t));   %回波去载频
% end
% tt=linspace(0,Tr,N);
% hk=exp(1j*pi*K*tt.^2);    %距离向匹配滤波器
% juya = zeros(size(range_fft));  %距离向压缩之后
% for i=1:frames
%     juya(i,:) = ifft(fft(range_fft(i,:)).*conj(fft(hk,N)));
% end
% range_fft=juya;
% figure(2);
% imagesc(abs(range_fft));
% title('距离向压缩'),xlabel('距离向'),ylabel('方位向');
%% 方位向傅里叶变换
for k=1:N
    range_fft(:,k) = fftshift(fft(range_fft(:,k)));
end
figure(3);
imagesc(abs(range_fft));
title('方位向傅里叶变换'),xlabel('距离向'),ylabel('方位向');
%% 距离徙动校正
RCMC = zeros(size(range_fft));
f =linspace(-PRF/2,PRF/2,frames);
h = waitbar(0,'Sinc插值中......');  %生成一个进度条
for k = 1:frames
    deltaR = ceil((lambda*f(k)/V)^2*rmin*8/8);  %注意这个f，到底是哪个f
    range_fft1(1,:)=[range_fft(k,:),zeros(1,deltaR)];
    for n = 1:N  %快时间
        RCMC(k,n) = range_fft1(1,n+deltaR);
    end
    clear range_fft1;
    waitbar(k/frames);
end
close(h);  %关闭进度条

%RD算法,sinc插值
% RCMC = zeros(size(range_fft));
% N1 = 6;  %插值点数
% Rp = sqrt(sum(([Y,0]-[0,H]).^2));  %目标到雷达的最近距离
% h = waitbar(0,'Sinc插值中......');  %生成一个进度条
% for k = 1:frames
%     f =linspace(-PRF/2,PRF/2,frames);
%     deltaR = (lambda*f(k)/V)^2*rmin/8;  %注意这个f，到底是哪个f
%     DU = 2*deltaR*Fr/C;
%     du = DU-floor(DU);
%     %     kernel_norm = sum(sinc(du-(-N/2:N/2-1)));
%     for n = N1/2+1:N  %快时间
%         for m = -N1/2:N1/2-1
%             if n+floor(DU)+m>N
%                 RCMC(k,n) = RCMC(k,n)+range_fft(k,N)*sinc(DU-m);
%             else
%                 RCMC(k,n) = RCMC(k,n)+range_fft(k,n+floor(DU)+m)*sinc(du-m);
%             end
%         end
%     end
%     waitbar(k/frames);
% end
% close(h);  %关闭进度条

% 图形四
figure(4)
imagesc(abs(RCMC));
title('RCMC后的信号幅度'),xlabel('距离向'),ylabel('方位向');

%% 方位压缩
out=zeros(frames,N);
% ttt=linspace(0,Ta,Ta*PRF);%方位向照射到目标的短时间范围
ttt=linspace(-Ta/2,Ta/2,256);%方位向照射到目标的短时间范围
win2=hamming(length(ttt))';
ha = win2 .* exp(-1j*pi*2*(V^2)/(lambda*rmin)*(ttt).^2);%方位向匹配滤波器
for ii = 1:N
    out(:,ii) = (ifft(fftshift(RCMC(:,ii)).*fft(ha,frames)'));%方位向压缩
end
% 图形五
figure(5);
imagesc(abs(out));
% mesh(abs(out));
title('方位压缩'),xlabel('距离向'),ylabel('方位向');
