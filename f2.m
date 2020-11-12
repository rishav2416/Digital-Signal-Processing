%% %unnormalised
m=54;                                               
q=floor(0.1*m);
r=m-10*q;
BL=25+1.9*q+4.1*r;
BH=BL+20;
stopband=[BL,BH];
passband=[BL-4,BH+4];
D1=0.384; D2=43.44;

%% %normalised

wp=passband*2*pi/260;                               
ws=stopband*2*pi/260;

%% %analog

Op=tan(wp/2);                                       
Os=tan(ws/2);

%% %frequency transform

O0=sqrt(Op(1)*Op(2));                               
B=Op(2)-Op(1);
Olp1=(B*Op)./(O0.*O0-Op.*Op);
Olp=1;
Ols1=(B*Os)./(O0.*O0-Os.*Os);
Ols=(min(abs(Ols1)));

%% %analog lowpass 

N=ceil(acosh(sqrt(D2/D1))/acosh(Ols));              
poles_lp=zeros(1,2*N);
for k=0:2*N-1
    R=sin((2*k+1)*pi/(2*N))*sinh(asinh(1/sqrt(D1))/N);
    I=cos((2*k+1)*pi/(2*N))*cosh(asinh(1/sqrt(D1))/N);
    poles_lp(k+1)=R+1i*I;
end
poles_lp=sort(poles_lp,'ComparisonMethod','real');

%% %analog bandstop

poles_bs=zeros(1,4*N);                              
i=1;
for k=1:2*N
    poles_bs(i)=(B/poles_lp(k)+sqrt((B/poles_lp(k))^2-4*O0^2))/2;
    poles_bs(i+1)=(B/poles_lp(k)-sqrt((B/poles_lp(k))^2-4*O0^2))/2;
    i=i+2;
end
zeros_bs=[-1i*O0,1i*O0]; %order N each

%% %discrete bandstop

poles=zeros(1,2*N);                                
for k=1:2*N
    poles(k)=(1+poles_bs(k))/(1-poles_bs(k));
end
zeros_disc=[(1-1i*O0)/(1+1i*O0),(1+1i*O0)/(1-1i*O0)]; %order N each

%% Paramters for BandStop to LPF Transformation
Omega_o = O0;
% B = (omega_p2-omega_p1);

%% Creating the Transfer function for the Analog Low Pass Filter
k = (poles_lp(1)*poles_lp(2)*poles_lp(3)*poles_lp(4))/sqrt(1+D1);
[numerator, denominator] = zp2tf([], poles_lp(1,1:N), k); % Transfer Function is multiplied with a constant to make DC gain=1 at s=0

%% Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(numerator,s)/poly2sym(denominator,s);
analog_bsf(s) = analog_lpf((B*s)/((s*s)+(Omega_o*Omega_o)));
discrete_bsf(z) = analog_bsf((z-1)/(z+1));

%% Coefficients of Analog LPF
[nums, dens] = numden(analog_lpf(s));
num_lpf = sym2poly(expand(nums));
den_lpf = sym2poly(expand(dens));
num_lpf = num_lpf./den_lpf(1);
den_lpf = den_lpf./den_lpf(1);
[H_lpf,f_lpf] = freqs(num_lpf, den_lpf, 10000);
pos_stop = 0;
pos_start = 0;
for i=1:length(f_lpf)
    if(f_lpf(i)>=Ols)
        pos_stop = i;
        break;
    end
end
for i=1:length(f_lpf)
    if(f_lpf(i)>=1)
        pos_start = i;
        break;
    end
end
figure
hold on;
plot(f_lpf, abs(H_lpf)); 
title('H_{analog, LPF}(s_L) on normalized frequency axis'); 
xlim([0.1,3.5]); ylim([0, 1.2]); xlabel('s_L'); ylabel('|H_{analog,LPF}(s_L)|'); 
plot(f_lpf(pos_stop),abs(H_lpf(pos_stop)),'r*');
plot(f_lpf(pos_start),abs(H_lpf(pos_start)),'r*');
grid on;

%% Coefficients of Analog BSF
[nums_b, dens_b] = numden(analog_bsf(s));
num_bsf = sym2poly(nums_b);
den_bsf = sym2poly(dens_b);
num_bsf = num_bsf./den_bsf(1);
den_bsf = den_bsf./den_bsf(1);
[H_bsf,f_bsf] = freqs(num_bsf, den_bsf, 10000);
critical_points_bsf = zeros(5);
check_points_bsf = [Op(1), Os(1), O0, Os(2), Op(2)];
for i=1:5
    for l=1:length(f_bsf)
        if(f_bsf(l)>=check_points_bsf(i))
            critical_points_bsf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_bsf, abs(H_bsf)); 
title('H_{analog, BSF}(s) on normalized frequency axis'); 
xlim([0.1,3.5]); ylim([0, 1.2]); xlabel('s'); ylabel('|H_{analog,BSF}(s)|'); 
for i=1:5
    plot(f_bsf(critical_points_bsf(i)), abs(H_bsf(critical_points_bsf(i))),'r*');
end
grid on;
%% Discrete Filter Coefficients

[nums_b2, dens_b2] = numden(discrete_bsf(z));
num_bsf2 = sym2poly(nums_b2);
den_bsf2 = sym2poly(dens_b2);
num_bsf2 = num_bsf2./den_bsf2(1);
den_bsf2 = den_bsf2./den_bsf2(1);
[H_bsf2,f_bsf2] = freqz(num_bsf2, den_bsf2, 1024*1024, 260e3);
critical_points_bsf = zeros(4);
check_points_bsf = [stopband(1)*1000, passband(1)*1000, passband(2)*1000, stopband(2)*1000];
for i=1:length(critical_points_bsf)
    for l=1:length(f_bsf2)
        if(f_bsf2(l)>=check_points_bsf(i))
            critical_points_bsf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_bsf2, abs(H_bsf2)); 
title('H_{discrete, BSF}(z) on un-normalized frequency axis'); 
ylim([0, 1.2]); xlabel('Frequency in Hz'); ylabel('|H_{discrete,BSF}(z)|'); 
for i=1:length(critical_points_bsf)
    plot(f_bsf2(critical_points_bsf(i)), abs(H_bsf2(critical_points_bsf(i))),'r*');
end
grid on;
%%
fvtool(num_bsf2, den_bsf2);




%%
f_samp = 260e3;
delta = 0.15; %stopband tolerance
A = -20*log10(delta);

%% Shape parameter
if(A<21)
    alpha=0;
elseif((A>=21) && (A<=50))
    alpha = 0.5842*((A-21)^0.4)+0.07886*(A-21);
else
    alpha = 0.1102*(A-8.7);
end
%% Length of FIR filter
delta_omega = 0.03076*pi;
Nk = ceil((A-8)/(2*2.285*delta_omega));

%% Kaiser Window formation
N_corrected = (Nk+7); %N is corrected due to poor lower bound derived previously
n = (2*(N_corrected)+1);
beta = alpha/Nk;
kaiser_window = kaiser(n,beta);

%% FIR BPF Filter formation
passband_1 = wp(1)+0.0154*pi;
passband_2 = wp(2)-0.0154*pi;
bsf = lpf_FIR(pi,n) - lpf_FIR(passband_2,n) + lpf_FIR(passband_1,n);
bsf_FIR = bsf.*kaiser_window';
fvtool(bsf_FIR);
%% Plotting FIR filter on unnormalized frequencies
[H_FIR, f_FIR] = freqz(bsf_FIR,1,1024,f_samp);
critical_points_bsf = zeros(4);
check_points_bsf = [stopband(1)*1000, passband(1)*1000, passband(2)*1000, stopband(2)*1000];
for i=1:length(critical_points_bsf)
    for l=1:length(f_FIR)
        if(f_FIR(l)>=check_points_bsf(i))
            critical_points_bsf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_FIR, abs(H_FIR)); 
title('H_{discrete, BPF}(z) on un-normalized frequency axis'); 
ylim([0, 1.2]); xlabel('Frequency in Hz'); ylabel('|H_{discrete,BPF}(z)|'); 
for i=1:length(critical_points_bsf)
    plot(f_FIR(critical_points_bsf(i)), abs(H_FIR(critical_points_bsf(i))),'r*');
end
grid on;

