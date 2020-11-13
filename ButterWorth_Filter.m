%%  %unnormalised

m=56;                                              
q=floor(0.1*m);
r=m-10*q;
BL=25+1.7*q+6.1*r;
BH=BL+20;
passband=[BL,BH];
stopband=[BL-4,BH+4];
D1=0.384; D2=43.44;

%% %normalised

wp=passband*2*pi/330;                               
ws=stopband*2*pi/330;

%%  %analog

Op=tan(wp/2);                                      
Os=tan(ws/2);

%%  %frequency transform

O0=sqrt(Op(1)*Op(2));                              
B=Op(2)-Op(1);
Olp1=(Op.*Op-O0.*O0)./(B*Op);
Olp=1;
Ols1=(Os.*Os-O0.*O0)./(B*Os);
Ols=(min(abs(Ols1)));

%% %analog lowpass 

N=ceil(log(sqrt(D2/D1))/log(Ols));                  
Oc=0.5*(Olp/(D1)^(1/(2*N))+Ols/(D2)^(1/(2*N)));
poles_lp=zeros(1,2*N);
for k=1:2*N
    poles_lp(k)=Oc*exp(1i*(pi*(2*k+1)/(2*N)+pi/2));
end
poles_lp=sort(poles_lp,'ComparisonMethod','real');

%%  %analog bandpass

poles_bp=zeros(1,4*N);                             
i=1;
for k=1:2*N
    poles_bp(i)=(poles_lp(k)*B+sqrt((B*poles_lp(k))^2-4*O0^2))/2;
    poles_bp(i+1)=(poles_lp(k)*B-sqrt((B*poles_lp(k))^2-4*O0^2))/2;
    i=i+2;
end
zeros_bp=0; %order N

%%  %discrete bandpass

poles=zeros(1,2*N);                                 
for k=1:2*N
    poles(k)=(1+poles_bp(k))/(1-poles_bp(k));
end
zeros_disc=[1,-1];  %order N each







%% Paramters for BandStop to LPF Transformation
Omega_o = O0;
% B = (omega_p2-omega_p1);

%% Creating the Transfer function for the Analog Low Pass Filter
[numerator, denominator] = zp2tf([], poles_lp(1,1:N), Oc^N); % Transfer Function is multiplied with a constant to make DC gain=1 at s=0

%% Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(numerator,s)/poly2sym(denominator,s);
analog_bpf(s) = analog_lpf(((s*s)+(Omega_o*Omega_o))/(B*s));
discrete_bpf(z) = analog_bpf((z-1)/(z+1));

%% Coefficients of Analog LPF
[nums, dens] = numden(analog_lpf(s));
num_lpf = sym2poly(expand(nums));
den_lpf = sym2poly(expand(dens));
num_lpf = num_lpf./den_lpf(1);
den_lpf = den_lpf./den_lpf(1);
[H_lpf,f_lpf] = freqs(num_lpf, den_lpf, 1000);
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

%% Coefficients of Analog LPF
[nums_b, dens_b] = numden(analog_bpf(s));
num_bpf = sym2poly(nums_b);
den_bpf = sym2poly(dens_b);
num_bpf = num_bpf./den_bpf(1);
den_bpf = den_bpf./den_bpf(1);
[H_bpf,f_bpf] = freqs(num_bpf, den_bpf, 1000);
critical_points_bpf = zeros(5);
check_points_bpf = [Os(1), Op(1), O0, Op(2), Os(2)];
for i=1:5
    for l=1:length(f_bpf)
        if(f_bpf(l)>=check_points_bpf(i))
            critical_points_bpf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_bpf, abs(H_bpf)); 
title('H_{analog, BPF}(s) on normalized frequency axis'); 
xlim([0.1,3.5]); ylim([0, 1.2]); xlabel('s'); ylabel('|H_{analog,BPF}(s)|'); 
for i=1:5
    plot(f_bpf(critical_points_bpf(i)), abs(H_bpf(critical_points_bpf(i))),'r*');
end
grid on;
%% Discrete Filter Coefficients
[nums_b2, dens_b2] = numden(discrete_bpf(z));
num_bpf2 = sym2poly(nums_b2);
den_bpf2 = sym2poly(dens_b2);
num_bpf2 = num_bpf2./den_bpf2(1);
den_bpf2 = den_bpf2./den_bpf2(1);
[H_bpf2,f_bpf2] = freqz(num_bpf2, den_bpf2, 1024*1024, 330e3);
critical_points_bpf = zeros(4);
check_points_bpf = [stopband(1)*1000, passband(1)*1000, passband(2)*1000, stopband(2)*1000];
for i=1:length(critical_points_bpf)
    for l=1:length(f_bpf2)
        if(f_bpf2(l)>=check_points_bpf(i))
            critical_points_bpf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_bpf2, abs(H_bpf2)); 
title('H_{discrete, BPF}(z) on un-normalized frequency axis'); 
ylim([0, 1.2]); xlabel('Frequency in Hz'); ylabel('|H_{discrete,BPF}(z)|'); 
for i=1:length(critical_points_bpf)
    plot(f_bpf2(critical_points_bpf(i)), abs(H_bpf2(critical_points_bpf(i))),'r*');
end
grid on;
%%
fvtool(num_bpf2, den_bpf2);





%%
f_samp = 330e3;
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
delta_omega = 0.0242*pi;
Nk = ceil((A-8)/(2*2.285*delta_omega));

%% Kaiser Window formation
N_corrected = (Nk+10); %N is corrected due to poor lower bound derived previously
n = (2*(N_corrected)+1);
beta = alpha/Nk;
kaiser_window = kaiser(n,beta);

%% FIR BPF Filter formation
passband_1 = wp(1)-0.0121*pi;
passband_2 = wp(2)+0.0121*pi;
bpf = lpf_FIR(passband_2,n) - lpf_FIR(passband_1,n);
bpf_FIR = bpf.*kaiser_window';
fvtool(bpf_FIR);
%% Plotting FIR filter on unnormalized frequencies
[H_FIR, f_FIR] = freqz(bpf_FIR,1,1024,f_samp);
critical_points_bpf = zeros(4);
check_points_bpf = [stopband(1)*1000, passband(1)*1000, passband(2)*1000, stopband(2)*1000];
for i=1:length(critical_points_bpf)
    for l=1:length(f_FIR)
        if(f_FIR(l)>=check_points_bpf(i))
            critical_points_bpf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_FIR, abs(H_FIR)); 
title('H_{discrete, BPF}(z) on un-normalized frequency axis'); 
ylim([0, 1.2]); xlabel('Frequency in Hz'); ylabel('|H_{discrete,BPF}(z)|'); 
for i=1:length(critical_points_bpf)
    plot(f_FIR(critical_points_bpf(i)), abs(H_FIR(critical_points_bpf(i))),'r*');
end
grid on;
