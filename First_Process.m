syms x n H k y z w t Xa;


interval = linspace(pi/10,5*pi,250);
% z = exp(j*w);

%First Question
t2 = -2*pi:0.01:2*pi;
n = 0:11;
H1=(1.3*z-1.3*exp(i*(-8*pi/40)))*(1.3*z-1.3*exp(i*(8*pi/40)))*(1.3*z-0.98*exp(i*(-6*pi/40)))*(1.3*z-0.98*exp(i*(6*pi/40)))*(1.3*z-0.98*exp(i*(-9*pi/40)))*(1.3*z-0.98*exp(i*(9*pi/40)))*(1.3*z-1.2*exp(i*(-10*pi/40)))*(1.3*z-1.2*exp(i*(10*pi/40)))*(1.3*z-1.2*exp(i*(-5*pi/40)))*(1.3*z-1.2*exp(i*(5*pi/40)));
H(z) =(1.3*z)*(1.3*z-0.2*exp(i*(pi/10)))*(1.3*z-0.2*exp(-i*(pi/10)))*(1.3*z-0.2*exp(i*(2*pi/4)))*(1.3*z-0.2*exp(-i*(2*pi/4)))/H1;
filterResult = double(abs(H(j*t2)));
%*(z-exp(i*(7*pi/8)))*(z-exp(i*(-7*pi/8)))*(z-exp(i*(11*pi/12)))*(z-exp(i*(-11*pi/12)))*(z-exp(i*(0)));
% hn=coeffs(H1);
% H(z) = H1/(z^(11));
% figure(2)
% stem(n,hn)
% title("Impulse response of the system function")
% xlabel("n")
% ylabel("h[n]")
% figure(2)
% plot(double(H1))
figure(1)
plot(t2,filterResult);
title("System function of the given Band-Pass filter")
ylabel("|H(e^j^w)|")
xlabel("w")
clear n;

% figure(3)
% plot(t2,angle(H(exp(j*t2))))
% title("System function of the given Band-Pass filter's phase plot")
% ylabel("Phase angle of H(e^j^w)")
% xlabel("w")


%Plotting the zeros on z-plane, Second Question
z0 = [-pi/8,1];
z1 = [pi/8,1];
z2 = [-pi*60/100,1];
z3 = [pi*60/100,1];
z4 = [-pi*3/4,1]; 
z5 = [pi*3/4,1]; 
z6 = [-pi*7/8,1];
z7 = [pi*7/8,1]; 
z8 = [-pi*11/12,1];
z9 = [pi*11/12,1];
z10 = [0,1]; 
figure(4);
polarscatter([z0(1),z1(1),z2(1),z3(1),z4(1),z5(1),z6(1),z7(1),z8(1),z9(1),z10(1)],1);
title("Zeros Plot of H(e^j^w)")
% ylabel("Im{Z}")
% xlabel("Re{Z}")

% %Second Question
% syms n;
% Xa(t) = cos(1000.*(t^2));
% 
% interval2 = linspace(1,1023,1023);
% 
% X(n) = cos((pi/512)*n.^(2));
% Xinp = double(X(interval2));
% figure(5)
% stem(interval2,X(interval2),'filled','.k')
% title("X_1[n], similar to X_a(t)")
% xlabel("n")
% ylabel("X_1[n]")
% save("X(n).mat","Xinp")
% clear n;
% 
% %Third Question
% syms n;
% interval3 = linspace(1,8191,8191);
% X2(n) = cos((pi/8191)*n.^(2));
% Xinp2 = double(X2(interval3));
% figure(6)
% stem(interval3,X2(interval3))
% title("Sampled X_a(t) with designated Ts")
% ylabel("X_g[n]")
% xlabel("n")
% save("X2(n).mat","Xinp2")
% 
% %Fourth Question
% interval4 = linspace(1,1033,1034);
% hn2 = double(real(hn));
% y1 = conv(Xinp,hn2);
% figure(7)
% stem(interval4,y1)
% title("Frequency response of the FIR Filter with input X_f(n]")
% ylabel("Y_1[n]")
% xlabel("n")
% 
% interval5 = linspace(1,8201,8202);
% % XaFFT = abs(fft(Xa(interval5)));
% hn2 = double(real(hn));
% y2 = conv(Xinp2,hn2);
% figure(8)
% stem(interval5,y2)
% title("Frequency response of the FIR Filter with input X_g[n]")
% ylabel("Y_2[n]")
% xlabel("n")
% save("y2.mat","y2")
% 
% 
% %Fifth Question
% alpha = 1000;
% len1 = 8192;
% load("X2(n).mat");
% Tsa = sqrt(pi/(len1.*alpha));
% SampRate = 1./Tsa;
% audioOut1 = audioplayer(Xinp2,SampRate);
% SampPeriod = Tsa.*length(Xinp2);
% 
% 
% % while(1)
% %     play(audioOut1);
% %     pause(SampPeriod-0.1);
% %     stop(audioOut1);
% % end
% 
% load("y2.mat");
% audioOut2 = audioplayer(y2,SampRate);
% SampPeriod2 = Tsa.*length(y2);
% 
% % 
% % while(1)
% %     play(audioOut2);
% %     pause(SampPeriod2-0.1);
% %     stop(audioOut2);
% % end
% 
% %Seventh Question
% yrt = pchip(interval5,y2,interval5);
% figure(9)
% plot(interval5,yrt)
% title("Interpolated continuous version of Y_2[n]")
% xlabel("w")
% ylabel("Y_R(t)")
% save("continuousSound.mat","yrt")
% 
% load("continuousSound.mat");
% audioOut3 = audioplayer(yrt,SampRate);
% SampPeriod3 = Tsa.*length(yrt);
% % 
% % while(1)
% %     play(audioOut3);
% %     pause(SampPeriod3-0.1);
% %     stop(audioOut3);
% % end
% % 
% % T = 1/SampRate;           % Sampling period
% % t3 = 0:T:2*pi;
% 
% %Seventh Question Cont'd
% intervalLast = -2*pi:Tsa:2*pi;
% figure(10)
% axis tight;
% plot(intervalLast,abs(H(exp(j*intervalLast))))
% title("Magnitude of |H_e_q(jw)|")
% xlabel("w")
% ylabel("|H_e_q(jw)|")
% 
% 
% 
% %Eighth Question
% [y,Fs] = audioread("Rammstein - Deutschland.m4a");
% cropMusic = (0.3)*y(1:2000000);
% OutputMusic3 = conv(hn2,cropMusic);
% audiowrite("OutputMusic3.mp4",OutputMusic3,Fs)
% 
% %Nineth Question
% [y9,Fs] = audioread("WhatsApp Audio 2022-12-09 at 11.19.35.mp4");
% OutputMusic4 = conv(hn2,y9(1:158340));
% audiowrite("OutputMusic4.mp4",OutputMusic4,Fs)









