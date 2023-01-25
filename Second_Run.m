clearvars;
syms H z t k x;



%First Question
t2 = -2*pi:0.01:2*pi;

Hz(z)=(z-exp(i*(-15*pi/40)))*(z-exp(i*(15*pi/40)))*(z-exp(i*(-2*pi/4)))*(z-exp(i*(2*pi/4)))*(z-1)*(z+1)*(z-exp(i*(30*pi/40)))*(z-exp(i*(-30*pi/40)));
H(z) = Hz*(1/100)/((z-(0.95)*exp(i*(-(0.75))))*(z-(0.95)*exp(i*((0.75))))*(z-(0.95)*exp(i*(-(0.692))))*(z-(0.95)*exp(i*((0.692))))*(z-(0.9)*exp(i*(-(0.592))))*(z-(0.9)*exp(i*((0.592))))*(z-(0.9)*exp(i*(-(0.492))))*(z-(0.9)*exp(i*((0.492))))*(z-(0.95)*exp(i*(-(0.4))))*(z-(0.95)*exp(i*((0.4)))));
Hp(z) = 100*((z-(0.95)*exp(i*(-(0.75))))*(z-(0.95)*exp(i*((0.75))))*(z-(0.95)*exp(i*(-(0.692))))*(z-(0.95)*exp(i*((0.692))))*(z-(0.9)*exp(i*(-(0.592))))*(z-(0.9)*exp(i*((0.592))))*(z-(0.9)*exp(i*(-(0.492))))*(z-(0.9)*exp(i*((0.492))))*(z-(0.95)*exp(i*(-(0.4))))*(z-(0.95)*exp(i*((0.4)))));


k= 1/100;
zeroz = [exp(i*(-15*pi/40)); exp(i*(15*pi/40)); exp(i*(-2*pi/4)); exp(i*(2*pi/4)); 1; -1; exp(i*(30*pi/40)); exp(i*(-30*pi/40));0;0];

poles = [(0.95)*exp(i*(-(0.75))); (0.95)*exp(i*((0.75))); (0.95)*exp(i*(-(0.692))); (0.95)*exp(i*((0.692))); (0.9)*exp(i*(-(0.592))); (0.9)*exp(i*((0.592))); (0.9)*exp(i*(-(0.492))); (0.9)*exp(i*((0.492))); (0.95)*exp(i*(-(0.4))); (0.95)*exp(i*((0.4)))];

figure(1)
plot(t2,abs(H(exp(j*t2))));
title("System function of the given Band-Pass filter")
ylabel("|H(e^j^w)|")
xlabel("w")
clear n;

figure(2)
plot(t2,angle(H(exp(j*t2))));
title("Phase of H(e^j^w)")
xlabel("w")
ylabel("<H(e^j^w)")


figure(3)
[B,A] = zp2tf(zeroz,poles,1.05);
zplane(zeroz,poles)



% Hzero = double(coeffs(Hz,z,"All"));
% Hpole = double(coeffs(Hp,z,'All'));
figure(4);
[hn,n] = impz(B,A,-40:130);
stem(n,hn)
title("h[n] for designated X(jw)")
xlabel("n (samples)")
ylabel("h[n]")
clear n;
syms n;
int3rval1 = 1:1:1034;
int3rval2 = 1:1:8204;
chirp1 = cos((int3rval1.^2)*pi/512);
chirp2 = cos((int3rval2.^2)*pi/8191);

%difference eq'ns
clear x n y;
syms x n;

A = A./2000;
y = zeros(1,1034);
chirp1 = [0 0 0 0 0 0 0 0 0 0 chirp1];
for i = -10:1023
    if i>1
    y(i+10) = ((-y(i+9))*A(1)-(y(i+8))*A(2)-(y(i+7))*A(3)-(y(i+6))*A(4)-(y(i+5))*A(5)-(y(i+4))*A(6)-(y(i+3))*A(7)-(y(i+2))*A(8)-(y(i+1))*A(9)-(y(i))*A(10)-(y(i-1))*A(11)) + (chirp1(i+11)*B(1)+chirp1(i+10)*B(2)+chirp1(i+9)*B(3)+chirp1(i+8)*B(4)+chirp1(i+7)*B(5)+chirp1(i+6)*B(6)+chirp1(i+5)*B(7)+chirp1(i+4)*B(8)+chirp1(i+3)*B(9)+chirp1(i+2)*B(10)+chirp1(i+1)*B(11));
    else
        y(i+11) = 0;
    end

end
output1 = 20*y;

y = zeros(1,8204);
chirp2 = [0 0 0 0 0 0 0 0 0 0 chirp2];
for i = -10:8191
    if i>1
    y(i+10) = ((-y(i+9))*A(1)-(y(i+8))*A(2)-(y(i+7))*A(3)-(y(i+6))*A(4)-(y(i+5))*A(5)-(y(i+4))*A(6)-(y(i+3))*A(7)-(y(i+2))*A(8)-(y(i+1))*A(9)-(y(i))*A(10)-(y(i-1))*A(11)) + (chirp2(i+11)*B(1)+chirp2(i+10)*B(2)+chirp2(i+9)*B(3)+chirp2(i+8)*B(4)+chirp2(i+7)*B(5)+chirp2(i+6)*B(6)+chirp2(i+5)*B(7)+chirp2(i+4)*B(8)+chirp2(i+3)*B(9)+chirp2(i+2)*B(10)+chirp2(i+1)*B(11));
    else
        y(i+11) = 0;
    end
end
output2 = 20*y;
%store .mat

figure(5)
stem(int3rval1,output1,'filled','.k')
axis tight
title("Plot of X_1[n] created by iteration algorythm")
xlabel("n (samples)")
ylabel("X_1[n]")
save("X1[n].mat","output1")

 figure(6)
 stem(int3rval2,output2,'filled')
 axis tight
 title("Plot of X_2[n] created by iteration algorythm")
xlabel("n (samples)")
ylabel("X_2[n]")
save("X2[n].mat","output2")

%Fourth Q

alpha = 1000;
len1 = 8192;

Tsa = sqrt(pi/(len1.*alpha));
SampRate = 1./Tsa;
audioOut1 = audioplayer(output2,SampRate); SampPeriod = Tsa.*length(output2);
while(1)
  play(audioOut1);
  pause(SampPeriod-0.1);

  stop(audioOut1);
end

yrt = pchip(int3rval2,output2,int3rval2);
figure(7)
plot(int3rval2,yrt)
title("Interpolated continuous version of Y_2[n]")
xlabel("t")
ylabel("Y_R(t)")
save("continuousSound.mat","yrt")
save("yrt.mat","yrt")

audioOut2 = audioplayer(yrt,SampRate);
SampPeriod2 = Tsa.*length(yrt);

intervalLast = -2*pi:Tsa:2*pi;
figure(8)
axis tight; plot(intervalLast,abs(H(exp(j*intervalLast))))
title("Magnitude of |H_e_q(jw)|")
xlabel("w")
ylabel("|H_e_q(jw)|")


% while(1)
%   play(audioOut2);
%   pause(SampPeriod2-0.1);
% 
%   stop(audioOut2);
% end

hn = hn./1000;

%Eighth Question
[y,Fs] = audioread("Rammstein - Deutschland.m4a");
cropMusic = (0.05)*y(1:2000000);
OutputMusic3 = conv(hn,cropMusic);
audiowrite("OutputMusic3.mp4",OutputMusic3,Fs)

%Nineth Question
[y9,Fs] = audioread("WhatsApp Audio 2022-12-09 at 11.19.35.mp4");
OutputMusic4 = conv(hn,y9(1:158340));
audiowrite("OutputMusic4.mp4",OutputMusic4,Fs)