% Defining Params
Distance = 4e-2;        %In Meters
SoundSpeed = 340;       %In Meters per second
%AoA = 40;                % Angle of attack in degrees
maxError = 0.1;
Estimated_angle = zeros(1, 180);
AoA = (0:179);
Fsa = [100, 200, 500, 1000]*1e3;         % Sampling frequency% Sampling frequency
for j = 1:4
Fs = Fsa(j);
tend = 0.1;
t = 0:1/Fs:tend*(1-1/Fs);
for k = 1:180
% Generating the two Microphone arrays using angle data.
Distance_Difference =  Distance*cos(AoA(k)*pi/180);
Actual_Time_delay = Distance_Difference/SoundSpeed;
MicA = sin(2*pi*10*(t));
MicB = sin(2*pi*10*(t+Actual_Time_delay));
 %Plotting
 %subplot(2,1,1);
 %plot(MicA);
 %subplot(2,1,2);
 %plot(MicB);
% Adding noise
maxError = 0.01;
Error_signal_MicA = 2*(randn(size(MicA))-0.5)*maxError;
Error_signal_MicB = 2*(rand(size(MicA))-0.5)*maxError;
%plot(Error_signal_MicA);
MicA = MicA+Error_signal_MicA;
MicB = MicB+Error_signal_MicB;

% Decoding Time Delay
c = xcorr(MicA, MicB);
L = length(t);
[maxnum, id] = max(c);
diff = (id - L)/L * tend;
error = (Actual_Time_delay - diff)/Actual_Time_delay;

% Getting angle Data
Estimated_angle(k) = 180*abs(acos(SoundSpeed*diff/Distance))/pi;


end
subplot(2,2,j);
plot(AoA);
hold on
plot(Estimated_angle);
ylabel('Angle');
title(strcat("sampling frequncy: ", mat2str(Fs), " Hertz"));
legend('AoA', 'Est Angle');
hold off
end
display(strcat("With maximum error of ",  mat2str(maxError*100) , "%"));