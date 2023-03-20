T=1.;  % pulse length in seconds
fs = 3000;  % sampling rate in samples per second
A = 1;  % Pulse Amplitude
fc = 920;  % Pulse Start Frequency (Hz)
B = 150;  % Pulse Sweep Bandwidth (Hz)

RR = 100;  % Range Rate in kts (+ = opening, - = closing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NT = round(T*fs);
t = (0:(NT-1))/fs;

x = A*sin( 2*pi*( fc*t + (B/T/2)*(t.^2) ) );

x_alpha = arbitrary_alpha( x, 1 - RR/3000 ); % C is about 3000 kts

t_alpha = (0:(length(x_alpha)-1))/fs;


plot( t, x, t_alpha, x_alpha+3, '--r')

X = abs(fft(x,4096));
X = X(1:2048);
X_alpha = abs(fft(x_alpha,4096));
X_alpha = X_alpha(1:2048);
f = (0:2047)*(fs/4096);
figure
plot( f, X, f, X_alpha, '--r')