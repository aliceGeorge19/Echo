[y,Fs] = audioread('hello.mp3'); %Reading the audio file and storing it in y with sampling frequency Fs.

ylen=length(y);                 % length of the audio in samples.
T = ylen/Fs;                    % Time needed to play the sound completely. 
time_axis_input = ((1:ylen)-1)/Fs;      % Define the time axis of input.
t1 = zeros(size(y,1),1);         % Initializing an empty array to generate an impsule train for convolution.
t2 = zeros(size(y,1),1);
t3 = zeros(size(y,1),1);

Delay1 = [0.2 0.8 ; 0.3 0.7 ; 0.4 0.6];  % The delay matrix consists of time of impusle (first column)
                                        % and signal amplitude (second column). 
                                        
% To generate strong/weak signals we only need to change the value of the
% second column. However, it must lie between 1 (strongest) and 0(weakest).
Delay2 = [0.2 0.5 ; 0.3 0.3 ; 0.4 0.1]; % Medium strength signal.
Delay3 = [0.2 0.1 ; 0.3 0.15 ; 0.4 0.0]; % Weak signal.


t1(1)=1;                                 % Generating the impulse train.
t2(1)=1;
t3(1)=1;
for i = 1:size(Delay1,1)
    t1(Delay1(i,1)*Fs) = Delay1(i,2);
end

for i = 1:size(Delay2,1)
    t2(Delay2(i,1)*Fs) = Delay2(i,2);
end

for i = 1:size(Delay3,1)
    t3(Delay3(i,1)*Fs) = Delay3(i,2);
end

% Now we need to convolove the audio signal y with the time axis t to find
% the generated echo signal y_out.

ly = length(y);
lt1 = length(t1);
lt2 = length(t2);
lt3 = length(t3);

outlength1 = ly + lt1 - 1;               % This is the total length of the new echo signal y_out.
outlength2 = ly + lt2 - 1;
outlength3 = ly + lt3 - 1;

%y_out = conv(y(1),t(1)); writing this statement as Fourier transform.
%(Matlab online doesn't have the full features of the offline version)

y_out1 = ifft(fft(y, outlength1) .* fft(t1, outlength1));
y_out2 = ifft(fft(y, outlength2) .* fft(t2, outlength2));
y_out3 = ifft(fft(y, outlength3) .* fft(t3, outlength3));

% In the previous statement, Inverse Fourier Transform is applied to the
% Fourier Transform of y & t which ressembles the convolution of the two
% signals.
        
sound(y_out1,Fs)  % Playing the echo signal y_out.
audiowrite('output1.wav',y_out1,Fs); % Saving the echo signal y_out "output1" file in 'wav' format.

sound(y_out2,Fs)
audiowrite('output2.wav',y_out2,Fs);

sound(y_out3,Fs)
audiowrite('output3.wav',y_out3,Fs);

y_original = ifft(fft(y_out1, outlength1) ./ fft(t1, outlength1)); 

% Now we can easily get the original hello file
% by dividing the echoed file by the time axis signal.

sound(y_original, Fs)                % The original signal.
audiowrite('output4.wav',y_original,Fs);    % Saving the original signal "output2" file in 'wav' format.