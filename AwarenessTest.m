run_connection

FC3_data = data(:,4);
FC3_data_OPEN = FC3_data(1:1280);                       % First 10 seconds
FC3_data_CLOSED = FC3_data(1281:2560);                  % Second 10 seconds
FC3_data_TEST = FC3_data(2561:end);                     % Last 10 seconds

Ts = 1/128;                                             % Sampling Interval (s)
Fs = 1/Ts;                                              % Sampling Frequency
Fn = Fs/2;                                              % Nyquist Frequency

%%%% ALL %%%%

FC3_data(isnan(FC3_data))=[];                           % Eliminate ‘NaN’ Values First
LF = size(FC3_data,1);           
T = (0:LF-1)*Ts;                                        % Time
 
FC3ref_voltage = FC3_data(2);
FC3_plot = (FC3_data - FC3ref_voltage)*0.51;            % Remove DC offset

fcutlow = 1;                                            % low cut frequency in Hz
fcuthigh = 15;                                          % high cut frequency in Hz
[b,a] = butter(2,[fcutlow,fcuthigh]/(Fs/2));            % 2nd Order Butterworth
Filt_filt_signal = filtfilt(b,a, FC3_data);             % Filtered Signal

FF = fft(Filt_filt_signal)/LF;                          % Fourier Series of Data, Freq Vector
FFNew = FF(2:end);
Fv = linspace(0,1,fix(LF/2)+1)*Fn;
FvNew = Fv(2:end);
Iv = 1:length(FvNew);                                   % Index Vector

%%%% OPEN %%%%

FC3_data_OPEN(isnan(FC3_data_OPEN))=[];                 
LF_OPEN = size(FC3_data_OPEN,1);           
T_OPEN = (0:LF_OPEN-1)*Ts;
 
FC3ref_voltage_OPEN = FC3_data_OPEN(2);
FC3_plot_OPEN = (FC3_data_OPEN - FC3ref_voltage_OPEN)*0.51;

fcutlow = 1;     
fcuthigh = 15;   
[b,a] = butter(2,[fcutlow,fcuthigh]/(Fs/2));
Filt_filt_signal_OPEN = filtfilt(b,a, FC3_data_OPEN);

FF_OPEN = fft(Filt_filt_signal_OPEN)/LF_OPEN;            
FFNew_OPEN = FF_OPEN(2:end);
Fv_OPEN = linspace(0,1,fix(LF_OPEN/2)+1)*Fn;
FvNew_OPEN = Fv_OPEN(2:end);
Iv_OPEN = 1:length(FvNew_OPEN);         

%%%% CLOSED %%%%

FC3_data_CLOSED(isnan(FC3_data_CLOSED))=[];     
LF_CLOSED = size(FC3_data_CLOSED,1);           
T_CLOSED = (0:LF_CLOSED-1)*Ts;
 
FC3ref_voltage_CLOSED = FC3_data_CLOSED(2);
FC3_plot_CLOSED = (FC3_data_CLOSED - FC3ref_voltage_CLOSED)*0.51;

fcutlow = 1;     
fcuthigh = 15;   
[b,a] = butter(2,[fcutlow,fcuthigh]/(Fs/2));
Filt_filt_signal_CLOSED = filtfilt(b,a, FC3_data_CLOSED);

FF_CLOSED = fft(Filt_filt_signal_CLOSED)/LF_CLOSED;             
FFNew_CLOSED = FF_CLOSED(2:end);
Fv_CLOSED = linspace(0,1,fix(LF_CLOSED/2)+1)*Fn;
FvNew_CLOSED = Fv_CLOSED(2:end);
Iv_CLOSED = 1:length(FvNew_CLOSED);          

%%%% TEST %%%%

FC3_data_TEST(isnan(FC3_data_TEST))=[];     
LF_TEST = size(FC3_data_TEST,1);           
T_TEST = (0:LF_TEST-1)*Ts;
 
FC3ref_voltage_TEST = FC3_data_TEST(2);
FC3_plot_TEST = (FC3_data_TEST - FC3ref_voltage_TEST)*0.51;

fcutlow = 1;     
fcuthigh = 15;   
[b,a] = butter(2,[fcutlow,fcuthigh]/(Fs/2));
Filt_filt_signal_TEST = filtfilt(b,a, FC3_data_TEST);

FF_TEST = fft(Filt_filt_signal_TEST)/LF_TEST;           
FFNew_TEST = FF_TEST(2:end);
Fv_TEST = linspace(0,1,fix(LF_TEST/2)+1)*Fn;
FvNew_TEST = Fv_TEST(2:end);
Iv_TEST = 1:length(FvNew_TEST);          

%%%% PLOT %%%%

% Plot All Data in time domain
figure(1)                        
plot(T, Filt_filt_signal)
xlabel('Time (s)')
ylabel('F3 (µV)')
axis([0  30    ylim])
grid

% Plot 3 frequency responses ( Open, Closed, Test)
figure(2)                   
plot(FvNew_OPEN, abs(FFNew_OPEN(Iv_OPEN)))
grid
xlabel('Frequency (Hz)')
ylabel('F3 (µV)')
axis([0  65    ylim])

figure(3)                  
plot(FvNew_CLOSED, abs(FFNew_CLOSED(Iv_CLOSED)))
grid
xlabel('Frequency (Hz)')
ylabel('F3 (µV)')
axis([0  65    ylim])

figure(4)                 
plot(FvNew_TEST, abs(FFNew_TEST(Iv_TEST)))
grid
xlabel('Frequency (Hz)')
ylabel('F3 (µV)')
axis([0  65    ylim])

%%%% Welch Periodogram %%%%

win = 4 * Fs;                                       % Hamming Window size

% OPEN
p_OPEN = pwelch(Filt_filt_signal_OPEN, win, Fs);    % Plot the power spectrum
figure(5)
plot(p_OPEN)
xlabel('Frequency (Hz)')
ylabel('Power spectral density (µV^2 / Hz)')
title("Welch's periodogram")
xlim([8 13])                                        % Alpha wave band

% CLOSED
p_CLOSED = pwelch(Filt_filt_signal_CLOSED, win, Fs);
figure(6)
plot(p_CLOSED)
xlabel('Frequency (Hz)')
ylabel('Power spectral density (µV^2 / Hz)')
title("Welch's periodogram")
xlim([8 13])

% TEST
p_TEST = pwelch(Filt_filt_signal_TEST, win, Fs);
figure(7)
plot(p_TEST)
xlabel('Frequency (Hz)')
ylabel('Power spectral density (µV^2 / Hz)')
title("Welch's periodogram")
xlim([8 13])

%%%% BAND POWER %%%%

bp_OPEN = bandpower(Filt_filt_signal_OPEN,Fs,[8 13]);       % Band power Open
bp_CLOSED = bandpower(Filt_filt_signal_CLOSED,Fs,[8 13]);   % Band Power Closed
bp_TEST = bandpower(Filt_filt_signal_TEST,Fs,[8 13]);       % Band power Test

bp_ref = sqrt(bp_OPEN*bp_CLOSED);                           % Reference band power

% If statement to determine whether eyes Open or Closed
if bp_TEST >= bp_ref 
    disp("Eyes closed")
else 
    disp("Eyes Open")
end

%%%% SAVE DATA %%%%
%{
save TESTDataAll data

save TESTF3DataT_ALL T
save TESTF3DataFV_ALL FvNew
save TESTF3Data_ALL Filt_filt_signal
save TESTF3DataFF_ALL FFNew

save TESTF3DataT_OPEN T_OPEN
save TESTF3DataFV_OPEN FvNew_OPEN
save TESTF3Data_OPEN Filt_filt_signal_OPEN
save TESTF3DataFF_OPEN FFNew_OPEN

save TESTF3DataT_CLOSED T_CLOSED
save TESTF3DataFV_CLOSED FvNew_CLOSED
save TESTF3Data_CLOSED Filt_filt_signal_CLOSED
save TESTF3DataFF_CLOSED FFNew_CLOSED

save TESTF3DataT_TEST T_TEST 
save TESTF3DataFV_TEST  FvNew_TEST 
save TESTF3Data_TEST  Filt_filt_signal_TEST 
save TESTF3DataFF_TEST  FFNew_TEST 
    %}