[Y,Fs] = audioread('./MFCC/gas_station.wav');

avg = mean(Y);

Y_new = [];

for i = 1:length(Y)
  Y_tmp = [Y(i)-avg];
  Y_new = [Y_new;Y_tmp];
endfor

# prob 4
logE = 0;
for i = 7841:7841+400-1
  logE = logE + Y_new(i)^2;
endfor
logE = log(logE);
logE = max(-50,logE);

# prob 5
Y_pe = [];
for i = 7841:7841+400-1
  Y_pe_tmp = [Y_new(i) - 0.97*Y_new(i-1)];
  Y_pe = [Y_pe;Y_pe_tmp];
endfor
#plot(Y_new(7841:7841+400-1))
#plot(Y_pe)
subplot(1,2,1)
plot(abs(fft(Y_new(7841:7841+400-1))))
subplot(1,2,2)
plot(abs(fft(Y_pe)))

# prob 6
ham = hamming(512);
ham = ham.*[Y_pe(1:400);zeros(112,1)];
y_fft = fft(ham);
#plot([0:31.25:512*31.25-1],abs(y_fft))
#ylabel('spectrum')
#xlabel('frequency (Hz)')

# prob 7
load('./MFCC/mel_filters.mat')
#plot((0:256)/256*8000, mel_filters)
#title('Mel filter banks')
#xlabel('frequency (Hz)')
#ylabel('weight')
energy = [];
for i = 1:23
  tmp_e = mel_filters(:,i).*abs(y_fft(1:257));
  tmp_e = sum(tmp_e);
  energy = [energy;tmp_e];
endfor
#plot(energy)

# prob 8
cs = [];
for i = 1:13
  tmp_c = 0;
  for j = 1:23
    tmp_c = tmp_c + max(-50, log(energy(j)))*cos((pi*i/23)*(j-0.5));
  endfor
  cs = [cs;tmp_c];
endfor
#cs
# prob 9
function features = compute_mfcc(file_path,window_size,shift_size)
  # input
  # window_size : number of sample in each window
  # shift_size : number of sample in each shift
  load('./MFCC/mel_filters.mat')
  features = [];
  [Y,Fs] = audioread(file_path);
  avg = mean(Y);
  Y_new = [];
  for i = 1:length(Y)
    Y_tmp = [Y(i)-avg];
    Y_new = [Y_new;Y_tmp];
  endfor
  
  frame_number = ceil((length(Y) - window_size)/shift_size + 1)
  
  # pad zeros
  Y_new = [Y_new;zeros(shift_size*(frame_number-1) + window_size - length(Y),1)];
  
  for frame_num = 1:frame_number
    feature = [];
    
    # find logE
    logE = 0;
    for i = shift_size*(frame_num-1) + 1 : shift_size*(frame_num-1) + window_size
      logE = logE + Y_new(i)^2;
    endfor
    logE = max(-50,log(logE));
    feature = [feature;logE];
    
    # Ype
    Y_pe = [];
    for i = shift_size*(frame_num-1) + 1 : shift_size*(frame_num-1) + window_size
      Y_pe_tmp = [];
      if (i == 1)
        Y_pe_tmp = [Y_new(i)];
      else
        Y_pe_tmp = [Y_new(i) - 0.97*Y_new(i-1)];
      endif
      Y_pe = [Y_pe;Y_pe_tmp];
    endfor
    
    # fft
    ham = hamming(400);
    ham = ham.*Y_pe(1:window_size);
    ham = [ham;zeros(112,1)];
    y_fft = fft(ham);
    
    # find energy
    energy = [];
    for i = 1:23
      tmp_e = mel_filters(:,i).*abs(y_fft(1:257));
      tmp_e = sum(tmp_e);
      energy = [energy;tmp_e];
    endfor
    
    # find cs
    cs = [];
    for i = 1:13
      tmp_c = 0;
      for j = 1:23
        tmp_c = tmp_c + max(-50, log(energy(j)))*cos((pi*(i-1)/23)*(j-0.5));
      endfor
      cs = [cs;tmp_c];
    endfor

    feature = [feature;cs];
    
    # append feature to features
    features = [features feature];
  endfor
endfunction

# section 2
#all_feature = compute_mfcc('./MFCC/gas_station.wav',400,160);
#all_mean = (mean(all_feature.')).';
#all_feature = all_feature.-all_mean;
#load('./MFCC/gas_station_mfcc.mat');
#subplot(1,2,1);
#plot(gas_station_mfcc(:,54));
#subplot(1,2,2);
#plot(all_feature(:,54));
#echo = compute_mfcc('./MFCC/gas_station_echo.wav',400,160);
#echo_mean = (mean(echo.')).';
#echo = echo.- echo_mean;
#distance = sum(sum((all_feature.-echo).^2))^0.5



  
