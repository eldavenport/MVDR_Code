%mvdr across time % MVDR Implementation

% R - spatial covariance matrix
% k - expected number of sources
% d0 - inter element spacing used to obtain R

% R should obviously be square

% ------------------------- OUR DATA ------------------------------------
N = 64; % num elements
fs = 1500; %hz
wavelength = 1500/250; 
spacing = 118/63;
d = spacing/wavelength;
 
data = load('vlaAcoustic64.mat');
samples = data.vlaAcoustic64.samples;

%% 

window_length = 3000;
nfft = 4096;
desired_frequency = 338; % hz
bin_number = ceil(desired_frequency / (fs/nfft)); % desF / (hz/bin)
start_time = 1;

% filter for a specific frequency, then use that data
j = 1;
for time_index = start_time:window_length:length(samples)-window_length
    data_window = samples(time_index:time_index+window_length-1, :)';

    % could add a kaiser window before doing this... 
    for i = 1:height(data_window)
        %*kaiser(window_length, 7.85)
        %data_window(i,:) = data_window(i,:);
        data_fft(i,:) = fft(data_window(i,:),nfft,2);
    end
%%
    % check this with a stem plot
    % add doppler shift compensation
    data_at_desired_bin = data_fft(:, bin_number); % 64x1
    
%    R = data_at_desired_bin*data_at_desired_bin'/N; % spatial covariance matrix 64x64 
    R = toeplitz(autocorr(data_at_desired_bin', N-1));
        
    % directions to look, if we know aperature is 120, can we do 60 to 60
    angles=(-90:.1:90);
    % steering vector to look
    a1=exp(-1i*2*pi*d*(0:N-1)'*(angles(:)'*pi/180));
    
    % inv(A)*b = A\b
    IR = inv(R);

    for k = 1:length(angles)
        mvdr(k,j) = 1/(a1(:,k)'*IR*a1(:,k));
    end

    j = j + 1;
end 

%% 
[row,col] = size(mvdr);

time_vector = 1:1:col;

for i = 1:col
    mvdr(:,i) = abs(mvdr(:,i)/max(mvdr(:,i)));
end

figure(1)
imagesc(time_vector, angles, 20*log10(mvdr))
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time'); ylabel('Angle');
colorbar;
set(gcf,'color','w')
title('MVDR')
ylim([-40 40])

