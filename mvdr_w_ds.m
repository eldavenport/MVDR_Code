%mvdr across time % MVDR Implementation

% R - spatial cov matrix
% d0 - element spacing used to obtain R


% time comparison mvdr vs MUSIC 
% implement accounting for doppler effect 
% number of signals 
%% ------------------------- OUR DATA ------------------------------------
N = 64; % num elements
fs = 1500; %hz
wavelength = 1500/250; 
spacing = 118/63;
d = spacing/wavelength;
 
data = load('vlaAcoustic64.mat');
samples = data.vlaAcoustic64.samples;



%% 


%% 
desired_frequency = [79, 130, 235, 338] ; % hz
for i = 1:4 
y1 = mvdr(desired_frequency(i), samples)  ;
y2 = mvdr_withds(desired_frequency(i), samples)  ;
end 


function y = mvdr_withds(desired_frequency, samples) 
nfft = 4096;
N = 64 ; 
fs = 1500; %hz
wavelength = 1500/250; 
spacing = 118/63;
d = spacing/wavelength;
bin_number = ceil(desired_frequency / (fs/nfft)); % desF / (hz/bin)
start_time = 1;

window_length = 3000;
% filter for a specific frequency, then use that data
j = 1;
for time_index = start_time:window_length:length(samples)-window_length
    
    data_window = samples(time_index:time_index+window_length-1, :)';

    for i = 1:size(data_window,1)
        data_window(i,:) = data_window(i,:).*kaiser(window_length, 7.85)';
        data_fft(i,:) = fft(data_window(i,:),nfft,2);
    end
%
    % check this with a stem plot
    % add doppler shift compensation
    % dont implement doppler shift, find the doppler shift 
    data_dp = data_fft(:, bin_number-13:bin_number+13) ; % approx 5 degrees in either direction
    [~, ind] = max(sum(abs(data_dp ).^2)); % calculate energy and find largest bin 
    
    data_at_desired_bin = data_dp(:, ind); % 64x1
    % select closest bin with largest magnitude
    
    R = toeplitz(autocorr(data_at_desired_bin', N-1));
        
    % look directions
    angles=(-90:.1:90);
    % steering vector to look
    a1=exp(-1i*2*pi*d*(0:N-1)'*(angles(:)'*pi/180));
    
    % inv(A)*b = A\b

    for k = 1:length(angles)
        y(k,j) = 1/(a1(:,k)'*(R\a1(:,k)));
    end

    j = j + 1;
end 

[~,col] = size(y);

time_vector = 1:1:col;

for i = 1:col
    y(:,i) = abs(y(:,i)/max(y(:,i)));
end

figure()
imagesc(time_vector, angles, (y))
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time'); ylabel('Angle');
colorbar;
set(gcf,'color','w')
title('MVDR with Doppler shift'+ string(desired_frequency))
imwrite(255*y,jet,'MVDRwDS' + string(desired_frequency)+ '.jpg')
ylim([-40 40])
end 

function y = mvdr(desired_frequency, samples) 
nfft = 4096;
N = 64 ; 
fs = 1500; %hz
wavelength = 1500/250; 
spacing = 118/63;
d = spacing/wavelength;
bin_number = ceil(desired_frequency / (fs/nfft)); % desF / (hz/bin)
start_time = 1;

window_length = 3000;
% filter for a specific frequency, then use that data
j = 1;
for time_index = start_time:window_length:length(samples)-window_length
    
    data_window = samples(time_index:time_index+window_length-1, :)';

    for i = 1:size(data_window,1)
        data_window(i,:) = data_window(i,:).*kaiser(window_length, 7.85)';
        data_fft(i,:) = fft(data_window(i,:),nfft,2);
    end
%
    % check this with a stem plot
    % add doppler shift compensation
    % dont implement doppler shift, find the doppler shift 
    
    data_at_desired_bin = data_fft(:, bin_number); % 64x1
    % select closest bin with largest magnitude
    
    R = toeplitz(autocorr(data_at_desired_bin', N-1));
        
    % look directions
    angles=(-90:.1:90);
    % steering vector to look
    a1=exp(-1i*2*pi*d*(0:N-1)'*(angles(:)'*pi/180));
    
    % inv(A)*b = A\b

    for k = 1:length(angles)
        y(k,j) = 1/(a1(:,k)'*(R\a1(:,k)));
    end

    j = j + 1;
end 

[~,col] = size(y);

time_vector = 1:1:col;

for i = 1:col
    y(:,i) = abs(y(:,i)/max(y(:,i)));
end

%figure()
%A = imagesc(time_vector, angles, (y)) ; 
%set(gca,'ydir','normal'); colormap(jet);
%xlabel('Time'); ylabel('Angle');
%colorbar;
%set(gcf,'color','w')
%title('MVDR' + string(desired_frequency))
imwrite(255*y,jet,'MVDR' + string(desired_frequency)+ '.jpg')
%ylim([-40 40])
end 