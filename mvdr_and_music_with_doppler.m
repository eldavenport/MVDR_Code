%mvdr across time % MVDR Implementation

% R - spatial covariance matrix
% d0 - inter element spacing used to obtain R


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
% desired_frequency = [79, 130, 235, 338] ; % hz
% for i = 1:4 
% y1 = mvdr(desired_frequency(i), samples)  ;
% y2 = mvdr_withds(desired_frequency(i), samples)  ;
% ym1 = music(desired_frequency(i), samples) ; 
% ym2 = music_withds(desired_frequency(i), samples) ; 
% end 

ym1 = mvdr(338, samples) ; 
ym2 = mvdr_withds(338, samples) ; 
y1 = music(338, samples)  ;
y2 = music_withds(338, samples)  ;


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

    % could add a kaiser window before doing this... 
    for i = 1:size(data_window,1)
        data_window(i,:) = data_window(i,:).*kaiser(window_length, 7.85)';
        data_fft(i,:) = fft(data_window(i,:),nfft,2);
    end
%
    % check this with a stem plot
    % add doppler shift compensation
    % dont implement doppler shift, find the doppler shift 
    data_dp = data_fft(:, bin_number-5:bin_number+5) ; % approx 5 degrees in either direction
    [~, ind] = max(sum(abs(data_dp ).^2)); % calculate energy and find largest bin 
    
    data_at_desired_bin = data_dp(:, ind); % 64x1
    % select closest bin with largest magnitude
    
    R = toeplitz(autocorr(data_at_desired_bin', N-1));
        
    % directions to look, if we know aperature is 120, can we do 60 to 60
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
imagesc(time_vector*2, angles, 20*log10(y))
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time (s)'); ylabel('Angle (deg)');
colorbar;
set(gcf,'color','w')
title('MVDR, Doppler shift '+ string(desired_frequency) + 'Hz')
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

    % could add a kaiser window before doing this... 
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
        
    % directions to look, if we know aperature is 120, can we do 60 to 60
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
A = imagesc(time_vector*2, angles, 20*log10(y)) ; 
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time (s)'); ylabel('Angle (deg)');
colorbar;
set(gcf,'color','w')
title('MVDR ' + string(desired_frequency) + 'Hz')
imwrite(255*y,jet,'MVDR' + string(desired_frequency)+ '.jpg')
ylim([-40 40])
end 


function y = music_withds(desired_frequency, samples) 
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

    % could add a kaiser window before doing this... 
    for i = 1:size(data_window,1)
        data_window(i,:) = data_window(i,:).*kaiser(window_length, 7.85)';
        data_fft(i,:) = fft(data_window(i,:),nfft,2);
    end

    % check this with a stem plot
    % add doppler shift compensation
    % dont implement doppler shift, find the doppler shift 
    data_dp = data_fft(:, bin_number-5:bin_number+5) ; % approx 5 degrees in either direction
    [~, ind] = max(sum(abs(data_dp ).^2)); % calculate energy and find largest bin 
    
    data_at_desired_bin = data_dp(:, ind); % 64x1
    % select closest bin with largest magnitude
    
    
    R = toeplitz(autocorr(data_at_desired_bin', N-1));
        
    % directions to look, if we know aperature is 120, can we do 60 to 60
    angles=(-90:.1:90);
    % steering vector to look
    a1=exp(-1i*2*pi*d*(0:N-1)'*(angles(:)'*pi/180));

    [Q, D] = eig(R); %eigenvalues and vectors of cov matrix
    [D, I] = sort(diag(D),1,'descend');
    
    r = 6; % total number of signals
    Q = Q(:,I); % sorts the eigenvectors to get signal first
    Qs = Q(:, 1:r); % signal eigenvectors
    Qn = Q(:,r+1:N); % noise eigenvectors

    for k=1:length(angles)  %Compute MUSIC 
        y(k,j)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
    end

    j = j + 1;
end 

[~,col] = size(y);

time_vector = 1:1:col;

for i = 1:col
    y(:,i) = abs(y(:,i)/max(y(:,i)));
end

figure()
A = imagesc(time_vector*2, angles, 20*log10(y)) ; 
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time (s)'); ylabel('Angle (deg)');
colorbar;
set(gcf,'color','w')
title('MUSIC, Doppler Shift ' + string(desired_frequency) + 'Hz')
imwrite(255*y,jet,'MUSICwithDS' + string(desired_frequency)+ '.jpg')
ylim([-40 40])
end 

function y = music(desired_frequency, samples) 
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

    % could add a kaiser window before doing this... 
    for i = 1:size(data_window,1)
        data_window(i,:) = data_window(i,:).*kaiser(window_length, 7.85)';
        data_fft(i,:) = fft(data_window(i,:),nfft,2);
    end

    % check this with a stem plot
    % add doppler shift compensation
    data_at_desired_bin = data_fft(:, bin_number); % 64x1
    
    R = toeplitz(autocorr(data_at_desired_bin', N-1));
        
    % directions to look, if we know aperature is 120, can we do 60 to 60
    angles=(-90:.1:90);
    % steering vector to look
    a1=exp(-1i*2*pi*d*(0:N-1)'*(angles(:)'*pi/180));

    [Q, D] = eig(R); %eigenvalues and vectors of cov matrix
    [D, I] = sort(diag(D),1,'descend');
    
    r = 6; % total number of signals
    Q = Q(:,I); % sorts the eigenvectors to get signal first
    Qs = Q(:, 1:r); % signal eigenvectors
    Qn = Q(:,r+1:N); % noise eigenvectors

    for k=1:length(angles)  %Compute MUSIC 
        y(k,j)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
    end

    j = j + 1;
end 

[~,col] = size(y);

time_vector = 1:1:col;

for i = 1:col
    y(:,i) = abs(y(:,i)/max(y(:,i)));
end

figure()
A = imagesc(time_vector*2, angles, 20*log10(y)) ; 
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time (s)'); ylabel('Angle (deg)');
colorbar;
set(gcf,'color','w')
title('MUSIC ' + string(desired_frequency) + 'Hz')
imwrite(255*y,jet,'MUSIC' + string(desired_frequency)+ '.jpg')
ylim([-40 40])
end 