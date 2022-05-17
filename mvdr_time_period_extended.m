%mvdr across time % MVDR Implementation

% R - spatial covariance matrix
% d0 - inter element spacing used to obtain R

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
        data_window(i,:) = data_window(i,:).*kaiser(window_length, 7.85)';
        data_fft(i,:) = fft(data_window(i,:),nfft,2);
    end
%%
    % check this with a stem plot
    % add doppler shift compensation
    data_at_desired_bin = data_fft(:, bin_number); % 64x1
    
    R = toeplitz(autocorr(data_at_desired_bin', N-1));
        
    % directions to look, if we know aperature is 120, can we do 60 to 60
    angles=(-90:.1:90);
    % steering vector to look
    a1=exp(-1i*2*pi*d*(0:N-1)'*(angles(:)'*pi/180));
    
    % inv(A)*b = A\b

    for k = 1:length(angles)
        mvdr(k,j) = 1/(a1(:,k)'*(R\a1(:,k)));
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
imagesc(time_vector, angles, (mvdr))
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time'); ylabel('Angle');
colorbar;
set(gcf,'color','w')
title('MVDR')
ylim([-40 40])

%%

figure(2)
subplot(1,3,1)
imagesc(time_vector*2, angles, 20*log10(mvdr_338))
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time (s)'); ylabel('Angle (deg)');
a = colorbar;
set(gcf,'color','w')
ylabel(a,'MVDR Power (dB)','FontSize',10,'Rotation',270);
a.Label.Position(1) = 3;
title('MVDR, 338 Hz')
ylim([-40 40])
subplot(1,3,2)
imagesc(time_vector*2, angles, 20*log10(mvdr_235))
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time (s)'); ylabel('Angle (deg)');
a = colorbar;
set(gcf,'color','w')
ylabel(a,'MVDR Power (dB)','FontSize',10,'Rotation',270);
a.Label.Position(1) = 3;
title('MVDR, 235 Hz')
ylim([-40 40])
subplot(1,3,3)
imagesc(time_vector*2, angles, 20*log10(mvdr_112))
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time (s)'); ylabel('Angle (deg)');
a = colorbar;
set(gcf,'color','w')
ylabel(a,'MVDR Power (dB)','FontSize',10,'Rotation',270);
a.Label.Position(1) = 3;
title('MVDR, 112 Hz')
ylim([-40 40])

figure(3)
subplot(1,2,1)
imagesc(time_vector*2, angles, 20*log10(mvdr_112))
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time (s)'); ylabel('Angle (deg)');
a = colorbar;
set(gcf,'color','w')
ylabel(a,'MVDR Power (dB)','FontSize',10,'Rotation',270);
a.Label.Position(1) = 3;
title('MVDR, 112 Hz')
ylim([-40 40])
subplot(1,2,2)
imagesc(time_vector*2, angles, 20*log10(mvdr_64))
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time (s)'); ylabel('Angle (deg)');
a = colorbar;
set(gcf,'color','w')
ylabel(a,'MVDR Power (dB)','FontSize',10,'Rotation',270);
a.Label.Position(1) = 3;
title('MVDR, 64 Hz')
ylim([-40 40])