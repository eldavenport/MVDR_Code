% MVDR Implementation

% R - spatial cov matrix
% k - expected number of sources
% d0 - element spacing

% ------------------------- OUR DATA ------------------------------------
N = 64; % num elements
fs = 1500; %hz
wavelength = 1500/250; 
spacing = 118/63;
d = spacing/wavelength;
 
data = load('vlaAcoustic64.mat');
samples = data.vlaAcoustic64.samples;


window_length = 3000;
nfft = 4096;
desired_frequency = 338; % hz
bin_number = ceil(desired_frequency / (fs/nfft)); % desF / (hz/bin)
start_time = floor(length(samples)/2);

% filter for a specific frequency, then use that data
data_window = samples(start_time:start_time+window_length-1, :)';

for i = 1:height(data_window)
    %data_window(i,:) = data_window(i,:)*kaiser(window_length, 7.85);
    data_fft(i,:) = fft(data_window(i,:),nfft,2);
end

%%

% add doppler shift compensation
data_at_desired_bin = data_fft(:, bin_number); % 64x1

r = 6; % total number of signals
% ------------------------- OUR DATA ------------------------------------

%R = data_at_desired_bin*data_at_desired_bin'/N;
R = toeplitz(autocorr(data_at_desired_bin', N-1));

[Q, D] = eig(R); %eigenvalues and vectors of cov matrix
[D, I] = sort(diag(D),1,'descend');

Q = Q(:,I); % sorts the eigenvectors to get signal first
Qs = Q(:, 1:r); % signal eigenvectors
Qn = Q(:,r+1:N); % noise eigenvectors
    
% directions to look, if we know aperature is 120, can we do 60 to 60
angles=(-90:0.1:90);
% steering vector to look
a1=exp(-1i*2*pi*d*(0:N-1)'*(angles(:)'*pi/180));

% inv(A)*b = A\b

for k = 1:length(angles)
    mvdr(k) = (a1(:,k)'*a1(:,k))/(a1(:,k)'*(R\a1(:,k)));
    %mvdr(k) = 1/(a1(:,k)'*pinv(R)*a1(:,k));
    music(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
end

figure(2)
hold on
plot(angles,20*log10(abs(music/max(music))))
plot(angles,20*log10(abs(mvdr/max(mvdr))))
hold off
xlim([-40 40])
xlabel('Arrival Angle (deg)')
ylabel('Beamformer output (dB)')
set(gcf,'color','w')
legend('MUSIC','MVDR')

