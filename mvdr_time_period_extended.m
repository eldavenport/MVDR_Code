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


% filter for a specific frequency, then use that data
j = 1;
for time_index = 3000*100:500:length(samples)-500
    X = samples(time_index:time_index+500, :)';
    for i = 1:height(X)
        X(i,:) = bandpass(X(i,:), [233 236], fs, 'Steepness',0.9);
    end
    % changing steepness has what affect? why? does it makes magnitude of mvdr
    % output small, also causes singular/badly conditioned R. steepness of 0.9
    % seems good
    
    
    % TODO - filter this data better (using a window fxn) to determine how many signals we are
    % getting - can we isolate from the signals at 235, 238, and 241? our
    % ability to separate these signals will tell us how many signals there are
    % in the MVDR portion - the answer is we can't really .... so should we
    % treat all 5 signals as one signal?? because they are so close in
    % frequency...
    % TODO - look at how the angle(s) of arrival change(s) over time (it does)
    % by running this for successive sections
    % TODO - how does this work if we use another (set of) frequency(ies)?
    % Note - the frequencies are in sets of 5, spaced 3 hz apart, why? and can
    % we legitimately distinguish them (that will take some fourier analysis
    % before MVDR)
    
    % windowing, how to do this in a vectorized way? 
    % for i = 1:height(X)
    %     X(i,:) = X(i,:).*kaiser(length(X),7.85)';
    % end
    
    % % FFT of filtered data
    
    % freq_x = fft(X,2048,2);
    % normalized_f = 0:1/2048:(1 - 1/2048);
    % 
    % figure(1)
    % hold on
    % % plot(normalized_f, 20*log10(abs(freq_x(1,:))))
    % % plot(normalized_f, 20*log10(abs(freq_x(15,:))))
    % plot(normalized_f, 20*log10(abs(freq_x(30,:))))
    % % plot(normalized_f, 20*log10(abs(freq_x(45,:))))
    % % plot(normalized_f, 20*log10(abs(freq_x(60,:))))
    % hold off
    % xlim([0 1])
    % xlabel('Normalized Frequency (rad/sample)')
    % ylabel('FFT Magnitude, After BPF 235 Hz')
    
    K = 500;
    r = 1; % total number of signals
    % ------------------------- OUR DATA ------------------------------------
    
    R = X*X'/K; % spatial covariance matrix 10x10
    
    [Q, D] = eig(R); %eigenvalues and vectors of cov matrix
    [D, I] = sort(diag(D),1,'descend');
    
    Q = Q(:,I); % sorts the eigenvectors to get signal first
    Qs = Q(:, 1:r); % signal eigenvectors
    Qr = Q(:,r+1:N); % noise eigenvectors
    
    IR = inv(R);
    
    % directions to look
    angles=(-90:0.1:90);
    % steering vector to look
    a1=exp(-1i*2*pi*d*(0:N-1)'*sin(angles(:)'*pi/180));
    
    for k = 1:length(angles)
        mvdr(k,j) = 1/(a1(:,k)'*IR*a1(:,k));
    end
    j = j + 1;
end 

time_vector = 1:1:21;

figure(2)
imagesc(time_vector, angles, abs(mvdr))
set(gca,'ydir','normal'); colormap(jet);
xlabel('Time'); ylabel('Angle');
colorbar;
set(gcf,'color','w')
title('MVDR')

