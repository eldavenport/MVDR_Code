% MVDR Implementation

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
P = [1 1];
doas = [-30 22]*pi/180;

K = 3000;
r = 2; % total number of signals

A=exp(-1i*2*pi*d*(0:N-1)'*sin(doas(:)'));
sig=round(rand(r,K))*2-1;
noise=sqrt(1/2)*(randn(N,K)+1i*randn(N,K));  %Uncorrelated noise

% ------------------------- OUR DATA ------------------------------------
X=A*diag(sqrt(P))*sig+noise;   %Generate data matrix

R = X*X'/K; % spatial covariance matrix 10x10
%R = toeplitz(autocorr(X, N-1));

IR = inv(R);

% directions to look
angles=(-90:0.01:90);
% steering vector to look
a1=exp(-1i*2*pi*d*(0:N-1)'*sin(angles(:)'*pi/180));

for k = 1:length(angles)
    mvdr(k) = 1/(a1(:,k)'*IR*a1(:,k));
end

figure(1)
plot(angles,abs(mvdr))
