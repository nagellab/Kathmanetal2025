function A = autocorrMatrix(S,nfft,w);
% A = autocorrMatrix(S,nfft);
% calculate autocorrelation matrix at frequency omega (w)
% use fft and resampling
% A is calculated for stimulus segments of nfft, every nfft/2
% stimulus frequnecy bands should be in rows

ind = [1:nfft/2:length(S)];
n = length(ind)-1;
 
%len = 500*32;
%n = fix(size(S,1)/len); 

% % normalize S
% meanS = mean(S);
% for i=1:size(S,2),
%     S(:,i) = S(:,i)/meanS(i); 
% end

A = zeros(size(S,2));
h = repmat(hanning(nfft),1,size(S,2));
for k=1:n-1
    s = S(ind(k):ind(k)+nfft-1,:); 
    Fs = fft(s.*h,nfft);
    A = A + Fs(w,:)'*(Fs(w,:));
end
s = S(ind(n):end,:);
Fs = fft(s,nfft);
A = A + Fs(w,:)'*(Fs(w,:)); 
    
    
A = A/n;

% imagesc(abs(A)); axis square;

