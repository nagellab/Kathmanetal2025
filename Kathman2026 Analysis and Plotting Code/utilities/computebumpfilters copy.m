function res = computebumpfilters(data,pre,post,lambda)

% close all

res.amp = [];
res.col = [];
res.odor = [];
res.fvel = [];
res.avel = [];
res.heading = [];

% trim bad data from ends and concatenate trials to form single vector for each parameter
j = 0;
for i=1:length(data) %add trials together
%     bump_col = data(i).bump_col(:,1); %bump position
    bump_amp = data(i).bump_amp(:,1); 
    if data(i).closed_loop && data(i).odor_on && data(i).wind_on && any(~isnan(data(i).plume_t)) && ~isempty(data(i).bump_amp) %only plume trials
            amp = bump_amp(401:end);
%             col = bump_col(401:end);
            odor = 10*data(i).odor(401:end);
            fvel = data(i).calc_deltapitch(401:end);%.*9.52/2;
            avel = data(i).calc_deltaz(401:end);
            heading = data(i).calc_windpos(401:end);

            ind = find(isnan(amp));
            amp(ind) = [];
%             col(ind) = [];
            odor(ind) = [];
            fvel(ind) =[];
            avel(ind) = [];
            heading(ind) = [];

            res.amp = [res.amp;amp];
%             res.col = [res.col;col];
            res.odor = [res.odor;odor];
            res.fvel = [res.fvel;fvel];
            res.avel = [res.avel;avel];
            res.heading = [res.heading;heading];
            j = j + 1;
    end
end

%% compute STAs
[res.co t] = quickfftxcorr(res.amp-mean(res.amp),res.odor-mean(res.odor),100,-pre,post);
[res.cf t] = quickfftxcorr(res.amp-mean(res.amp),(res.fvel-mean(res.fvel))/std(res.fvel),100,-pre,post);
[res.ca t] = quickfftxcorr(res.amp-mean(res.amp),(abs(res.avel)-mean(abs(res.avel)))/std(abs(res.avel)),100,-pre,post);

res.t = t;

%% decorrelate

% downsample raw STAs and inputs; form matrices
c = downsample([res.co res.cf res.ca],10); %filters
S = downsample([res.odor res.fvel abs(res.avel)],10);

% take filters to the Fourier domain
nfft = 512;
Fc = fft(c,nfft);
Fd = zeros(size(Fc'));

for i=1:257
    A = autocorrMatrix(S,nfft,i);
    [u,s,v] = svd(A);
    is = 1./diag(s) + lambda;
    Fd(:,i) = v*diag(is)*u'*Fc(i,:)';
    if i>=2 & i < nfft/2
        Fd(:,nfft+2-i) = conj(Fd(:,i));
    end
end

D = real(ifft(Fd'));
D = D(1:length(c),1:3);

res.c = c; %downsampled filters
res.S = S; %downsampled inputs
res.D = D; %decorrelated 