function[trans_img] = ProcessFreqDomain(wHc, tau_step, tau_len)

N = 1000/tau_step;
sz = size(wHc);
trans_img = zeros(sz(1),sz(2),tau_len);

for m = 1:sz(1);
    for n = 1:sz(2);
    pix = wHc(m,n,:);
    pix = pix(:);
    A = [abs(pix(1));pix];
    a = ifft(A,N);
    a = 2*abs(real(a(1:tau_len)));
    trans_img(m,n,:) = permute(a,[3 2 1]);
    end
end

trans_img = trans_img/tau_step;

end