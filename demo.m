load ar_data
lambda=0.1;
W=RSR(data',lambda);
figure
wc=sum(W.*W,2);
wc=reshape(wc,[60 43]);
imagesc(wc)
colormap(gray)
