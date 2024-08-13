x = rand(1,100)*0.8;
beta = 1;
y1 = beta*x + (rand(1,length(x))-0.5)*0.2;
purity = rand(1,length(x))*0.9+0.1;
y2 = max(x)*beta -x.*beta + (rand(1,length(x))-0.5)*0.2;
xp = purity.*x;
y1p = y1.*purity;
y2p = y2.*purity;
% false example 1
xf1 = randn(1,100)*0.1 + 0.55;
yf1 = (randn(1,100))*0.1 + 0.55;
xf1p = xf1.*purity;
yf1p = yf1.*purity;
% false example 2
xf2 = [x(1:80),randn(1,20)*0.07 + 0.1];
yf2 = [y2(1:80), randn(1,20)*0.07 + 0.1];
xf2p = xf2.*purity;
yf2p = yf2.*purity;

% false example 2
xf3 = [x(1:20),randn(1,80)*0.07 + 0.1];
yf3 = [y2(1:20), randn(1,80)*0.07 + 0.1];
xf3p = xf3.*purity;
yf3p = yf3.*purity;


figure(11);
cc = [0.9923    0.8034    0.1939];
sz = 6;
f = purity > 0.35 & purity < 0.85;
subplot(3,3,1);scatter(x,y1,sz,cc,"filled","markeredgecolor","k",'LineWidth',0.2); title("W/O interception, pure");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.8]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(x',y1');
text(0.35,0.1,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar

subplot(3,3,2);scatter(xp,y1p,sz,purity,"filled","markeredgecolor","k",'LineWidth',0.2); title("W/O interception, mix");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.7]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xp',y1p');
text(0.35,0.1,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar

subplot(3,3,4);scatter(xf1,yf1,sz,cc,"filled","markeredgecolor","k",'LineWidth',0.2); title("W/O interception, pure");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.8]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xf1',yf1');
text(0.3,0.1,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar

subplot(3,3,5);scatter(xf1p,yf1p,sz,purity,"filled","markeredgecolor","k",'LineWidth',0.2); title("W/O interception, mix");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.7]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xf1p',yf1p');
text(0.3,0.1,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar




figure(12);
subplot(3,3,1);scatter(x,y2,sz,cc,"filled","markeredgecolor","k",'LineWidth',0.2); title("W/ interception, pure");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.8]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(x',y2');
text(0.3,0.7,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar

subplot(3,3,2);scatter(xp,y2p,sz,purity,"filled","markeredgecolor","k",'LineWidth',0.2); title("W/ interception, mix");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.7]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xp',y2p');
text(0.3,0.7,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar


subplot(3,3,3);scatter(xp(f),y2p(f),sz,purity(f),"filled","markeredgecolor","k",'LineWidth',0.2); title("W/ interception, mix");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.7]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xp(f)',y2p(f)');
text(0.3,0.7,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar





subplot(3,3,4);scatter(xf2,yf2,sz,cc,"filled","markeredgecolor","k",'LineWidth',0.2); title("W/ interception, pure");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.8]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xf2',yf2');
text(0.3,0.7,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar

subplot(3,3,5);scatter(xf2p,yf2p,sz,purity,"filled","markeredgecolor","k",'LineWidth',0.2); title("W/ interception, mix");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.7]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xf2p',yf2p');
text(0.3,0.7,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar

subplot(3,3,6);scatter(xf2p(f),yf2p(f),sz,purity(f),"filled","markeredgecolor","k",'LineWidth',0.2); title("W/ interception, mix");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.7]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xf2p(f)',yf2p(f)');
text(0.3,0.7,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar



subplot(3,3,7);scatter(xf3,yf3,sz,cc,"filled","markeredgecolor","k",'LineWidth',0.2); title("W/ interception, pure");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.8]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xf3',yf3');
text(0.3,0.7,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar

subplot(3,3,8);scatter(xf3p,yf3p,sz,purity,"filled","markeredgecolor","k",'LineWidth',0.2); title("W/ interception, mix");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.7]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xf3p',yf3p');
text(0.3,0.7,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar

subplot(3,3,9);scatter(xf3p(f),yf3p(f),sz,purity(f),"filled","markeredgecolor","k",'LineWidth',0.2); title("W/ interception, mix");
xlabel("\Delta Variable 1");ylabel("\Delta Variable 2");xlim([-0.05,0.7]);ylim([-0.05,0.8]);colormap(parula)
R1 = corr(xf3p(f)',yf3p(f)');
text(0.3,0.7,strcat("R=",num2str(R1)));caxis([0.2 1]);
colorbar