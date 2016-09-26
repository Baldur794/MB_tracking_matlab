idx_seq = 100; % frame index-----------------------------------------
fore_grnd_img = load_img(idx_seq); % creates a temporary img
fore_grnd_img = fore_grnd_img(544:578,81:160);
figure(7), imagesc(fore_grnd_img(:,:)), colormap(gray); colorbar;
xlabel('Lateral (mm)'); ylabel('Axial (mm)'); % title('Micro-Bubble image');
set(gca,'Xtick',linspace(0,80,9)); set(gca, 'XTickLabel',linspace(0,0.8,9));
set(gca,'Ytick',linspace(0,30,4)); set(gca, 'YTickLabel',linspace(0,0.3,4));
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])

figure(8)
plot(fore_grnd_img(20,:), 'color','black');
xlim([1 80]); ylim([0 700]);
xlabel('Lateral (mm)'); ylabel('Intensity');
set(gca,'Xtick',linspace(0,80,9)); set(gca, 'XTickLabel',linspace(0,0.8,9));

figure(9)
plot(fore_grnd_img(:,46), 'color','black')
xlim([1 35]); ylim([0 700]);
xlabel('Axial (mm)'); ylabel('Intensity');
set(gca,'Xtick',linspace(0,30,4)); set(gca, 'XTickLabel',linspace(0,0.3,4));

%%
r_index = linspace(-3,3,100);
r_gaus = normpdf(r_index,0,1);

y_l = zeros(1,200);
y_r = zeros(1,200);
y_l(1:100) = r_gaus;
y_r(101:200) = r_gaus;

figure(10); clf
plot(y_l(1:100), 'color','k', 'LineWidth', 2)
xlim([1 200]);
hold on
plot(101:200, y_r(101:200), 'color','k', 'LineWidth', 2)
%xlim([1 35]); ylim([0 700]);
xlabel('Distance'); ylabel('Intensity');
set(gca,'Xtick',linspace(0,200,3)); set(gca, 'XTickLabel',' ' );
xlim([1 200]);
set(gca, 'YTickLabel',' ' );
%%
r_index = linspace(-3,3,100);
r_gaus = normpdf(r_index,0,1);

y_sum = zeros(1,200);
y_sum(21:120) = r_gaus;
y_sum(81:180) = y_sum(81:180) + r_gaus;

y_l = zeros(1,200);
y_l(21:120) = r_gaus;

y_r = zeros(1,200);
y_r(81:180) = r_gaus;


figure(11); clf
plot(21:180, y_sum(21:180),'k','LineWidth', 2)
%plot(21:120,y_l(21:120),'-b','LineWidth', 2)
hold on
%plot(81:180,y_r(81:180),'-g','LineWidth', 2)
%xlim([1 35]); ylim([0 700]);
xlabel('Distance'); ylabel('Intensity');
set(gca,'Xtick',linspace(0,200,3)); set(gca, 'XTickLabel',' ' );
set(gca, 'YTickLabel',' ' );
%%
y_sum = zeros(1,200);
y_sum(31:130) = r_gaus;
y_sum(71:170) = y_sum(71:170) + r_gaus;

y_l = zeros(1,200);
y_l(31:130) = r_gaus;

y_r = zeros(1,200);
y_r(71:170) = r_gaus;


figure(12); clf;
plot(31:170, y_sum(31:170),'k','LineWidth', 2)
%plot(31:130,y_l(31:130),'-b','LineWidth', 2)
hold on
%plot(71:170,y_r(71:170),'-g','LineWidth', 2)
%xlim([1 35]); ylim([0 700]);
xlabel('Distance'); ylabel('Intensity');
set(gca,'Xtick',linspace(0,200,3)); set(gca, 'XTickLabel',' ' );
set(gca, 'YTickLabel',' ' );
ylim([0 0.4246]);

%%

fs=100; %sampling frequency
sigma=0.1;
t=-0.5:1/fs:0.5; %time base

variance=sigma^2;
x=1/(sqrt(2*pi*variance))*(exp(-t.^2/(2*variance)));
figure(1); clf
plot(t,[x(21:end) zeros(1,20)],'r');
hold on
plot(t,[zeros(1,20) x(1:81)],'b');
plot(t,[zeros(1,5) x(1:96)],'b');
plot(t,[x(26:end) zeros(1,25)],'r');
legend('MB_{v1} at t_1 & t_2', 'MB_{v2} at t_1 & t_2');
%set(gca,'Xtick',linspace(0,0.5,20)); 
set(gca, 'XTickLabel',' ' );
set(gca,'Ytick',linspace(0,4,11));
set(gca, 'YTickLabel',linspace(0,1,11));
xlabel('Distance');
ylabel('Amplitude (Normalized)');
 
L=length(x);
NFFT = 1024;
X = fftshift(fft(x,NFFT));
Pxx=X.*conj(X)/(NFFT*NFFT); %computing power with proper scaling
f = fs*(-NFFT/2:NFFT/2-1)/NFFT; %Frequency Vector
 
figure(2); clf
plot(f,abs(X)/(L),'k');
set(gca,'Xtick',linspace(0,0,1)); 
set(gca, 'XTickLabel','0' );
xlabel('Frequency')
ylabel('Magnitude |X(f)| (Normalized)');
xlim([-10 10])

%%
fs=100; %sampling frequency
%clf(3);
%clf(4);
for i = 0.5:1:4.5
    sigma=0.1;
    t=-0.5:1/fs:0.5; %time base

    variance=sigma^2;
    x=1/(sqrt(2*pi*variance))*(exp(-(t*i).^2/(2*variance)));
    figure(3);
    hold on
    plot(t,x);
    set(gca, 'XTickLabel',' ' );
    set(gca,'Ytick',linspace(0,4,11));
    set(gca, 'YTickLabel',linspace(0,1,11));
    xlabel('Time');
    ylabel('Amplitude (Normalized)');
    xlim([-0.2 0.2])

    L=length(x);
    NFFT = 1024;
    X = fftshift(fft(x,NFFT));
    Pxx=X.*conj(X)/(NFFT*NFFT); %computing power with proper scaling
    f = fs*(-NFFT/2:NFFT/2-1)/NFFT; %Frequency Vector

    figure(4)
    hold on
    plot(f,abs(X)/(L)*1/max(abs(X)));
    set(gca,'Xtick',linspace(0,0,1)); 
    set(gca, 'XTickLabel','0' );
    set(gca,'Ytick',linspace(0,0.01,11));
    set(gca, 'YTickLabel',linspace(0,1,11));
    xlabel('Frequency')
    ylabel('Magnitude |X(f)| (Normalized)');
    xlim([-10 10])
end

figure(3)
legend('v_{MB1} = 1','v_{MB2} = 2','v_{MB3} = 3','v_{MB4} = 4','v_{MB5} = 5');
figure(4)
legend('v_{MB1} = 1','v_{MB2} = 2','v_{MB3} = 3','v_{MB4} = 4','v_{MB5} = 5');

%%
figure(1);
s = scatter(10,10,100);
s.MarkerFaceColor = 'flat';
s.MarkerFaceColor = 'k';
s.MarkerEdgeColor = 'k';

xlim([0 30]); ylim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','on');
set(gca,'YGrid','on');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%set(gca,'Xtick',linspace(0,300,4)); 
set(gca, 'XTickLabel',' ');
%set(gca,'Ytick',linspace(0,900,10)); 
set(gca, 'YTickLabel', ' ');
xlabel('Lateral Distance'); ylabel('Axial Distance');
rectangle('position',[5 5 10 10],'EdgeColor','k','LineWidth', 2);

figure(2);
s = scatter(10,10,100);
s.MarkerFaceColor = 'flat';
s.MarkerFaceColor = 'k';
s.MarkerEdgeColor = 'k';

xlim([0 30]); ylim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','on');
set(gca,'YGrid','on');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%set(gca,'Xtick',linspace(0,300,4)); 
set(gca, 'XTickLabel',' ');
%set(gca,'Ytick',linspace(0,900,10)); 
set(gca, 'YTickLabel', ' ');
xlabel('Lateral Distance'); ylabel('Axial Distance');
rectangle('position',[5 5 10 10],'EdgeColor','k','LineWidth', 2);



%%
figure(1); clf(1);
s = scatter(10,10,100);
s.MarkerFaceColor = 'flat';
s.MarkerFaceColor = 'g';
s.MarkerEdgeColor = 'g';
hold on
c = scatter(10,10,10000,'LineWidth',1.5);
c.MarkerEdgeColor = 'k';

xlim([0 30]); ylim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','off');
set(gca,'YGrid','off');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%set(gca,'Xtick',linspace(0,300,4)); 
set(gca, 'XTickLabel',' ');
%set(gca,'Ytick',linspace(0,900,10)); 
set(gca, 'YTickLabel', ' ');
xlabel('Lateral'); ylabel('Axial');
a = quiver(10,10,4.5,4.5,'k');
a.MaxHeadSize = 1;
a.LineWidth = 1.5
text(12,11, '$d_{max}$','Interpreter','latex','Color','k','FontSize',14);
%rectangle('position',[5 5 10 10],'EdgeColor','k','LineWidth', 2);

%%
figure(1); clf(1);
s1 = scatter(10,10,100,'LineWidth',2);
% s1.MarkerFaceColor = 'flat';
% s1.MarkerFaceColor = 'g';
s1.MarkerEdgeColor = 'g';
hold on

s2 = scatter(14,10,100);
s2.MarkerFaceColor = 'flat';
s2.MarkerFaceColor = 'g';
s2.MarkerEdgeColor = 'g';

c = scatter(10,10,10000,'LineWidth',1.5);
c.MarkerEdgeColor = 'k';

xlim([0 30]); ylim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','off');
set(gca,'YGrid','off');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%set(gca,'Xtick',linspace(0,300,4)); 
set(gca, 'XTickLabel',' ');
%set(gca,'Ytick',linspace(0,900,10)); 
set(gca, 'YTickLabel', ' ');
xlabel('Lateral'); ylabel('Axial');
%a = quiver(10,10,4.5,4.5,'k');
%a.MaxHeadSize = 1;
%a.LineWidth = 1.5
%text(12,11, '$d_{max}$','Interpreter','latex','Color','k','FontSize',14);
%rectangle('position',[5 5 10 10],'EdgeColor','k','LineWidth', 2);


%%
figure(1); clf(1);
s1 = scatter(10,10,100);
s1.MarkerFaceColor = 'flat';
s1.MarkerFaceColor = 'g';
s1.MarkerEdgeColor = 'g';
hold on
s2 = scatter(6,9,100,'LineWidth',2);
%s2.MarkerFaceColor = 'flat';
%s2.MarkerFaceColor = 'g';
s2.MarkerEdgeColor = 'g';
c = scatter(14,11,10000,'LineWidth',1.5);
c.MarkerEdgeColor = 'k';

xlim([0 30]); ylim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','off');
set(gca,'YGrid','off');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%set(gca,'Xtick',linspace(0,300,4)); 
set(gca, 'XTickLabel',' ');
%set(gca,'Ytick',linspace(0,900,10)); 
set(gca, 'YTickLabel', ' ');
xlabel('Lateral'); ylabel('Axial');
a1 = quiver(14,11,4.5,4.5,'k');
a1.MaxHeadSize = 1;
a1.LineWidth = 1.5
text(16,12, '$d_{max}$','Interpreter','latex','Color','k','FontSize',14);
a2 = quiver(10,10,4.5,1.3,'k');
a2.MaxHeadSize = 1;
a2.LineWidth = 1.5
text(12,9, '$\hat{d}$','Interpreter','latex','Color','k','FontSize',14);

%rectangle('position',[5 5 10 10],'EdgeColor','k','LineWidth', 2);
%%
figure(1); clf(1);
s1 = scatter(10,10,100,'LineWidth',2);
% s1.MarkerFaceColor = 'flat';
% s1.MarkerFaceColor = 'g';
s1.MarkerEdgeColor = 'g';
hold on
s2 = scatter(6,9,100,'LineWidth',2);
%s2.MarkerFaceColor = 'flat';
%s2.MarkerFaceColor = 'g';
s2.MarkerEdgeColor = 'g';

s3 = scatter(14,10,100);
s3.MarkerFaceColor = 'flat';
s3.MarkerFaceColor = 'g';
s3.MarkerEdgeColor = 'g';

c = scatter(14,11,2000,'LineWidth',1.5);
c.MarkerEdgeColor = 'k';

xlim([0 30]); ylim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','off');
set(gca,'YGrid','off');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%set(gca,'Xtick',linspace(0,300,4)); 
set(gca, 'XTickLabel',' ');
%set(gca,'Ytick',linspace(0,900,10)); 
set(gca, 'YTickLabel', ' ');
xlabel('Lateral'); ylabel('Axial');
% a1 = quiver(14,11,4.5,4.5,'k');
% a1.MaxHeadSize = 1;
% a1.LineWidth = 1.5
% text(16,12, '$d_{max}$','Interpreter','latex','Color','k','FontSize',14);
% a2 = quiver(10,10,4.5,1.3,'k');
% a2.MaxHeadSize = 1;
% a2.LineWidth = 1.5
% text(12,9, '$\hat{d}$','Interpreter','latex','Color','k','FontSize',14);

%rectangle('position',[5 5 10 10],'EdgeColor','k','LineWidth', 2);


%%
figure(1); clf(1);
s1 = scatter(10,10,100);
s1.MarkerFaceColor = 'flat';
s1.MarkerFaceColor = 'k';
s1.MarkerEdgeColor = 'k';
hold on
c1 = scatter(10,10,10000);
c1.MarkerEdgeColor = 'k';

%s2 = scatter(22,18,100);
s2.MarkerFaceColor = 'flat';
s2.MarkerFaceColor = 'k';
s2.MarkerEdgeColor = 'k';
hold on
c2 = scatter(22,18,10000);
c2.MarkerEdgeColor = 'k';

xlim([0 30]); ylim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','off');
set(gca,'YGrid','off');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%set(gca,'Xtick',linspace(0,300,4)); 
set(gca, 'XTickLabel',' ');
%set(gca,'Ytick',linspace(0,900,10)); 
set(gca, 'YTickLabel', ' ');
xlabel('Lateral'); ylabel('Axial');
a = quiver(10,10,4.5,4.5,'k');
a.MaxHeadSize = 1;
a.LineWidth = 0.8
text(12,11, '$d_{max}$','Interpreter','latex','Color','k','FontSize',14);
%rectangle('position',[5 5 10 10],'EdgeColor','k','LineWidth', 2);


%%

fs=100; %sampling frequency
sigma=0.1;
t=-0.5:1/fs:0.5; %time base

variance=sigma^2;
x=1/(sqrt(2*pi*variance))*(exp(-(t).^2/(2*variance)));
figure(1); clf
plot(t,x,'k', 'LineWidth', 2,'Color','k');
hold on
plot(t(20:end-20),x(20:end-20),'k', 'LineWidth', 2,'Color','g');
plot([0,0], [0,0.2],'-b', 'LineWidth', 2);
plot([t(39),t(39)], [0,max(x(:))-0.1],'--r', 'LineWidth', 2);
plot([t(63),t(63)], [0,max(x(:))-0.1],'--r', 'LineWidth', 2);%legend('MB');
%set(gca,'Xtick',linspace(0,0.5,20)); 
set(gca, 'XTickLabel',' ' );
set(gca,'Ytick',linspace(0,4,11));
set(gca, 'YTickLabel',linspace(0,1,11));
xlabel('Distance');
ylabel('Amplitude');
xlim([-0.4 0.4]);
%ylim([0 5]);

h1 = annotation('arrow');
h1.X = [0.518 0.518];
h1.Y = [0.22 0.16];
h1.LineWidth = 1.5;

g1 = annotation('textbox');
g1.String = 'Centroid';
g1.EdgeColor = [1 1 1]
g1.Position = [0.439 0.29 0.1 0.0];
g1.FontSize = 12;

h = annotation('doublearrow');
h.X = [0.405 0.63];
h.Y = [0.52 0.52];
h.LineWidth = 1.5;

g = annotation('textbox');
g.String = 'FWHM*';
g.EdgeColor = [1 1 1]
g.Position = [0.46 0.6 0.1 0.0];
g.FontSize = 12;

%fs=100; %sampling frequency
%sigma=0.1;
%t=-0.5:1/fs:0.5; %time base
%%
variance=sigma^2;
x=0.99/(sqrt(2*pi*variance))*(exp(-(t).^2/(2*variance)))+0.8/(sqrt(2*pi*variance))*(exp(-((t/0.7)+0.30).^2/(2*variance)));
figure(2); clf
plot(t,x,'k', 'LineWidth', 2);
hold on
% plot([t(23),t(23)], [0,max(x(:))],'--r', 'LineWidth', 2);
% plot([t(62),t(62)], [0,max(x(:))],'--r', 'LineWidth', 2);
% plot([t(43),t(43)], [0,(x(43))],'--b', 'LineWidth', 2);
%legend('MB');
%set(gca,'Xtick',linspace(0,0.5,20)); 
set(gca, 'XTickLabel',' ' );
set(gca,'Ytick',linspace(0,4,11));
set(gca, 'YTickLabel',linspace(0,1,11));
xlabel('Distance');
ylabel('Amplitude');
xlim([-0.4 0.3]);
% 
% fs=100; %sampling frequency
% sigma=0.1;
% t=-0.5:1/fs:0.5; %time base

variance=sigma^2;
x=0.99/(sqrt(2*pi*variance))*(exp(-(t).^2/(2*variance)))+0.8/(sqrt(2*pi*variance))*(exp(-((t/0.7)+0.30).^2/(2*variance)));
figure(3); clf
plot(t,x,'k', 'LineWidth', 2);
hold on
% plot([t(41),t(41)], [0,max(x(:))],'--r', 'LineWidth', 2);
% plot([t(61),t(61)], [0,max(x(:))],'--r', 'LineWidth', 2);
% plot([t(48),t(48)], [0,(x(48))],'--b', 'LineWidth', 2);
%legend('MB');
%set(gca,'Xtick',linspace(0,0.5,20)); 
set(gca, 'XTickLabel',' ' );
set(gca,'Ytick',linspace(0,4,11));
set(gca, 'YTickLabel',linspace(0,1,11));
xlabel('Distance');
ylabel('Amplitude');
xlim([-0.4 0.3]);

%%

figure(1); clf
ylim([0 30]);
xlim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','off');
set(gca,'YGrid','off');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
hold on
plot([5,10], [7, 20],'--k');
plot([8,13], [7, 20],'--k');

s1 = scatter([7],[10],100);
s2 = scatter([10],[14],100);
%s.MarkerFaceColor = 'flat';
%s.MarkerFaceColor = 'k';
s1.MarkerEdgeColor = 'r';
s2.MarkerEdgeColor = 'b';

plot([15,20], [13, 25],'-k');
plot([18,23], [13, 25],'-k');
s1 = scatter([19],[20],100);
s2 = scatter([21],[23],100);


s1.MarkerEdgeColor = 'r';
s2.MarkerEdgeColor = 'b';
s1.MarkerFaceColor = 'flat';
s2.MarkerFaceColor = 'flat';
s1.MarkerFaceColor = 'r';
s2.MarkerFaceColor = 'b';

a1 = quiver(18,10,-11.2,-7,'-k');
a1.MaxHeadSize = 0.5;
a1.LineWidth = 0.8

text(12,8, '$\bar{d}_c$','Interpreter','latex','Color','k','FontSize',14);

%set(gca,'Xtick',linspace(0,0.5,20)); 
set(gca, 'XTickLabel',' ' );
% set(gca,'Ytick',linspace(0,4,11));
set(gca, 'YTickLabel', ' ');
xlabel('Lateral Distance'); ylabel('Axial Distance');


figure(2); clf
ylim([0 30]);
xlim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','off');
set(gca,'YGrid','off');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
hold on
plot([5,10], [7, 20],'-k');
plot([8,13], [7, 20],'-k');

s1 = scatter([9],[15],100);
s2 = scatter([11],[18],100);
s3 = scatter([7],[10],100);
s4 = scatter([10],[14],100);
%s.MarkerFaceColor = 'flat';
%s.MarkerFaceColor = 'k';
s1.MarkerFaceColor = 'r';
s2.MarkerFaceColor = 'b';
s1.MarkerEdgeColor = 'r';
s2.MarkerEdgeColor = 'b';

s3.MarkerEdgeColor = 'r';
s4.MarkerEdgeColor = 'b';



a1 = quiver(18,10,-11.2,-7,'-k');
a1.MaxHeadSize = 0.5;
a1.LineWidth = 0.8

text(12,8, '$\bar{d}_c$','Interpreter','latex','Color','k','FontSize',14);

%set(gca,'Xtick',linspace(0,0.5,20)); 
set(gca, 'XTickLabel',' ' );
% set(gca,'Ytick',linspace(0,4,11));
set(gca, 'YTickLabel', ' ');
xlabel('Lateral Distance'); ylabel('Axial Distance');
%%
figure(1); clf
ylim([0 30]);
xlim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','off');
set(gca,'YGrid','off');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
hold on
plot([5,10], [7, 20],'--k');
plot([8,13], [7, 20],'--k');

s1 = scatter([7],[10],100);
s2 = scatter([10],[14],100);
%s.MarkerFaceColor = 'flat';
%s.MarkerFaceColor = 'k';
s1.MarkerEdgeColor = 'r';
s2.MarkerEdgeColor = 'b';

plot([15,20], [13, 25],'-k');
plot([18,23], [13, 25],'-k');
s1 = scatter([19],[20],100);
s2 = scatter([21],[23],100);


s1.MarkerEdgeColor = 'r';
s2.MarkerEdgeColor = 'b';
s1.MarkerFaceColor = 'flat';
s2.MarkerFaceColor = 'flat';
s1.MarkerFaceColor = 'r';
s2.MarkerFaceColor = 'b';

%a1 = quiver(18,10,-11.2,-7,'-k');
%a1.MaxHeadSize = 0.5;
%a1.LineWidth = 0.8

%text(12,8, '$\bar{d}_c$','Interpreter','latex','Color','k','FontSize',14);

%set(gca,'Xtick',linspace(0,0.5,20)); 
set(gca, 'XTickLabel',' ' );
% set(gca,'Ytick',linspace(0,4,11));
set(gca, 'YTickLabel', ' ');
xlabel('Lateral Distance'); ylabel('Axial Distance');
%%

figure(2); clf
ylim([0 30]);
xlim([0 30]);
set(gca,'Box','on');
set(gca,'XGrid','off');
set(gca,'YGrid','off');
set(gca,'GridColorMode','manual')
set(gca,'GridColor',[0 0 0])
set(gca,'GridAlpha',0.3)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
hold on
% plot([5,15], [0, 30],'-k','LineWidth',1.5);
% plot([12,22], [0, 30],'-k','LineWidth',1.5);

s1 = scatter(linspace(7,13,4),linspace(3,21,4),100);
s2 = scatter(linspace(9,18,5),linspace(1,28,5),100);
s3 = scatter(linspace(12,20,7),linspace(4,28,7),100);

s4 = scatter([7 5 9 4],[21 28 26 23],100);
%s5 = scatter([3 18 25 22 14],[10 7 11 20 19],100);

%  c1 = scatter([7],[21],7000);
%  c1.MarkerEdgeColor = 'k';
 
c2 = scatter([9],[ 9],7000);
c2.MarkerEdgeColor = 'k';
% 
% c11 = scatter([7],[21],100);
% c11.MarkerFaceColor = 'r';
% c11.MarkerEdgeColor = 'r';
% 
c22 = scatter([9],[9],100);
c22.MarkerFaceColor = 'm';
c22.MarkerEdgeColor = 'm';

%s.MarkerFaceColor = 'flat';
%s.MarkerFaceColor = 'k';
% s1.MarkerFaceColor = 'flat';
% s2.MarkerFaceColor = 'flat';
% s3.MarkerFaceColor = 'flat';
% s4.MarkerFaceColor = 'flat';
% s5.MarkerFaceColor = 'flat';
% % 
% s1.MarkerFaceColor = 'k';
% s2.MarkerFaceColor = 'k';
% s3.MarkerFaceColor = 'k';
% s4.MarkerFaceColor = 'k';
% s5.MarkerFaceColor = 'k';
% 
% s1.MarkerEdgeColor = 'k';
% s2.MarkerEdgeColor = 'k';
% s3.MarkerEdgeColor = 'k';
% s4.MarkerEdgeColor = 'k';
% s5.MarkerEdgeColor = 'k';
% 
s1.MarkerFaceColor = 'm';
s2.MarkerFaceColor = 'b';
s3.MarkerFaceColor = 'g';
s4.MarkerFaceColor = 'r';
%s5.MarkerFaceColor = 'k';

s1.MarkerEdgeColor = 'm';
s2.MarkerEdgeColor = 'b';
s3.MarkerEdgeColor = 'g';
s4.MarkerEdgeColor = 'r';
%s5.MarkerEdgeColor = 'k';


%set(gca,'Xtick',linspace(0,0.5,20)); 
set(gca, 'XTickLabel',' ' );
% set(gca,'Ytick',linspace(0,4,11));
set(gca, 'YTickLabel', ' ');
xlabel('Lateral Distance'); ylabel('Axial Distance');
%%




img_B_mode = load_img_B_mode(1000);
norm = max(abs(img_B_mode(:)));
limg=20*log10(abs((img_B_mode))/norm);
figure(1);
%h_b = imagesc(limg,[-50 0]);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
h_b = imagesc(abs(img_B_mode), [0 80]);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
colormap(gca,'gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)');
set(gca, 'DataAspectRatio',[1 3.46 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca,'Xtick',linspace(0,280,5)); set(gca, 'XTickLabel',linspace(0,12,5));
set(gca,'Ytick',linspace(0,1960,6)); set(gca, 'YTickLabel',linspace(0,25,6));

%%

img_size_contrast = [489,61];
img_size_B_mode = [1960,280];

k = 1;
i = 10000;
img_contrast = zeros(img_size_contrast(1),img_size_contrast(2));
img = load_img_contrast(i,k);
%h_b = imagesc(img);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
norm = max(abs(img(:)));
limg=20*log10(abs((img))/norm);
figure(1);
h_b = imagesc(limg,[-10 0]);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
colormap(gca,'gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)');
set(gca, 'DataAspectRatio',[1 3.46 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca,'Xtick',linspace(0,61,5)); set(gca, 'XTickLabel',linspace(0,12,5));
set(gca,'Ytick',linspace(0,489,6)); set(gca, 'YTickLabel',linspace(0,25,6));

img_size_contrast = [489,61];
img_size_B_mode = [1960,280];

k = 1;
i = 900;
img_contrast = zeros(img_size_contrast(1),img_size_contrast(2));
img = load_img_contrast(i,k);
%h_b = imagesc(img);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
norm = max(abs(img(:)));
limg=20*log10(abs((img))/norm);
figure(1);
h_b = imagesc(limg,[-12 0]);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
colormap(gca,'gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)');
set(gca, 'DataAspectRatio',[1 3.46 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca,'Xtick',linspace(0,61,5)); set(gca, 'XTickLabel',linspace(0,12,5));
set(gca,'Ytick',linspace(0,489,6)); set(gca, 'YTickLabel',linspace(0,25,6));

%%
figure(1)

k = 1;
i = 900;

img = load_img_contrast(i,k);
plot(img);
xlim([1 489]);
xlabel('Axial (mm)'); ylabel('Intensity');
set(gca,'Xtick',linspace(0,489,6)); set(gca, 'XTickLabel',linspace(0,25,6));

figure(1)

k = 1;
i = 900;

img = load_img_contrast(i,k);
plot(img);
xlim([1 489]);
xlabel('Axial (mm)'); ylabel('Intensity');
set(gca,'Xtick',linspace(0,489,6)); set(gca, 'XTickLabel',linspace(0,25,6));

%%
figure(1)

k = 1;
i = 60;
img = load_img_contrast(i,k);
hi = histogram(img(img(:)<100),92);
hi.FaceColor = 'k';
hi.FaceAlpha  = 1;
xlim([0 100]);
xlabel('Amplitude'); ylabel('Frequency');
mean(img(:))
std(img(:))
%%

hold on;
gca
plot([40 40], [0 36],'--g', 'LineWidth', 2)
plot([2 81], [18 18],'--b', 'LineWidth', 2)

%%
hold on;
viscircles([180 230],40,'EnhanceVisibility',false,'LineStyle','--','Color','g');
viscircles([130 520],40,'EnhanceVisibility',false,'LineStyle','--','Color','g');
viscircles([100 740],40,'EnhanceVisibility',false,'LineStyle','--','Color','g');
viscircles([165 310],30,'EnhanceVisibility',false,'LineStyle','--','Color','b');
viscircles([35 180],30,'EnhanceVisibility',false,'LineStyle','--','Color','r');
viscircles([280 600],30,'EnhanceVisibility',false,'LineStyle','--','Color','r');

%%
hold on
xlim([260 300]);
ylim([580 620]);
viscircles([280 600],18,'EnhanceVisibility',false,'LineStyle','--','Color','r');

%%
hold on
imagesc(zeros(999,301));
s1 = scatter([180 130 100],[230 520 740],100,'+g','LineWidth',1.5);
s2 = scatter([282 46],[65 445],100,'+r','LineWidth',1.5);

%%

hold on;
caxis([30 650]);
s1 = scatter([50],[50],100,'+r','LineWidth',1.5);
s1.MarkerFaceColor = 'flat';
s1.MarkerFaceColor = 'r';


%%
hold on
 plot([-10.7,-10.7], [0,140],'--r', 'LineWidth', 2);
 plot([10,10], [0,140],'--r', 'LineWidth', 2);
 
%%
%%
figure(8); clf
psf = scatter(50,50);
psf.SizeData = 50000;
psf.MarkerFaceColor = 'g';
psf.MarkerEdgeColor = 'none';
psf = plot(50, 50, '.g', 'MarkerSize',600);
%viscircles([50 50],50,'EnhanceVisibility',false,'LineStyle','-','Color','k');
hold on
plot(50, 50, '.b', 'MarkerSize',69)
set(gca,'Box','off');
set(gca,'Visible','off');
set(gca, 'DataAspectRatio',[1 1 1])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
set(gca,'Xtick',[]); set(gca, 'XTickLabel',linspace(0,10,6));
set(gca,'Ytick',[]); set(gca, 'YTickLabel',linspace(0,20,11));

h1 = annotation('arrow');
h1.X = [0.15 0.15];
h1.Y = [0.9 0.7];
h1.LineWidth = 1.5;

h2 = annotation('arrow');
h2.X = [0.15 0.3];
h2.Y = [0.9 0.9];
h2.LineWidth = 1.5;


g1 = text(-29,88,'Axial');
set(g1,'Rotation',90)
set(g1,'FontSize',12)
g2 = text(-20,112,'Lateral');
set(g2,'FontSize',12)
g3 = text(38,60,'Centroid');
set(g3,'FontSize',12)
xlim([-10 110]);
ylim([-10 110]);

%% AM sequence
x = -5:.01:5;
y = zeros(1,4000);
y(500:1500)=sinc(x);
y(2500:3500)=1/2*sinc(x);
figure();
plot(y(1:2000),'k');
hold on;
%plot([2000 2000], [0.03 -0.03],'k','LineWidth',1);
% xlim([1 430]); 
ylim([-1 2]);
ylim([-0.5 1]);
set(gca,'Xtick',linspace(1,4000,3)); set(gca, 'XTickLabel',[]);
%set(gca,'Ytick',linspace(-2,2,5)); 
set(gca, 'YTickLabel',[]);
% xlabel('Time'); 
ylabel('Amplitude'); % title('Micro-Bubble image');


