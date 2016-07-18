fprf = 1000/6;
Tprf = 1/fprf;
frame_skip = 1;
frame_rate = fprf/frame_skip;
T = 1/frame_rate;
f_nyquist = frame_rate/2;
f0 = 5*10^6;
c = 1540;
%% Filter
N = 5;
N_skip = 10;
b = zeros(1,N*N_skip);
for i = N_skip:N_skip:N*N_skip
   b(i) = 1/N; 
end
b = [1 -b];

[h,w] = freqz(b,1,512);
figure(1);
freqz(b)

%% Input
dr_lateral = 10*10^(-6);
MB_size_metric = 0.5*10^(-3);
MB_size = MB_size_metric/dr_lateral;

v = 1*10^(-3);

t_pix = dr_lateral/v;
t_total = t_pix*MB_size;

r_length = ceil(t_total/T);
r_index = linspace(-3,3,r_length);

r_gaus = normpdf(r_index,0,1);

r_sin_I = sin(2*pi*2*v/c*f0*linspace(0,r_length*T,r_length));
r_sin_Q = sin(2*pi*2*v/c*f0*linspace(0,r_length*T,r_length)+pi/2);


r_in = sqrt((r_gaus.*r_sin_I).^2+(r_gaus.*r_sin_Q).^2);
r_in_I = r_gaus.*r_sin_I;
r_in_Q = r_gaus.*r_sin_Q;
% Output;
r_echo_env = abs(filter(b,1,r_in));
r_filter_I = filter(b,1,r_in_I);
r_filter_Q = filter(b,1,r_in_Q);
r_echo_IQ = sqrt((r_filter_I).^2+(r_filter_Q).^2);

figure(1); clf;
hold on
plot(r_index,r_in);
plot(r_index,r_echo_env);
legend('r_{env}(t)', 'r_{env\_echo}(t)');
set(gca,'Xtick',linspace(-3,3,7)); 
set(gca, 'XTickLabel',' ' );
set(gca,'Ytick',linspace(0,0.45,12));
set(gca, 'YTickLabel',linspace(0,1.1,12));
xlabel('Time');
ylabel('Amplitude (Normalized)');
ylim([0 0.45]);


figure(2); clf;
hold on
plot(r_index,r_in);
plot(r_index,r_echo_IQ);
legend('r_{env}(t)', 'r_{RF\_echo}(t)');
set(gca,'Xtick',linspace(-3,3,7)); 
set(gca, 'XTickLabel',' ' );
set(gca,'Ytick',linspace(0,0.45,12));
set(gca, 'YTickLabel',linspace(0,1.1,12));
xlabel('Time');
ylabel('Amplitude (Normalized)');
ylim([0 0.45]);

figure(3); clf;
hold on

plot(r_index,r_in);
plot(r_index,r_in_I);
plot(r_index,r_in_Q);
legend('r_{env}(t)', 'r_{I}(t)','r_{Q}(t)');
set(gca,'Xtick',linspace(-3,3,7)); 
set(gca, 'XTickLabel',' ' );
set(gca,'Ytick',linspace(-0.4,0.4,21));
set(gca, 'YTickLabel',linspace(-1,1,21));
xlabel('Time');
ylabel('Amplitude (Normalized)');


%% Filter f0
fsh = c/(2*v)*1/T;
%fsh = c/(2*1)*3200;
N = 1;
N_skip = 1;
b = zeros(1,N);
for i = 1:N
   b(i) = i*N_skip;
end
b = [0 -b];
b = diag(b);
f = repmat(linspace(0,10*10^6,10000)',1,size(b,2));
b = sum(exp(-j*2*pi*(f*b)/fsh)*diag([1/2 -1/2*ones(1,size(b,1)-1)]),2)';
hf = figure(2);
hl = plot(f(:,1),20*log10(abs(b)));
ha = gca;
set(hl,'color','black');
xlabel('Frequency [MHz]'); ylabel('Magnitude [dB]');
set(gca,'Xtick',linspace(0,10*10^6,6))
set(gca, 'XTickLabel',linspace(0,10,6))
ylim([-40 0]);
%set(gca,'Ytick',linspace(0,684,6))
%set(gca, 'YTickLabel',{'0', '6.8', '13.7','20.5','27.36','35.1'})
%b = 1/2*(1-z^(1*2*(1/3200)/1540))

%% Filter contrast
n_bck_grnd = 5;
n_bck_grnd_skip = 10;
b_back = zeros(1,n_bck_grnd*n_bck_grnd_skip);
for i = n_bck_grnd_skip:n_bck_grnd_skip:n_bck_grnd*n_bck_grnd_skip
   b_back(i) = 1/n_bck_grnd; 
end
b_back = [-b_back];

n_fore_grnd = 3;
b_fore = zeros(1,n_fore_grnd);
for i = 1:n_fore_grnd
   b_fore(i) = 1/n_fore_grnd; 
end
b = b_back;
b(1:size(b_fore,2)) = b_fore;


[h,w] = freqz(b,1,512);
figure(1);
freqz(b,1,512)


