%% WSOLA implementation
clear all; close all; clc;


%% Read input signal

%[input_signal, Fs] = audioread('Thieves_who_rob_friends_deserve_jail.wav');
[input_signal, Fs] = audioread('myspeech.wav');

%% Setup

L = Fs/1000 * 20;       % 20 ms
Nfft = L;
S = L/2;                % 50% overlap (10 ms)
win = hanning (L, 'periodic');

delta = round(Fs/1000 * 5);    % 5 ms

TSR = 1.6;

Sout = S;
Sin = round(Sout/TSR);

%% Overlap-add with time-scaling with WSOLA

% COLA_ratio = sum( win.*win ) / S;
% win = win/sqrt(COLA_ratio);

pin=0; pout=0; nseg = 1; 

inlen = length(input_signal);
outlen = ceil(TSR*inlen+delta);
output_signal = zeros(outlen,1);

synthesis_frame = zeros(L,1);
deltas = [];

% il primo frame viene inserito in posizione pin+Sin,
% quindi ne cerco uno simile a quello che inizia a pin+Sin
output_signal(pout+1:pout+L) = input_signal(pin+1:pin+L);
pref = pin + Sout; 
pin = pin + Sin;
pout = pout + Sout;


while ( (pref + L) < inlen ) && ( (pin+L+delta) < inlen )
    
    reference_frame = input_signal(pref+1:pref+L);

    % check min(pin+1-delta,1)
    analysis_frame = input_signal(pin+1-delta:pin+L+delta); %.* win;
    
    [ xc, lags ] = xcorr(analysis_frame, reference_frame, 2*delta);
    % devo trovare il massimo di xcorr nella zona -delta,+delta rispetto a
    % pin, cio? tra 0 e 2*delta di analysis frame, valori che si trovano in
    % xcorr a partire dalla posizione delta
    aligned = 2*delta+1;
    xc_delta = xc(aligned:aligned+2*delta);
    [ ~, i ] = max(abs(xc_delta));
 
    % normalized cross-correlation
    xc_rms = rms(buffer(analysis_frame(1:L+2*delta),L,L-1,'nodelay'));
    [ ~, i ] = max(abs(xc_delta./xc_rms'));
    
    idx = i-1;
    deltas = [ deltas idx ];
    
    synthesis_frame = analysis_frame(idx+1:idx+L);

	output_signal(pout+1:pout+L) = ...
        output_signal(pout+1:pout+L) + synthesis_frame .* win;
   
    %fprintf('nseg: %d, pref: %d, Sout: %d - pin: %d, Sin: %d, delta: %d - pok: %d\n', nseg, pref, Sout, pin, Sin, delta, pin - delta + idx);

    %{
    if ( nseg > 50 && nseg < 100 )

    fprintf('nseg: %d, pref: %d, Sout: %d - pin: %d, Sin: %d, delta: %d - pok: %d\n', nseg, pref, Sout, pin, Sin, delta, pin - delta + idx);
 
    m1 = pref+1:pref+L;
    m2 = pin+1-delta:pin+L+delta;
    m3 = pin+1-delta:pin+1-delta+2*delta-1;
    m4 = pin - delta + idx + 1: pin - delta + idx + L;
    n  = min(pref+1-L/2, pin+1-delta) : pref+L+L/2;
    
    subplot(4,1,1);
    plot(n,input_signal(n), m1,reference_frame);
    legend('Input','Reference');
    axis tight;
    
    subplot(4,1,2);
    plot(n,input_signal(n), m2,analysis_frame, m3, analysis_frame(1:2*delta));  
    legend('Input','Analysis', 'Delta');
    axis tight;
    
    subplot(4,1,3);
    plot(n,input_signal(n), m4,synthesis_frame);  
    legend('Input','Match');
    axis tight;
    
    n = pout+1-L/2:pout+L+L/2;
    k = pout+1:pout+L;
    subplot(4,1,4);
    plot(n, output_signal(n), k, output_signal(k));
    axis tight;
    
    fprintf('  input: %d-%d, reference: %d:%d, analysis: %d:%d\n', pref+1-L/2, pref+L+L/2, pref+1, pref+L, pin+1-delta,pin+L+delta);

    pause;
    
    end
    %}

    % reference frame starts at Sout after the beginning of the
    % newly added synthesis frame (pin+1-delta+idx)
    pref = pin - delta + idx + Sout;
    pin = pin + Sin;
    pout = pout + Sout;

    nseg = nseg + 1;
    
end;

n_seg = nseg - 1;

    
%% Listening

soundsc(output_signal(1:outlen),Fs);
%soundsc(output_signal(1:outlen),Fs*TSR);

