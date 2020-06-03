%% Modelfase 

% Kies in het menu het signaal. De output geeft het verloop van het
% inputsignaal en het outputsignaal over de tijd.

clear, clc, clf, close all

% Input
    t_max=20;         % feel free to choose a different end time of your simulation
    t=[0:1e-3:t_max]; % do NOT change the time-step 1e-3

% Eigen signaal

        % Invoeren condities TRIANGLE
        % Stap 1: Voer de maximale angle in {maximaal 34 graden, enkel positief}
        % Stap 2: Voer de gewenste hoeksnelheid in

        maxAngle_invoer = 15;
        Hoeksnelheid = 60;                      %graden per seconde
        snelheid = Hoeksnelheid*1e-3;           %graden per millisec

        % Invoeren condities STEP
        % Stap 1: Voer het gewenste aantal pulsen per seconde in (integer)
        % Stap 2: Voer de gewenste pulsduur in

        pulsen_per_seconde = 1;
        pulsduur = .2; %sec
        frequentie_puls = length(t)/((length(t)*1e-3)*pulsen_per_seconde);
        pulslengte = pulsduur/1e-3; %vectoreenheden

        % Signaal
        F(:,1) = t;
        F(:,2) = 0;
        list = {'Triangle','Steps','Cancel'};
        [idx,tf1] = listdlg("ListString",list);
        switch idx
            case 1
            Ft = 0;
            a = Hoeksnelheid*1e-3;       % snelheid
            maxAngle = maxAngle_invoer;
            peakToPeak = calculatePeakToPeak(a,maxAngle);     %tijd tussen pieken [*10 s]
            peaks = [0.5:1:1/peakToPeak-0.5];
            peaks = peaks * peakToPeak;
            for i = 1:peaks(1)*(length(t)-1)
                F(i,2) = Ft;
                Ft = Ft + a;
            end
            for k = 1:length(peaks)-1
                if mod(k,2) ~= 0
                    for i = uint16(peaks(k)*(length(t)-1)):uint16(peaks(k+1)*(length(t)-1))
                        F(i,2) = Ft;
                        Ft = Ft - a;
                    end
                else
                    for i = uint16(peaks(k)*(length(t)-1)):uint16(peaks(k+1)*(length(t)-1))
                        F(i,2) = Ft;
                        Ft = Ft + a;
                    end
                end
            end
            for i = uint16(peaks(end)*(length(t)-1)):uint16((length(t)))
                F(i,2) = Ft;
                if Ft > 0
                    Ft = Ft - a;
                else
                    Ft = Ft + a;
                end
            end
            case 2
                possibleAnglesMinus = [-34:1:-5];
                possibleAnglesPlus = [5:1:34];
                possibleAngles = [possibleAnglesMinus possibleAnglesPlus];
                Ft = 0;
                s = RandStream("dsfmt19937",'Seed',20);
                for i = 1:length(t)
                    if abs(mod(i,frequentie_puls)) <= 0.9
                        angle = randsample(s,possibleAngles,1);
                        F(i-pulslengte:i,2) = angle;
                    end
                end
                AllAngles = nonzeros(F(:,2));
                uniqueAngles = unique(AllAngles, 'stable');


            case 3
                return
        end

% Analyse aangereikt model
    % Genereer laser stimulus en output EOG en Gyroscoop
    phi_head=zeros(1,length(t));
    phi_laser_space = F; 
    phi_head_space=[t;phi_head]';
    sim('BMTM8_vExp2_prot2');
    EOG = EOG.Data';
    Gyr = Gyr.Data';

% EOG_filtering
    % Detrend haalt het baseline-artefact uit het EOG-signaal
    EOG_filter_detrend = detrend (EOG);

    % Low-pass filter zorgt voor verwijdering hoog frequente ruis
    % Parameters
    ts              = 1e-3;
    fs              = 1/ts;
    f_Nyquist       = fs/2;
    n               = 2; %2e orde filter voor nu
    fkantel         = 12.5; %hz

    % Maak een filter aan met butter
    wn              = fkantel/f_Nyquist; %genormaliseerde afkapfrequentie
    [b,a]           = butter(n,wn,'low'); %lowpass filter, verander low in high voor highpass

    EOG_filter = filtfilt(b,a,EOG_filter_detrend);

    % Kalibratiematrix [hellingsgetal, startgetal]
    B = [5.3901e-6 -6.3992e-6];

    % Output
    Return_angle_matrix = EOG_filter./(B(1,1))-(B(1,1)/B(1,2));

% Snelheidsvectoren berekenen
    speed_eyes = diff(Return_angle_matrix(70:end))./diff(t(70:end));
    speed_signal = diff(phi_laser_space(70:end,2))./diff(phi_laser_space(70:end,1));
    t_a = t(71:end);

    % Berekenen van minima en maxima in speed_eyes-array
    [maxtab, mintab] = peakdet(speed_eyes, 70, t_a);
    goodmax = [];
    goodmin = [];

    % Deze loop zorgt ervoor dat de toppen geplot kunnen worden in subplot
    % (2,1,2), ook indien de matrix leeg is
    if ~isempty(maxtab)&& ~isempty(mintab)
    [maxtab_filtered, ~] = find(maxtab(:,2)>100);
    goodmax = maxtab(maxtab_filtered,:);
    [mintab_filtered, ~] = find(mintab(:,2)<-100);
    goodmin = mintab(mintab_filtered,:);
    end

    % Defenitie van aantal saccades bij lege matrix, en bij gevulde matrix
    if isempty(goodmax)&& isempty(goodmin)
    n_saccades = 0;
    else
        n_saccades = length(goodmax')+length(goodmin');
    end

% Plotten van resultaten
    figure(1)
    subplot (2,2,1), hold on
    plot (phi_laser_space(:,1), phi_laser_space(:,2))
    plot (t, Return_angle_matrix) 
    legend ('Inputsignaal', 'Outputsignaal')
    title ({['Inputsignaal en outputsignaal met snelheid ',num2str(Hoeksnelheid),' graden/s over tijd']});
    ylabel('Hoek [graden]');
    ylim ([-50 50])
    xlabel('Tijd [seconden]')
    hold off

    subplot (2,2,2), hold on
        if isempty(maxtab)& isempty(mintab)
            plot(t_a,speed_eyes)
            plot(t_a,speed_signal)
        else
            plot (t_a,speed_eyes), hold on
            plot (goodmax(:,1), goodmax(:,2), 'r*', goodmin(:,1), goodmin(:,2), 'g*')
            plot(t_a,speed_signal), hold off
        end

    title ({'Hoeksnelheden ogen en signaal',[num2str(n_saccades),' saccades']})
    legend ('Hoeksnelheid ogen', 'Saccades naar rechts', 'Saccades naar links', 'Hoeksnelheid signaal');
    ylabel('Hoeksnelheid [graden/seconde]')
    xlim ([0 5])
    ylim([ -1000 1000]);
    xlabel('Tijd [seconden]')
    hold off
%%
% Analyse eigen model
    % Model en parameters gebaseerd op McSpadden (1998)
    % teta = hoek in graden
    % Differentiaalvergelijking: J_g teta'' + B_g teta' + k_g teta = dF[Newton]
    % gt = grams tension
    % 1 gt = 1 gram x 980 cm/s^2 = 980 dynes = 980 *10^5 Newton

% McSpadden parameters
    J_g = 6e-5;             % gt*s^2/degree (Moment of inertia)
    B_g = 0.0158;           % gt*s/degree   (Friction)       
    K_g = 0.79;             % gt/degree     (Stiffness) 
    M = 0.748;              % gram

%Mechanisch model ogen
%     wn = sqrt(K_g/J_g);                 % Natural frequency
    wn = 10;
%     mu = B_g/(2*sqrt(K_g*J_g));         % Damping ratio
    mu = 1.5;
    wdamped = wn*sqrt(1-(mu^2));          % Damped natural frequency
    k_1 = 1/K_g;                        % System gain eyes

% Westheimer parameters
     wn = 120;
     mu = 5;

% Eigenschappen transferfunctie ogen
    Transferfunction_eyes = tf([k_1],[1/wn (2*mu)/wn 1]);
    figure()
    title ('Overdracht ogen')
    subplot (1,3,1), hold on
    rlocus(Transferfunction_eyes)
    subplot (1,3,2), hold on
    bode(Transferfunction_eyes)
    subplot (1,3,3), hold on
    nyquist (Transferfunction_eyes)
    hold off

% Omzetten naar constanten model
    Time_delay = 0.00;

% Hersenen
% Beschrijving geven van:

% 1. Saccade magnitude/amplitude
    % saccade_magnitude = 20;             %degrees

% 2. Time to peak overshoot (geeft indicatie van duur saccade) (Enderle)
    % Uit Enderle (experimenteel)
        % saccade_duur_pos = 0.00181*saccade_magnitude_pos+0.025;     % [sec]
        % saccade_duur_neg = -0.00243*saccade_magnitude_pos+0.020;    % [sec]

    % Uit Westheimer (constante) (Zie Enderle hfst 3.2.1)
        %Tp = pi/(wn*sqrt(1-mu^2));
        Tp = 0.037;

% 3. Peak velocity
    % Parameters alfa & beta uit Enderle voor fit v_max in positieve (naar
    % rechts) en negatieve (naar links) richting
        % alfa_pos = 825;             % steady-state peak velocity-saccade magnitude curve
        % alfa_neg = 637;             % steady-state peak velocity-saccade magnitude curve
        % beta_pos = 9.3;             % the “time constant” for the peak velocity-saccade magnitude curve
        % beta_neg = 6.9;             % the “time constant” for the peak velocity-saccade magnitude curve
        % 
        % v_max_pos = alfa_pos*(1-exp(-saccade_magnitude_pos/beta_pos));    % Peak velocity uit Enderle (1.1) voor positieve saccade (maximaal 400 degree/sec)
        % v_max_neg = alfa_neg*(1-exp(-saccade_magnitude_neg/beta_neg));    % Peak velocity uit Enderle (1.1) voor negatatieve saccade (maximaal 400 degree/sec)

% Uit Westheimer (constante) (Zie Enderle hfst 3.2.2)
    wd = pi/Tp;
    Tmv = (1/wd)*atan(sqrt(1-mu^2)/mu);         % Time at peak velocity

% Latente periode na een saccade (zou onafhankelijk van saccade magnitude
    % zijn, wordt niet verder op ingegaan in Enderle)
    % t_latent = 0.17

% Gain signaal vanuit Brein
    Gain_brain = 20;

% Simulink Simulatie
    sim('testmodel_3');

    % Berekenen van minima en maxima in speed_eyes-array
    [maxtab, mintab] = peakdet(Speed_eyes_simulink, 70, t);
    goodmax = [];
    goodmin = [];

    % Deze loop zorgt ervoor dat de toppen geplot kunnen worden in subplot
    % (2,1,2), ook indien de matrix leeg is
    if ~isempty(maxtab)&& ~isempty(mintab)
    [maxtab_filtered, ~] = find(maxtab(:,2)>100);
    goodmax = maxtab(maxtab_filtered,:);
    [mintab_filtered, ~] = find(mintab(:,2)<-100);
    goodmin = mintab(mintab_filtered,:);
    end

    % Defenitie van aantal saccades bij lege matrix, en bij gevulde matrix
    if isempty(goodmax)&& isempty(goodmin)
    n_saccades = 0;
    else
        n_saccades = length(goodmax')+length(goodmin');
    end

    % Plotten van resultaten
    figure(1)
    subplot (2,2,3), hold on
    plot (Angle_laser_simulink)
    plot (Angle_eyes_simulink) 
    legend ('Inputsignaal', 'Outputsignaal')
    title ({['Inputsignaal en outputsignaal met snelheid ',num2str(Hoeksnelheid),' graden/s over tijd']});
    ylabel('Hoek [graden]');
    ylim ([-50 50])
    xlabel('Tijd [seconden]')
    hold off

    subplot (2,2,4), hold on
        if isempty(maxtab)& isempty(mintab)
            plot(t,Speed_eyes_simulink)
            plot(t,Speed_laser_simulink)
        else
            plot (t,Speed_eyes_simulink), hold on
            plot (goodmax(:,1), goodmax(:,2), 'r*', goodmin(:,1), goodmin(:,2), 'g*')
            plot(t,Speed_laser_simulink), hold off
        end

    title ({'Hoeksnelheden ogen en signaal',[num2str(n_saccades),' saccades']})
    legend ('Hoeksnelheid ogen', 'Saccades naar rechts', 'Saccades naar links', 'Hoeksnelheid signaal');
    ylabel('Hoeksnelheid [graden/seconde]')
    xlim ([0 5])
    ylim([ -1000 1000]);
    xlabel('Tijd [seconden]')
    hold off

%% Gebruikte functies binnen script

function peakToPeak = calculatePeakToPeak(a,maxAngle)
totalAngle = 2*maxAngle;
speed = a * 10^4;                   
peakToPeak = totalAngle/(2*speed);
end

function [maxtab, mintab]=peakdet(v, delta, x)
%PEAKDET Detect peaks in a vector
% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.

maxtab = [];
mintab = [];

v = v(:); % Just in case this wasn't a proper vector

if nargin < 3
  x = (1:length(v))';
else 
  x = x(:);
  if length(v)~= length(x)
    error('Input vectors v and x must have same length');
  end
end
  
if (length(delta(:)))>1
  error('Input argument DELTA must be a scalar');
end

if delta <= 0
  error('Input argument DELTA must be positive');
end

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;

lookformax = 1;

for i=1:length(v)
  this = v(i);
  if this > mx, mx = this; mxpos = x(i); end
  if this < mn, mn = this; mnpos = x(i); end
  
  if lookformax
    if this < mx-delta
      maxtab = [maxtab ; mxpos mx];
      mn = this; mnpos = x(i);
      lookformax = 0;
    end  
  else
    if this > mn+delta
      mintab = [mintab ; mnpos mn];
      mx = this; mxpos = x(i);
      lookformax = 1;
    end
  end
end
end



















