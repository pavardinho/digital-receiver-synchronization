clear
close all
clc

%% Étape 1 : Premier affichage
signal_path = fullfile(fileparts(mfilename('fullpath')), 'signal_recu.mat');
if ~isfile(signal_path)
    error('MissingSignalFile:Place signal_recu.mat beside tests.m before running this analysis script.');
end
data = load(signal_path, 'signal_recu');
signal = data.signal_recu;

figure;
plot(signal);
grid on;

%% Étape 2 : Affichage réel/imaginaire
% ici preière visualisation di signal, on observe qu'il était sur
% échantilloné, et on visiualise les constellations aussi pour constater à
% quel signalon a affaire
signal_recu = signal(:).';

N = 1000; % On ne prend que les 2000 premiers points
figure;
plot(real(signal_recu(1:N)), imag(signal_recu(1:N)), '.-');
grid on; axis equal;
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
title(['Zoom sur les ' num2str(N) ' premiers points']);

M = N; % On regarde les 500 premiers échantillons
t = 1:M;
figure;
subplot(2,1,1);
plot(t, real(signal_recu(1:M)));
title('Partie Réelle (I) en fonction du temps'); grid on;
subplot(2,1,2);
plot(t, imag(signal_recu(1:M)));
title('Partie Imaginaire (Q) en fonction du temps'); grid on;

% Paramètres du segment
start_idx = 7000;
end_idx = 8000;
segment = signal_recu(start_idx:end_idx);

% Affichage de la constellation (segment 10k-11k)
figure;
plot(real(segment), imag(segment), '.-', 'MarkerSize', 5);
grid on; axis equal;
xlabel('Partie Réelle (I)');
ylabel('Partie Imaginaire (Q)');
title(['Constellation : échantillons ' num2str(start_idx) ' à ' num2str(end_idx)]);

% Comparaison temporelle
figure;
subplot(2,1,1);
plot(start_idx:end_idx, real(segment));
ylabel('Réel (I)'); grid on;
title('Évolution temporelle (I et Q)');
subplot(2,1,2);
plot(start_idx:end_idx, imag(segment));
ylabel('Imaginaire (Q)'); xlabel('Échantillons'); grid on;

%% Étape 3 : Extraction de la raie d'horloge
% On cherche la fréquence symbole fs = 1/Ts
z = abs(signal_recu).^2; % Passage au carré pour faire ressortir le rythme
[pxx, f] = pwelch(z - mean(z), hanning(4096), 2048, 4096, 1, 'centered');

figure('Name', 'Détermination de Ts');
plot(f, 10*log10(pxx)); 
grid on;
title('Mise en évidence de la raie d''horloge');
xlabel('Fréquence normalisée (f/Fe)');
ylabel('Densité Spectrale de Puissance (dB/Hz)');

%% Étape 4 : Analyse du spectre et identification du roll-off (alpha)
% On cherche à identifier les paramètres du filtre de mise en forme (SRRC)
[pxx, f] = pwelch(signal_recu, hanning(2048), 1024, 2048, 1, 'centered');

figure('Name', 'Analyse du Filtre');
plot(f, 10*log10(pxx/max(pxx)), 'LineWidth', 1.5);
grid on;
title('Spectre du signal reçu (Zoom sur la bande utile)');
xlabel('Fréquence normalisée (f/Fe)');
ylabel('Amplitude normalisée (dB)');
xlim([-0.15 0.15]); % Focus sur le lobe principal

%% Étape 5 : Synthèse du filtre adapté h(t) et affichage
% On génère la réponse impulsionnelle du filtre de réception
Ts = 10;
alpha = 0.5;
span = 20; % On définit une durée de 20 symboles pour la troncature

% Génération du filtre SRRC (Square Root Raised Cosine)
h = rcosdesign(alpha, span, Ts, 'sqrt');

% Affichage de la réponse impulsionnelle
t_filtre = (-(span*Ts)/2 : (span*Ts)/2) / Ts;
figure('Name', 'Réponse impulsionnelle du filtre');
stem(t_filtre, h, 'filled', 'MarkerSize', 2);
grid on;
title(['Filtre adapté SRRC (\alpha = ' num2str(alpha) ', T_s = ' num2str(Ts) ')']);
xlabel('Temps normalisé (t / T_s)');
ylabel('Amplitude h(t)');

%% Étape 6 : Validation de l'alignement filtre/signal
% Calcul de la réponse fréquentielle du filtre h
[H_freq, f_h] = freqz(h, 1, 4096, 'whole', 1); 
f_h = f_h - 0.5; 
H_freq = fftshift(H_freq); 

% Calcul de la DSP du signal reçu
[pxx, f] = pwelch(signal_recu, hanning(4096), 2048, 4096, 1, 'centered');

figure('Name', 'Gabarit du Filtre');
plot(f, 10*log10(pxx/max(pxx)), 'b'); 
hold on;
% On descend un peu le tracé du filtre (-13dB) pour faciliter la lecture visuelle
plot(f_h, 20*log10(abs(H_freq)/max(abs(H_freq)))-13, 'r--', 'LineWidth', 2);
grid on; 
legend('DSP Signal Reçu', 'Réponse du Filtre Adapté');
xlim([-0.2 0.2]); 
title('Superposition du gabarit du filtre sur le spectre du signal');
xlabel('Fréquence normalisée'); ylabel('Magnitude (dB)');

%% Étape 7 : Filtrage adapté et recalage temporel
% Application du filtre SRRC
signal_filtre = filter(h, 1, signal_recu);

% Calcul et compensation du retard de groupe (Group Delay)
retard = (span * Ts) / 2;
signal_propre = signal_filtre(retard + 1 : end);

% Visualisation de l'effet du filtrage (Zoom sur 500 points)
n_zoom = 500; 
indices = 10000 : 10000 + n_zoom;

figure('Name', 'Effet du Filtrage Adapté');
subplot(2,1,1);
plot(indices, real(signal_recu(indices)), 'Color', [0.8 0.8 0.8]); hold on;
plot(indices, real(signal_propre(indices)), 'LineWidth', 1.5, 'Color', [0 0.447 0.741]);
title('Composante I : Avant vs Après filtrage');
legend('Brut (Bruité)', 'Filtré'); grid on;

subplot(2,1,2);
plot(indices, imag(signal_recu(indices)), 'Color', [0.8 0.8 0.8]); hold on;
plot(indices, imag(signal_propre(indices)), 'LineWidth', 1.5, 'Color', [0.85 0.32 0.09]);
title('Composante Q : Avant vs Après filtrage');
legend('Brut (Bruité)', 'Filtré'); grid on;
xlabel('Échantillons');

%% Étape 8 : Observation de la constellation filtrée, plus de bruit
% On vérifie la propreté du signal avant échantillonnage
figure('Name', 'Analyse du Plan Complexe après Filtrage');

% On affiche un segment restreint pour observer les transitions
segment_obs = signal_propre(10000:12000);
plot(real(segment_obs), imag(segment_obs), 'b-', 'LineWidth', 0.5); 
hold on;
plot(real(segment_obs), imag(segment_obs), 'r.', 'MarkerSize', 5); % Points de transition

axis equal; grid on;
xlabel('I (In-Phase)'); ylabel('Q (Quadrature)');
title('Constellation après Filtrage Adapté (Transitions visibles)');

%% Étape 9 : Comparaison des Densités Spectrales de Puissance (DSP)
% On vérifie l'efficacité du filtrage dans le domaine fréquentiel
[pxx_brut, f] = pwelch(signal_recu, hanning(2048), 1024, 2048, 1, 'centered');
[pxx_filt, ~] = pwelch(signal_propre, hanning(2048), 1024, 2048, 1, 'centered');

figure('Name', 'Efficacité du Filtrage Adapté');
semilogy(f, pxx_brut, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Signal Brut (Bruité)'); 
hold on;
semilogy(f, pxx_filt, 'LineWidth', 1.5, 'DisplayName', 'Signal Filtré (SRRC)');
grid on;
legend;
xlim([-0.3 0.3]); % Focus sur la zone utile et les premières lobes secondaires
xlabel('Fréquence normalisée (f/Fe)');
ylabel('Densité Spectrale de Puissance (V^2/Hz)');
title('Comparaison de la DSP avant et après filtrage adapté');

%% Étape 10 : Recherche de l'instant d'échantillonnage optimal
% On teste les 10 positions de phase possibles (Ts = 10)
Ts = 10;
energies_branches = zeros(1, Ts);

for k = 1:Ts
    % Sous-échantillonnage (décimation) pour chaque phase k
    flux_test = signal_propre(k : Ts : end);
    % Calcul de l'énergie moyenne de la branche
    energies_branches(k) = mean(abs(flux_test).^2);
end

% Identification du maximum
[~, k_opt] = max(energies_branches);

% Affichage du résultat de la synchronisation
figure('Name', 'Optimisation du Rythme');
stem(1:Ts, energies_branches, 'filled', 'MarkerFaceColor', [0 0.5 0]);
grid on;
title(['Synchronisation temporelle : k_{opt} = ' num2str(k_opt)]);
xlabel('Indice d''échantillonnage (k \in [1, Ts])');
ylabel('Énergie moyenne du flux');

%% Étape 11 : Décimation à l'instant optimal k_opt
% Passage à 1 échantillon par symbole (Symbole-rate sampling)
signal_echantillonne = signal_propre(k_opt : Ts : end);

% Visualisation de la constellation échantillonnée
figure('Name', 'Constellation Décimée');
% On affiche un nombre restreint de points pour voir la netteté
plot(real(signal_echantillonne(1000:2000)), imag(signal_echantillonne(1000:2000)), 'r.');
axis equal; grid on;
xlabel('I (In-Phase)'); ylabel('Q (Quadrature)');
title(['Constellation après décimation (k_{opt} = ' num2str(k_opt) ')']);

%% Étape 12 : Synchronisation de phase par boucle PLL
L_total = length(signal_echantillonne);
signal_pll = zeros(1, L_total);
theta_hat = zeros(1, L_total);
phi = 0; integrateur = 0;
G1 = 0.05; G2 = G1/50; % Gains de la boucle

for n = 1:L_total
    % Rotation inverse
    signal_pll(n) = signal_echantillonne(n) * exp(-1j * phi);
    
    % Estimation de l'erreur (Decision-Directed)
    I_dec = sign(real(signal_pll(n)));
    Q_dec = sign(imag(signal_pll(n)));
    erreur = imag(signal_pll(n) * (I_dec - 1j*Q_dec));
    
    % Mise à jour du filtre de boucle (PI)
    integrateur = integrateur + G2 * erreur;
    phi = phi + G1 * erreur + integrateur;
    theta_hat(n) = phi;
end

% Affichage de la convergence
figure('Name', 'Performance de la PLL');
subplot(1,2,1);
plot(theta_hat, 'LineWidth', 1.5); grid on;
title('Évolution de la phase estimée \phi');
xlabel('Symboles'); ylabel('Phase (rad)');

subplot(1,2,2);
plot(real(signal_pll(100:end)), imag(signal_pll(100:end)), '.');
axis equal; grid on; hold on;
line([-2 2], [0 0], 'Color', 'k', 'LineStyle', '--');
line([0 0], [-2 2], 'Color', 'k', 'LineStyle', '--');
title('Constellation QPSK synchronisée');


%% Étape 13 : Détection de la zone utile (Seuillage)
% L'objectif est d'aligner le premier échantillon sur le premier pixel de l'image
puissance_inst = abs(signal_echantillonne).^2;
% Lissage pour obtenir l'enveloppe du signal
puissance_lissee = filter(ones(1, 50)/50, 1, puissance_inst); 

% Détermination du seuil de détection
seuil_detection = 0.01 * max(puissance_lissee); % Valeur optimisée pour le tigre
indices_actifs = find(puissance_lissee > seuil_detection);

% Découpage du signal
idx_debut = indices_actifs(1);
idx_fin = indices_actifs(end);
signal_utile = signal_echantillonne(idx_debut : idx_fin);

% Affichage de la détection pour le rapport
figure('Name', 'Détection de Trame');
plot(puissance_lissee, 'LineWidth', 1.5); hold on;
yline(seuil_detection, 'r--', 'Seuil de détection (1%)', 'LineWidth', 2);
title(['Extraction de la zone utile : ' num2str(length(signal_utile)) ' symboles']);
xlabel('Indices des symboles'); ylabel('Puissance lissée');
grid on;

%% Étape 14 : Synchronisation de phase sur la zone utile
% Application de la PLL sur le segment de données identifié (signal_utile)
L_u = length(signal_utile);
signal_sync_final = zeros(1, L_u);
phi = 0; integrateur = 0;
G1 = 0.05; G2 = G1/50; 

for n = 1:L_u
    % Correction de phase
    signal_sync_final(n) = signal_utile(n) * exp(-1j * phi);
    
    % Boucle de décision (QPSK)
    I_dec = sign(real(signal_sync_final(n)));
    Q_dec = sign(imag(signal_sync_final(n)));
    erreur = imag(signal_sync_final(n) * (I_dec - 1j*Q_dec));
    
    % Mise à jour de la boucle
    integrateur = integrateur + G2 * erreur;
    phi = phi + G1 * erreur + integrateur;
end

% Visualisation finale avant conversion en bits
figure('Name', 'Signal prêt pour le décodage');
plot(real(signal_sync_final), imag(signal_sync_final), '.', 'MarkerSize', 1);
axis equal; grid on;
title(['Constellation finale synchronisée (' num2str(L_u) ' symboles)']);
