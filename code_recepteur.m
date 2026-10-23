clear
close all
clc

level = 1;

%% Chargement du fichier contenant le signal reçu
signal_path = fullfile(fileparts(mfilename('fullpath')), 'signal_recu.mat');
if ~isfile(signal_path)
    error('MissingSignalFile:Place signal_recu.mat beside code_recepteur.m before running this receiver.');
end
signal_data = load(signal_path, 'signal_recu');
if ~isfield(signal_data, 'signal_recu')
    error('InvalidSignalFile: signal_recu.mat must contain a variable named signal_recu.');
end
signal_recu = signal_data.signal_recu(:).';
if level > 1
  signal_recu = signal_recu(1:5:end);
end

%% Votre récepteur
% 1. Paramètres identifiés lors de tes tests
Ts = 10;                
alpha = 0.5;            
span = ceil(10/alpha);  
h = rcosdesign(alpha, span, Ts, 'sqrt'); % sqrt pour respecter le critère de Nyquist

% 2. Filtrage adapté et compensation du retard
signal_filtre = filter(h, 1, signal_recu);
retard = (span * Ts) / 2;
signal_propre = signal_filtre(retard + 1 : end);

% 3. Détection de la zone utile (Seuillage pour enlever le bruit début/fin)
puissance_lissee = filter(ones(1, 100)/100, 1, abs(signal_propre).^2);
seuil = 0.01 * max(puissance_lissee);
indices_utiles = find(puissance_lissee > seuil);
signal_utile = signal_propre(indices_utiles(1) : indices_utiles(end));

% 4. Synchronisation temporelle (recherche de l'instant optimal k_opt)
% On cherche l'instant où l'énergie est maximale sur les 10 positions possibles
energie = zeros(1, Ts);
for k = 1:Ts
    energie(k) = mean(abs(signal_utile(k:Ts:end)).^2);
end
[~, k_opt] = max(energie);
signal_echantillonne = signal_utile(k_opt : Ts : end);

% 5. Synchronisation de Phase (PLL) pour corriger la rotation
L_u = length(signal_echantillonne);
signal_sync = zeros(1, L_u);
phi = 0; 
integrateur = 0;
G1 = 0.05; % Corrige l'erreur de la phase instantanée, plus il est petit, plus la PLL est stable, mais pus c'est long
G2 = G1/50; % Estime et compense le décalage de fréquence, ajuste la vitesse de rotation

for n = 1:L_u
    signal_sync(n) = signal_echantillonne(n) * exp(-1j * phi);
    I_dec = sign(real(signal_sync(n)));
    Q_dec = sign(imag(signal_sync(n)));
    erreur = imag(signal_sync(n) * (I_dec - 1j*Q_dec));
    integrateur = integrateur + G2 * erreur;
    phi = phi + G1 * erreur + integrateur;
end

% 6. Génération de hatB
bit_LSB = real(signal_sync) > 0; 
bit_MSB = imag(signal_sync) > 0;
hatB = [bit_LSB ; bit_MSB]; % L'image est en niveaux de gris codée sur 8 bits
% par pixel. Mais comme hatB a 2 lignes (I et Q), chaque colonne
% de hatB contient 2 bits. Pour avoir 8 bits, il me faut 4 colonnes.

% Ajustement
% On a besoin de T*T pixels, et chaque pixel fait 8 bits.
% Nombre de bits nécessaires = 128 * 128 * 8 = 131072 bits.
nb_bits_voulus = 128 * 128 * 4; % car hatB a 2 lignes (2 * Ns = total bits)
if size(hatB, 2) > nb_bits_voulus
    hatB = hatB(:, 1:nb_bits_voulus); % On tronque le surplus
elseif size(hatB, 2) < nb_bits_voulus
    hatB(:, end+1 : nb_bits_voulus) = 0; % On complète avec des zéros si besoin
end

%% Décodage de source
hatMatBitImg = reshape(hatB(:),[],8);
matImg = bi2de(hatMatBitImg); 
T = 128; 
Img = reshape(matImg,T,T);

%% Affichage
figure
imagesc(Img)
colormap gray
