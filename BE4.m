clear
close all

partie1 = false;
partie2 = true;
figures = false;

% PARTIE 1
if partie1
    % Question 1 : Application des masques dérivés les directions X et Y
    I = imread('Images_HOG_1\car1.bmp');
    Ix = getXGradient(I);
    Iy = getYGradient(I);
    magnitude = getMagnitude(I);
    orientation = getOrientation(I);

    % Question 2 : Colorisation en fonction de l'orientation
    coloredOrientation = getOrientationColored(orientation);

    % Question 3 : Creation des cellules -> createCells()

    % Question 4 : Calcul et affichage des HOG
    HOG(magnitude, orientation, figures);
end

% Affichage des images
if partie1 && figures
    figure();
    subplot(3,3,1); imshow(I); title('I : image initiale');
    subplot(3,3,2); imshow(Ix,[]); title('Ix : dérivée en x');
    subplot(3,3,3); imshow(Iy,[]); title('Iy : dérivée en y');
    subplot(3,3,4); imshow(magnitude,[]); title('Magnitude du gradient');
    subplot(3,3,5); imshow(orientation,[]); title('Orientation du gradient');
    subplot(3,3,6); imshow(coloredOrientation,[]); title('Orientation du gradient colorée');
    subplot(3,3,7); createCells(magnitude);
end

% PARTIE 2
if partie2
    % Paramètres des blocs pour la normalisation
    blockHeight = 2;
    blockWidth = 2;

    I = imread('Images_HOG_2\hog_similar.bmp');
    imageGauche = I(150-127:150,1:63);
    imageDroite = I(150-127:150,90:153);
    disp("____ Image : hog_similar.bmp ____");
    calculSimilarite(imageGauche, imageDroite, blockHeight, blockWidth, figures);

    I = imread('Images_HOG_2\hog_different.bmp');
    imageGauche = I(150-127:150,1:64);
    imageDroite = I(150-127:150,91:154);
    disp("____ Image : hog_different.bmp ____");
    calculSimilarite(imageGauche, imageDroite, blockHeight, blockWidth, figures);
end

% Fonctions
function m = getMagnitude(I)
    % Calcul de l'amplitude à partir des composantes en X et Y du gradient
    Ix= getXGradient(I);
    Iy= getYGradient(I);
    m = sqrt(Ix.*Ix + Iy.*Iy);
end

function o = getOrientation(I)
    % Calcul de l'orientation à partir des composantes en X et Y du gradient
    Ix= getXGradient(I);
    Iy= getYGradient(I);
    o = atan2(Iy,Ix);
end

function Ix = getXGradient(I)
    % Convolution de chaque ligne de I par Dx
    Dx = [-1 0 1];
    [m,n]=size(I);
    for i = 1:m
       Ix(i,:) = conv(I(i,:),Dx);
    end
    % Suppression du premier et du dernier élément des colonnes de Ix pour
    % retrouver la taille initiale de l'image
    [m,n]=size(Ix);
    Ix = Ix(:,2:n-1);
end

function Iy = getYGradient(I)
    % Convolution de chaque colonne de I par Dy
    Dy = [-1 0 1]';
    [m,n]=size(I);
    for i = 1:n
       Iy(:,i) = conv(I(:,i),Dy);
    end
    % Suppression du premier et du dernier élément des lignes de Iy pour
    % retrouver la taille initiale de l'image
    [m,n]=size(Iy);
    Iy = Iy(2:m-1,:);
end

function oc = getOrientationColored(orientation)
    [m,n] = size(orientation);
    % Matrice RGB
    oc = zeros(m,n,3);
    % Coloration en Red et Green des pixels d'orientation positive donc
    % entre 0 et 180°
    oc(:,:,1) = (orientation>=0);
    oc(:,:,2) = (orientation>=0);
    % Coloration en Rouge et Green des pixels d'orientation négative donc
    % entre 0° et -180° i.e. entre 180° et 360°
    oc(:,:,3) = (orientation<0);
end

function [cells,cellHeight,cellWidth] = createCells(I)
    % Récupération des paramètres entrés par l'utilisateur
    prompt = {'Cell height:','Cell width:'};
    input = inputdlg(prompt);
    
    cellHeight = str2num(input{1});
    cellWidth = str2num(input{2});
    
    % Reprise de l'image initiale
    cells = imshow(I,[]);
    [m,n]=size(I);
    
    % Tracé des lignes verticales
    for i = 1:floor(n/cellWidth)
        l = line([i*cellWidth i*cellWidth],[1 m]);
        l.Color = "red";
    end
    
    % Tracé des lignes horizontales
    for i = 1:floor(m/cellHeight)
        l = line([1 n],[i*cellHeight i*cellHeight]);
        l.Color = "red";
    end
    
    title('Cellules');
end

function y = HOG(magnitude, orientation, figures)
    % Lecture des paramètres d'entrée
    prompt = {'Cell height:','Cell width:','Number of bins:'};
    input = inputdlg(prompt);
    
    cellHeight = str2num(input{1});
    cellWidth = str2num(input{2});
    nb_bins = str2num(input{3});
    
    % Gradient signé en degrés
    orientation = orientation*180/pi + (orientation<0)*360;
    
    if figures
        figure();
    end
    
    [H,L]=size(orientation);
    
    % Matrice de stockage des histogrammes des cellules : la trosième 
    % valeur de l'histogramme de la cellule en haut à gauche sera stockée
    % en y(1,1,3)
    y = zeros(round(H/cellHeight),round(L/cellWidth),nb_bins);
    
    % Parcours des cellules
    for i = 1:round(H/cellHeight)
        for j = 1:round(L/cellWidth)
            % Ordonnée et abscisse du coin supérieur gauche dans l'image et
            % longeur et hauteur de la cellule
            cellLimits = [(i-1)*cellHeight+1 (j-1)*cellWidth+1 cellWidth cellHeight];
            
            % Si aucun contact avec les bords de l'image
            if(cellLimits(1)+cellLimits(4) < H && cellLimits(2)+cellLimits(3)<L)
                hist = weightedhist(orientation(cellLimits(1):cellLimits(1)+cellLimits(4),cellLimits(2):cellLimits(2)+cellLimits(3)), magnitude(cellLimits(1):cellLimits(1)+cellLimits(4),cellLimits(2):cellLimits(2)+cellLimits(3)), nb_bins);
            
            % Si contact uniquement avec la droite de l'image
            elseif(cellLimits(1)+cellLimits(4) < H)
                hist = weightedhist(orientation(cellLimits(1):cellLimits(1)+cellLimits(4),cellLimits(2):L),magnitude(cellLimits(1):cellLimits(1)+cellLimits(4),cellLimits(2):L),nb_bins);
            
            % Si contact uniquement avec le bas de l'image
            elseif(cellLimits(2)+cellLimits(3)<L)
                hist = weightedhist(orientation(cellLimits(1):H,cellLimits(2):cellLimits(2)+cellLimits(3)),magnitude(cellLimits(1):H,cellLimits(2):cellLimits(2)+cellLimits(3)),nb_bins);
            
            % Si contact avec le bas et la droite de l'image (cellule dans
            % le coin inférieur droit)
            else
                hist = weightedhist(orientation(cellLimits(1):H,cellLimits(2):L),magnitude(cellLimits(1):H,cellLimits(2):L),nb_bins);
            end
            
            if figures
                % Position de l'histogramme sur la figure, coordonnées du
                % coin inférieur gauche, coordonnées entre 0 et 1
                position = [(j-1)/round(L/cellWidth) 1-i/round(H/cellHeight) 1/round(L/cellWidth) 1/round(H/cellHeight)];
                subplot('Position',position);
                x = linspace(0,360,nb_bins);
                bar(x,hist,1,'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5]);
                set(gca,'visible','off');
            end
            
            % Stockage des nb_bins composantes de l'histogramme dans y
            for a = 1:nb_bins
                y(i,j,a) = hist(a,1);
            end
        end
    end
end

function hist = weightedhist(values, weight, nb_bins)
    % Largeur des bandes de valeurs
    pas = 360/nb_bins;
    % Passage des valeurs nulles à une valeur de 1 pour éviter erreur lors
    % de la division
    values = values + (values==0)*1;
    hist = zeros(nb_bins);
    [h,l] = size(values);
    
    % Parcours de l'ensemble des valeurs
    for i=1:h
        for j=1:l
            % Indice dans l'histogramme donné par le quotient de la
            % division euclidienne de la valeur du pixel par le pas
            % (arrondi supérieur car indices commencent à 1 et non 0)
            index = ceil(values(i,j)/pas);
            hist(index) = hist(index) + weight(i,j);
        end
    end
end

function s = cosineSimilarity(HoG1, HoG2)
    [m,n,r]=size(HoG1);
    num = 0;
    norme1 = 0;
    norme2 = 0;
    
    % Parcours de l'ensemble des cellules
    for i=1:m
        for j=1:n
            % Parcours des "couches" de la cellule (k ième couche = k ième
            % valeur dans l'histogramme)
            for k=1:r
                % Ajout de la contribution de la k ième couche de la
                % cellule (i,j) au produit scalaire
                num = num + HoG1(i,j,k)*HoG2(i,j,k);
                % Ajout de la contribution de la k ième couche de la
                % cellule (i,j) aux normes des 2 HoG
                norme1 = norme1 + HoG1(i,j,k)*HoG1(i,j,k);
                norme2 = norme2 + HoG2(i,j,k)*HoG2(i,j,k);
            end
        end
    end
    % Calcul de la similarité
    s = num/(sqrt(norme1)*sqrt(norme2));
end

function s = euclideanSimilarity(HoG1, HoG2)
    [m,n,r]=size(HoG1);
    s = 0;
    % Parcours de l'ensemble des cellules
    for i=1:m
        for j=1:n
            % Parcours des "couches" de la cellule
            for k=1:r
                s = s + (HoG1(i,j,k)-HoG2(i,j,k))^2;
            end
        end
    end
    % Calcul de la similarité
    s = sqrt(s);
end

function y = RHOGnormalisationL2(HoG,bh,bw)
    % Paramètres d'entrée pour R-HoG
    blockHeight = bh; % in cells
    blockWidth = bw; % in cells
    e = 0.5; % constante faible (valeur sans importance)
    
    [m,n,r]=size(HoG);
    % Parcours de l'ensemble des blocs
    % Non chevauchement des blocs avec :
        %for i=1:m-blockHeight:blockHeight
        %    for j=1:n-blockWidth:blockWidth
        % ...
    for i=1:m-blockHeight
        for j=1:n-blockWidth
            % Bloc de coin supérieur gauche la cellule (i,j)
            
            % Matrice non-normalisée des histogrammes du bloc, contient k
            % couches (k valeurs dans l'histogramme)
            v = HoG(i:i+blockHeight-1,j:j+blockWidth-1,:);
            normev2 = 0;
            
            % Calcul de la norme 2 de v
            % Parcours des cellules du bloc
            for a=1:blockHeight
                for b=1:blockWidth
                    % Parcours des k valeurs de l'histogramme de la cellule
                    % (a,b)
                    for k=1:r
                        normev2 = normev2 + v(a,b,k)^2;
                    end
                end
            end
            normev2 = sqrt(normev2);
            % Application du facteur de normalisation
            v = v/(sqrt(e + normev2^2));
            % Remplacement des valeurs non-normalisées par celles
            % normalisées dans la matrice d'entrée
            HoG(i:i+blockHeight-1,j:j+blockWidth-1,:) = v;
        end
    end
    y = HoG;
end

function y = RHOGnormalisationL1(HoG,bh,bw)
    % Paramètres d'entrée pour R-HoG
    blockHeight = bh; % in cells
    blockWidth = bw; % in cells
    e = 0.5; % small constant
    
    [m,n,r]=size(HoG);
    % Parcours de l'ensemble des blocs
    for i=1:m-blockHeight
        for j=1:n-blockWidth
            % Matrice non-normalisée des histogrammes du bloc
            v = HoG(i:i+blockHeight-1,j:j+blockWidth-1,:);
            
            % Calcul de la norme 1 de v
            % Parcours des cellules du bloc
            normev1 = 0;
            for a=1:blockHeight
                for b=1:blockWidth
                    for k=1:r
                        normev1 = normev1 + abs(v(a,b,k));
                    end
                end
            end
            normev1 = sqrt(normev1);
            % Application du facteur de normalisation
            v = v/(e + normev1);
            % Remplacement des valeurs non-normalisées
            HoG(i:i+blockHeight-1,j:j+blockWidth-1,:) = v;
        end
    end
    y = HoG;
end

function y = RHOGnormalisationL1sqrt(HoG,bh,bw)
    % Paramètres d'entrée pour R-HoG
    blockHeight = bh; % in cells
    blockWidth = bw; % in cells
    e = 0.5; % small constant
    
    [m,n,r]=size(HoG);
    % Parcours de l'ensemble des blocs
    for i=1:m-blockHeight
        for j=1:n-blockWidth
            % Matrice non-normalisée des histogrammes du bloc
            v = HoG(i:i+blockHeight-1,j:j+blockWidth-1,:);
            
            % Calcul de la norme 1 de v
            % Parcours des cellules du bloc
            normev1 = 0;
            for a=1:blockHeight
                for b=1:blockWidth
                    for k=1:r
                        normev1 = normev1 + abs(v(a,b,k));
                    end
                end
            end
            normev1 = sqrt(normev1);
            % Application du facteur de normalisation
            v = v/(e + normev1);
            % Remplacement des valeurs non-normalisées
            HoG(i:i+blockHeight-1,j:j+blockWidth-1,:) = sqrt(v);
        end
    end
    y = HoG;
end

function calculSimilarite(imageGauche, imageDroite, blockHeight, blockWidth, figures)
    % Calcul des HoG pour les deux images
    HoGGauche = HOG(getMagnitude(imageGauche), getOrientation(imageGauche), figures);
    HoGDroite = HOG(getMagnitude(imageDroite), getOrientation(imageDroite), figures);

    % Application de la normalisation L2
    HoGGauchenormL2 = RHOGnormalisationL2(HoGGauche,blockHeight,blockWidth);
    HoGDroitenormL2 = RHOGnormalisationL2(HoGDroite,blockHeight,blockWidth);

    % Application de la normalisation L1
    HoGGauchenormL1 = RHOGnormalisationL1(HoGGauche,blockHeight,blockWidth);
    HoGDroitenormL1 = RHOGnormalisationL1(HoGDroite,blockHeight,blockWidth);

    % Application de la normalisation L1-sqrt
    HoGGauchenormL1sqrt = RHOGnormalisationL1sqrt(HoGGauche,blockHeight,blockWidth);
    HoGDroitenormL1sqrt = RHOGnormalisationL1sqrt(HoGDroite,blockHeight,blockWidth);

    % Calcul et affichage des similarités avec et sans normalisation
    similiarity = cosineSimilarity(HoGGauche,HoGDroite)
    similiaritynormL2 = cosineSimilarity(HoGGauchenormL2,HoGDroitenormL2)
    similiaritynormL1 = cosineSimilarity(HoGGauchenormL1,HoGDroitenormL1)
    similiaritynormL1sqrt = cosineSimilarity(HoGGauchenormL1sqrt,HoGDroitenormL1sqrt)
    similiarityEucli = euclideanSimilarity(HoGGauche,HoGDroite)
end
