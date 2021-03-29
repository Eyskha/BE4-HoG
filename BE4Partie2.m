clear
close all

affichage_images = false;

% PARTIE 2
% Question 1
% Calcul et affichage des HOG de l'image hog_similar.bmp
I = imread('Images_HOG_2\hog_similar.bmp');
imageGauche = I(150-127:150,1:63);
imageDroite = I(150-127:150,90:153);
HoGGauche = HOG(getMagnitude(imageGauche), getOrientation(imageGauche));
HoGDroite = HOG(getMagnitude(imageDroite), getOrientation(imageDroite));

HoGGauchenormL2 = RHOGnormalisationL2(HoGGauche);
HoGDroitenormL2 = RHOGnormalisationL2(HoGDroite);

HoGGauchenormL1 = RHOGnormalisationL1(HoGGauche);
HoGDroitenormL1 = RHOGnormalisationL1(HoGDroite);

HoGGauchenormL1sqrt = RHOGnormalisationL1sqrt(HoGGauche);
HoGDroitenormL1sqrt = RHOGnormalisationL1sqrt(HoGDroite);

similiarityQ1 = cosineSimilarity(HoGGauche,HoGDroite)
similiarityQ1normL2 = cosineSimilarity(HoGGauchenormL2,HoGDroitenormL2)
similiarityQ1normL1 = cosineSimilarity(HoGGauchenormL1,HoGDroitenormL1)
similiarityQ1normL1sqrt = cosineSimilarity(HoGGauchenormL1sqrt,HoGDroitenormL1sqrt)
%similiarityQ1E = euclideanSimilarity(HoGGauche,HoGDroite)

I = imread('Images_HOG_2\hog_different.bmp');
imageGauche = I(150-127:150,1:64);
imageDroite = I(150-127:150,91:154);
HoGGauche = HOG(getMagnitude(imageGauche), getOrientation(imageGauche));
HoGDroite = HOG(getMagnitude(imageDroite), getOrientation(imageDroite));

HoGGauchenormL2 = RHOGnormalisationL2(HoGGauche);
HoGDroitenormL2 = RHOGnormalisationL2(HoGDroite);

HoGGauchenormL1 = RHOGnormalisationL1(HoGGauche);
HoGDroitenormL1 = RHOGnormalisationL1(HoGDroite);

HoGGauchenormL1sqrt = RHOGnormalisationL1sqrt(HoGGauche);
HoGDroitenormL1sqrt = RHOGnormalisationL1sqrt(HoGDroite);

similiarityQ2 = cosineSimilarity(HoGGauche,HoGDroite)
similiarityQ2normL2 = cosineSimilarity(HoGGauchenormL2,HoGDroitenormL2)
similiarityQ2normL1 = cosineSimilarity(HoGGauchenormL1,HoGDroitenormL1)
similiarityQ2normL1sqrt = cosineSimilarity(HoGGauchenormL1sqrt,HoGDroitenormL1sqrt)
%similiarityQ2E = euclideanSimilarity(HoGGauche,HoGDroite)

% I = imread('Images_HOG_2\hog_similar2.bmp');
% imageGauche = I(142-127:142,1:64);
% imageDroite = I(142-127:142,91:154);
% HoGGauche = HOG(getMagnitude(imageGauche), getOrientation(imageGauche));
% HoGDroite = HOG(getMagnitude(imageDroite), getOrientation(imageDroite));
% similiarityQ1b = cosineSimilarity(HoGGauche,HoGDroite)
% 
% I = imread('Images_HOG_2\hog_different2.bmp');
% size(I)
% imageGauche = I(142-127:142,1:64);
% imageDroite = I(142-127:142,91:154);
% HoGGauche = HOG(getMagnitude(imageGauche), getOrientation(imageGauche));
% HoGDroite = HOG(getMagnitude(imageDroite), getOrientation(imageDroite));
% similiarityQ2b = cosineSimilarity(HoGGauche,HoGDroite)


% Affichage des images
if affichage_images
end

% Fonctions
function m = getMagnitude(I)
    Ix= getXGradient(I);
    Iy= getYGradient(I);
    m = sqrt(Ix.*Ix + Iy.*Iy);
end

function o = getOrientation(I)
    Ix= getXGradient(I);
    Iy= getYGradient(I);
    o = atan2(Iy,Ix);
end

function Ix = getXGradient(I)
    Dx = [-1 0 1];
    [m,n]=size(I);
    for i = 1:m
       Ix(i,:) = conv(I(i,:),Dx);
    end
    [m,n]=size(Ix);
    Ix = Ix(:,2:n-1);
end

function Iy = getYGradient(I)
    Dy = [-1 0 1]';
    [m,n]=size(I);
    for i = 1:n
       Iy(:,i) = conv(I(:,i),Dy);
    end
    [m,n]=size(Iy);
    Iy = Iy(2:m-1,:);
end

function oc = getOrientationColored(orientation)
    [m,n] = size(orientation);
    oc = zeros(m,n,3);
    oc(:,:,1) = (orientation>=0);
    oc(:,:,2) = (orientation>=0);
    oc(:,:,3) = (orientation<0);
end

function [cells,cellHeight,cellWidth] = createCells(I)
    prompt = {'Cell height:','Cell width:'};
    input = inputdlg(prompt);
    
    cellHeight = str2num(input{1});
    cellWidth = str2num(input{2});
    
    cells = imshow(I,[]);
    % Lignes verticales
    for i = 1:floor(length(I)/cellWidth)
        l = line([i*cellWidth i*cellWidth],[1 height(I)]);
        l.Color = "red";
    end
    
    % Lignes horizontales
    for i = 1:floor(height(I)/cellHeight)
        l = line([1 length(I)],[i*cellHeight i*cellHeight]);
        l.Color = "red";
    end
    
    title('Cells display');
end

function y = HOG(magnitude, orientation)
    % Lecture des paramètres d'entrée
    prompt = {'Cell height:','Cell width:','Number of bins:'};
    input = inputdlg(prompt);
    
    cellHeight = str2num(input{1});
    cellWidth = str2num(input{2});
    nb_bins = str2num(input{3});
    
    % Gradient signé
    orientation = orientation*180/pi + (orientation<0)*360;
    
    % Parcours des cellules
    % figure();
    [H,L]=size(orientation); 
    y = zeros(round(H/cellHeight),round(L/cellWidth),9);
    for i = 1:round(H/cellHeight)
        for j = 1:round(L/cellWidth)
            cellLimits = [(i-1)*cellHeight+1 (j-1)*cellWidth+1 cellWidth cellHeight];
            %position = [(j-1)/round(L/cellWidth) 1-i/round(H/cellHeight) 1/round(L/cellWidth) 1/round(H/cellHeight)];
            
            if(cellLimits(1)+cellLimits(4) < H && cellLimits(2)+cellLimits(3)<L)
                hist = weightedhist(orientation(cellLimits(1):cellLimits(1)+cellLimits(4),cellLimits(2):cellLimits(2)+cellLimits(3)), magnitude(cellLimits(1):cellLimits(1)+cellLimits(4),cellLimits(2):cellLimits(2)+cellLimits(3)), nb_bins);
            elseif(cellLimits(1)+cellLimits(4) < H)
                hist = weightedhist(orientation(cellLimits(1):cellLimits(1)+cellLimits(4),cellLimits(2):L),magnitude(cellLimits(1):cellLimits(1)+cellLimits(4),cellLimits(2):L),nb_bins);
            elseif(cellLimits(2)+cellLimits(3)<L)
                hist = weightedhist(orientation(cellLimits(1):H,cellLimits(2):cellLimits(2)+cellLimits(3)),magnitude(cellLimits(1):H,cellLimits(2):cellLimits(2)+cellLimits(3)),nb_bins);
            else
                hist = weightedhist(orientation(cellLimits(1):H,cellLimits(2):L),magnitude(cellLimits(1):H,cellLimits(2):L),nb_bins);
            end
            
            %subplot('Position',position);
            %x = linspace(20,340,nb_bins);
            %bar(x,hist,2);
            %set(gca,'visible','off');
            %xlim([-10 370]);
            
            for a = 1:nb_bins
                y(i,j,a) = hist(a,1);
            end
        end
    end
end

function hist = weightedhist(values, weight, nb_bins)
    pas = 360/nb_bins;
    values = values + (values==0)*1;
    hist = zeros(nb_bins,1);
    [h,l] = size(values);
    for i=1:h
        for j=1:l
            index = ceil(values(i,j)/pas);
            hist(index,1) = hist(index,1) + weight(i,j);
        end
    end
end

function s = cosineSimilarity(HoG1, HoG2)
    [m,n,r]=size(HoG1);
    num = 0;
    norme1 = 0;
    norme2 = 0;
    for i=1:m
        for j=1:n
            for k=1:r
                num = num + HoG1(i,j,k)*HoG2(i,j,k);
                norme1 = norme1 + HoG1(i,j,k)*HoG1(i,j,k);
                norme2 = norme2 + HoG2(i,j,k)*HoG2(i,j,k);
            end
        end
    end
    s = num/(sqrt(norme1)*sqrt(norme2));
end

function s = euclideanSimilarity(HoG1, HoG2)
    [m,n,r]=size(HoG1);
    s = 0;
    for i=1:m
        for j=1:n
            for k=1:r
                s = s + (HoG1(i,j,k)-HoG2(i,j,k))^2;
            end
        end
    end
    s = sqrt(s);
end

function y = RHOGnormalisationL2(HoG)
    % Paramètres d'entrée pour R-HoG
    blockHeight = 2; % in cells
    blockWidth = 2; % in cells
    e = 0.5; % small constant
    
    [m,n,r]=size(HoG);
    for i=1:m-blockHeight
        for j=1:n-blockWidth
            % Bloc de coin supérieur gauche en (i,j)
            % Vecteur non-normalisé des histogrammes du block
            v = HoG(i:i+blockHeight-1,j:j+blockWidth-1,:);
            normev2 = 0;
            for a=1:blockHeight
                for b=1:blockWidth
                    for k=1:r
                        normev2 = normev2 + v(a,b,k)^2;
                    end
                end
            end
            normev2 = sqrt(normev2);
            v = v/(sqrt(e + normev2^2));
            HoG(i:i+blockHeight-1,j:j+blockWidth-1,:) = v;
        end
    end
    y = HoG;
end

function y = RHOGnormalisationL1(HoG)
    % Paramètres d'entrée pour R-HoG
    blockHeight = 2; % in cells
    blockWidth = 2; % in cells
    e = 0.5; % small constant
    
    [m,n,r]=size(HoG);
    for i=1:m-blockHeight
        for j=1:n-blockWidth
            % Bloc de coin supérieur gauche en (i,j)
            % Vecteur non-normalisé des histogrammes du block
            v = HoG(i:i+blockHeight-1,j:j+blockWidth-1,:);
            normev1 = 0;
            for a=1:blockHeight
                for b=1:blockWidth
                    for k=1:r
                        normev1 = normev1 + abs(v(a,b,k));
                    end
                end
            end
            normev1 = sqrt(normev1);
            v = v/(e + normev1);
            HoG(i:i+blockHeight-1,j:j+blockWidth-1,:) = v;
        end
    end
    y = HoG;
end

function y = RHOGnormalisationL1sqrt(HoG)
    % Paramètres d'entrée pour R-HoG
    blockHeight = 2; % in cells
    blockWidth = 2; % in cells
    e = 0.5; % small constant
    
    [m,n,r]=size(HoG);
    for i=1:m-blockHeight
        for j=1:n-blockWidth
            % Bloc de coin supérieur gauche en (i,j)
            % Vecteur non-normalisé des histogrammes du block
            v = HoG(i:i+blockHeight-1,j:j+blockWidth-1,:);
            normev1 = 0;
            for a=1:blockHeight
                for b=1:blockWidth
                    for k=1:r
                        normev1 = normev1 + abs(v(a,b,k));
                    end
                end
            end
            normev1 = sqrt(normev1);
            v = v/(e + normev1);
            HoG(i:i+blockHeight-1,j:j+blockWidth-1,:) = sqrt(v);
        end
    end
    y = HoG;
end