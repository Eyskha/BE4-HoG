clear
close all

affichage_images = true;

% Question 1
% Application des masques dérivés les directions X et Y

I = imread('Images_HOG_1\car1.bmp');

Ix = getXGradient(I);
Iy = getYGradient(I);

magnitude = getMagnitude(I);
orientation = getOrientation(I);


% Question 2
% Colorisation en fonction de l'orientation
coloredOrientation = getOrientationColored(orientation);


% Question 3
% Creation des cellules


% Question 4
% Calcul et affichage des HOG
HOG(magnitude, orientation);


% Affichage des images
if affichage_images
    figure();
    % 3.1
    subplot(3,3,1); imshow(I); title('I : image initiale');
    subplot(3,3,2); imshow(Ix,[]); title('Ix : dérivée en x');
    subplot(3,3,3); imshow(Iy,[]); title('Iy : dérivée en y');
    subplot(3,3,4); imshow(magnitude,[]); title('Magnitude du gradient');
    subplot(3,3,5); imshow(orientation,[]); title('Orientation du gradient');
    subplot(3,3,6); imshow(coloredOrientation,[]); title('Orientation du gradient colorée');
    %subplot(3,3,7); createCells(magnitude);
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
    for i = 1:height(I)
       Ix(i,:) = conv(I(i,:),Dx);
    end
    Ix = Ix(:,2:length(Ix)-1);
end

function Iy = getYGradient(I)
    Dy = [-1 0 1]';
    for i = 1:length(I)
       Iy(:,i) = conv(I(:,i),Dy);
    end
    Iy = Iy(2:height(Iy)-1,:);
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
    figure();
    H = height(orientation);
    L = length(orientation);
    for i = 1:round(H/cellHeight)
        for j = 1:round(L/cellWidth)
            cellLimits = [(i-1)*cellHeight+1 (j-1)*cellWidth+1 cellWidth cellHeight];
            position = [(j-1)/round(L/cellWidth) 1-i/round(H/cellHeight) 1/round(L/cellWidth) 1/round(H/cellHeight)];
            
            subplot('Position',position);
            if(cellLimits(1)+cellLimits(4) < H && cellLimits(2)+cellLimits(3)<L)
                hist = histogram(orientation(cellLimits(1):cellLimits(1)+cellLimits(4),cellLimits(2):cellLimits(2)+cellLimits(3)),nb_bins);
                hist.Values
            elseif(cellLimits(1)+cellLimits(4) < H)
                histogram(orientation(cellLimits(1):cellLimits(1)+cellLimits(4),cellLimits(2):L),nb_bins);
            elseif(cellLimits(2)+cellLimits(3)<L)
                histogram(orientation(cellLimits(1):H,cellLimits(2):cellLimits(2)+cellLimits(3)),nb_bins);
            else
                histogram(orientation(cellLimits(1):H,cellLimits(2):L),nb_bins);
            end

            set(gca,'visible','off');
            xlim([-10 370]);
        end
    end
end














