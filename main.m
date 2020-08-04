clc;
clear all;
close all;

[filename,filepath]=uigetfile({'*.jpg'},'Selecione o primeiro olho');
img1 = imread(strcat(filepath, filename));

[filename,filepath]=uigetfile({'*.jpg'},'Selecione o segundo olho');
img2 = imread(strcat(filepath, filename));



new1 = img1;
[a1,b1]=size(img1);
[ci1,cp1,o1]=thresh(img1,100,160); %Baugman pra localização da iris e pupila 
R1 = cp1(3) + 50; %raio maximo de analise da pupila (menor que a iris pra evitar cilios e palpebras)
Rp1 = cp1(3)+10; %maior que a pupila pra evitar que ela seja processa por acidente

new2 = img2;
[a2,b2]=size(img2);
[ci2,cp2,o2]=thresh(img2,100,160); %Baugman pra localização da iris e pupila 
R2 = cp2(3) + 50;%raio maximo de analise da pupila (menor que a iris pra evitar cilios e palpebras)
Rp2 = cp2(3)+10;%maior que a pupila pra evitar que ela seja processa por acidente

Olho1 = calEye(ci1,R1,Rp1,a1,b1,new1);
Olho2 = calEye(ci2,R2,Rp2,a2,b2,new2);

IDNT1 = mat2gray(Olho1);
IDNT1 = imbinarize(IDNT1); %binarizando para a comparacao logica entre as identidades

IDNT2 = mat2gray(Olho2);
IDNT2 = imbinarize(IDNT2);

compara = hemming(IDNT1, IDNT2); %diferenca entre as iris


if compara > 0.465
    msgbox('Os olhos são diferentes');
else
    msgbox('Os olhos são iguais');
end


function H = Hfilter( x, y, P, Q, D, n) %Funcao Notch
    H = zeros(P,Q); %Filtro vazio
    for i = 1 : P
        for j = 1 : Q
               Dk1 = sqrt( (j  - x)^2 + (i  - y)^2 ); %Primeiro ponto
               H(i,j) =  (1/(1 + (D/Dk1)^2*n )); %tende a anular ao redor dos pontos dados, em um raio D                 
        end
    end
end

function ham = hemming(orig, comp)
    [a, ~] = size(orig);
    [c, ~] = size(comp);
    HD = 0;
    VD = 0;
    if a < c
        x = a;
    else
        x = c;
    end
    y = 360;
    numPix = (x*y);
    
    for i = 1:x %compara os pixels em ordem vertcial
        for j = 1:y
            VD = VD + xor(orig(i,j),comp(i,j));
        end
    end
    
    for j = 1:y %compara os pixels em ordem horizontal
        for i = 1:x
            HD = HD + xor(orig(i,j),comp(i,j));
        end
    end
    
    ham = (HD+VD)/(numPix*2); %divide a soma da comparação de pixels pelo dobro de pixels presentes na imagem
end

function [B, hor, ver] = padd(img)
    [a,b] = size(img);
    hor = 8*ceil(a/8) - a; %quantidade de blocos 8x8
    ver = 8*ceil(b/8) - b;
    A = img; %matriz clone para operações
    % Analisa o tamanho da imagem para a realizacao de um padding, caso necessario,
    % para a divisao da imagem em blocos de 8


    A = padarray(A, [hor, 0], 0, 'pre'); %padding para blocos 8x8 perfeitos
    A = padarray(A, [0, ver], 0, 'pre');

    [x, y] = size(A);

    B = double(zeros(8,8,( (y/8).*(x/8) ))); %Matriz que armazena os blocos 8x8
    Bquant = size(B);

    XX = 0;
    YY = 0;

    %Montagem dos Blocos 8x8 a partir da imagem com padding
    for bloco = 1 : Bquant(3)

        for i = 1 : 8

            for j = 1 : 8
                B(i,j, bloco) = A(i + XX,j + YY); 
            end
        end
        
        YY = YY + 8;
        
        if YY >= y
            YY = 0;
            XX = XX + 8;
            if XX >= x                
                break;            
            end
        end
    end
end

%Remonta os blocos na ordem da imagem original
function C = deblock(B, orig)
    [~,~,q] = size(B);
    [x, y] = size(orig);
    C = zeros(x,y);
    XX = 0;
    YY = 0;    
    for bloco = 1 : q
        for i = 1 : 8
            for j = 1 : 8
                C(i + XX, j + YY) = B(i, j, bloco);
            end
        end
        
        YY = YY + 8;
        
        if YY >= y
            YY = 0;
            XX = XX + 8;
            if XX >= x
                break;
            end
        end
    end
end


function resultado = calEye(ci, R, Rp, a, b, new) %funçcao que tem todo o codigo dos calculos gerais
    polar = zeros(abs(R - Rp),360); %matriz com valores da iris em coordenadas polares

    for i = 1:a
        for j = 1:b
            if ((i-ci(1))^2 + (j-ci(2))^2 <= (R)^2) && (i-ci(1))^2 + (j-ci(2))^2 > Rp^2 %dentro do circulo da iris e fora da pupila           
                r = (sqrt((i-ci(1))^2 + (j-ci(2))^2) - Rp); %raio -  raio pupila
                if r < 0
                    r = ceil(r);
                else
                    r = floor(r);
                end            
                tan = (j - ci(2))/(i - ci(1)); %tangente na circunferencia
                theta = atand(tan); %arco tangente

                if theta < 0 %quadrante 2 ou 4
                    if i < ci(1) %quadrante 2
                        theta = 180 + theta;

                    else %quadrante 4
                        theta = 360 + theta;                    
                    end

                else %quadrantes 1 e 3
                    if i < ci(1) %quadrante 3, 1 eh normal                
                        theta = theta + 180;              
                    end
                end
                theta = floor(theta); %inteiro
                polar(r+1, theta+1) = new(i,j); %copia pixel pra tabela polar            
            end
        end
    end

    [a,b]=size(polar);

    for i = 1:a
        for j = 2:b
            if polar(i,j)==0
                polar(i,j) = polar(i,j-1);
            end
        end
    end

    norm = histeq(uint8(polar));

    [B, hor, ver] = padd(norm);

    [m, n, q] = size(B);

    DCT = zeros(m,n,q);
    iDCT = DCT;
    lowDCT = DCT;
    for bloco = 1:q
        DCT(:,:,bloco) = dct2(B(:,:,bloco));
    %DCT(:,:,bloco) = B(:,:,bloco);

        D = 6;
        num = 4;

        H = Hfilter(m,n, m, n, D, num);

        lowDCT(:,:,bloco) = DCT(:,:,bloco).*H;
        iDCT(:,:,bloco) = idct2(lowDCT(:,:,bloco));
        iDCT(:,:,bloco) = histeq(uint8(iDCT(:,:,bloco)));
    end

    C = deblock(iDCT, norm);
    FINAL = C;

    for i = 1 : a %tirando padding
        for j = 1 : b
            FINAL(i , j) = floor(C(i + hor, j + ver)); 
        end
    end
    resultado = FINAL;
end