%% Pirma uzduotis
% Daugiasluoksnio neurono struktura:
% - Vienas iejimas x su duomenimis 0.1:1/22:1;
% - Aproksimuojama funkcija y = ((1 + 0.6*sin(2*pi*x/0.7)) + (0.3*sin(2*pi*x)))/2;
% - Paslepto sluoksnio aktyvavimo funkcija hiperbolinis tangentas arba
% - sigmoide.
% - Isejimo neurone tiesine funkcija
% - Mokymo algoritmas - atgalinio sklidimo
% - Pasleptas sluoksnis is 6 neuronu
clc
clear
close all


% Iejimo vektorius
x= 0.1:1/22:1;

d = ((1 + 0.6*sin(2*pi*x/0.7)) + (0.3*sin(2*pi*x)))/2;

% Parametrai
n = 0.15;
% Pirmas sluoksnis
w11_1 = rand(1);
w21_1 = rand(1);
w31_1 = rand(1);
w41_1 = rand(1);
w51_1 = rand(1);
w61_1 = rand(1);
b1_1 = rand(1);
b2_1 = rand(1);
b3_1 = rand(1);
b4_1 = rand(1);
b5_1 = rand(1);
b6_1 = rand(1);

% Antras sluoksnis
w11_2 = rand(1);
w12_2 = rand(1);
w13_2 = rand(1);
w14_2 = rand(1);
w15_2 = rand(1);
w16_2 = rand(1);
b1_2 = rand(1);

% Mokymosi epochos
for index=1:20000

    % Neuronų atsako ir aktyvavimo funkcijų skaičiavimai, bei svorių
    % atnaujinimas su kiekvienų x(i) pavyzdžiu.
    for i = 1:length(x)
        % Pirmo sluoksnio atsako ir aktyvavimo funkcijos
        v1_1 = w11_1*x(i) + b1_1;
        y1_1 = tanh(v1_1);
        
        v2_1 = w21_1*x(i) + b2_1;
        y2_1 = tanh(v2_1);
        
        v3_1 = w31_1*x(i) + b3_1;
        y3_1 = tanh(v3_1);

        v4_1 = w41_1*x(i) + b4_1;
        y4_1 = tanh(v4_1);
    
        v5_1 = w51_1*x(i) + b5_1;
        y5_1 = tanh(v5_1);

        % v6_1 = w61_1*x(i) + b6_1;
        % y6_1 = tanh(v6_1);

        % Antro sluoksnio atsako ir aktyvavimo funkcijos
        % v1_2 = y1_1*w11_2 + y2_1*w12_2 + y3_1*w13_2...
        %     +y4_1*w14_2 + y5_1*w15_2 + y6_1*w16_2 + b2_1;
                v1_2 = y1_1*w11_2 + y2_1*w12_2 + y3_1*w13_2...
            +y4_1*w14_2 + y5_1*w15_2 + b2_1;
        y1_2 = v1_2;

        % Klaidos skaiciavimas
        e = d(i) - y1_2;

        % Mokymas
        % Antro sluoksnio
            % delta = e*f'(v). Kadangi f(v) = v, tai f'(v) = 1.
            delta1_2 = e;   
            w11_2 = w11_2 + n*delta1_2*y1_1;
            w12_2 = w12_2 + n*delta1_2*y2_1;
            w13_2 = w13_2 + n*delta1_2*y3_1;
            w14_2 = w14_2 + n*delta1_2*y4_1;
            w15_2 = w15_2 + n*delta1_2*y5_1;
            % w16_2 = w16_2 + n*delta1_2*y6_1;
            b1_2 = b1_2 + n*delta1_2;

        % Pirmo sluoksnio
            % delta = 
            delta1_1 = (1 - tanh(v1_1)^2)*delta1_2*w11_2;
            delta2_1 = (1 - tanh(v2_1)^2)*delta1_2*w12_2;
            delta3_1 = (1 - tanh(v3_1)^2)*delta1_2*w13_2;
            delta4_1 = (1 - tanh(v4_1)^2)*delta1_2*w14_2;
            delta5_1 = (1 - tanh(v5_1)^2)*delta1_2*w15_2;
            % delta6_1 = (1 - tanh(v6_1)^2)*delta1_2*w16_2;

            w11_1 = w11_1 + n*delta1_1*x(i);
            b1_1 = b1_1 + n*delta1_1;

            w21_1 = w21_1 + n*delta2_1*x(i);
            b2_1 = b2_1 + n*delta2_1;

            w31_1 = w31_1 + n*delta3_1*x(i);
            b3_1 = b3_1 + n*delta3_1;

            w41_1 = w41_1 + n*delta4_1*x(i);
            b4_1 = b4_1 + n*delta4_1;

            w51_1 = w51_1 + n*delta5_1*x(i);
            b5_1 = b5_1 + n*delta5_1;

            % w61_1 = w61_1 + n*delta6_1*x(i);
            % b6_1 = b6_1 + n*delta6_1;
    end
end

% Funkcijos aproksimavimas su naujomis x vertėmis
apx_x = linspace(0.1,1,200);
apx_y = zeros(1,length(apx_x));
for i=1:length(apx_x)
            % Pirmo sluoksnio atsako ir aktyvavimo funkcijos
        v1_1 = w11_1*apx_x(i) + b1_1;
        y1_1 = tanh(v1_1);
        
        v2_1 = w21_1*apx_x(i) + b2_1;
        y2_1 = tanh(v2_1);
        
        v3_1 = w31_1*apx_x(i) + b3_1;
        y3_1 = tanh(v3_1);

        v4_1 = w41_1*apx_x(i) + b4_1;
        y4_1 = tanh(v4_1);
    
        v5_1 = w51_1*apx_x(i) + b5_1;
        y5_1 = tanh(v5_1);

        % v6_1 = w61_1*apx_x(i) + b6_1;
        % y6_1 = tanh(v6_1);

        % Antro sluoksnio atsako ir aktyvavimo funkcijos
        % v1_2 = y1_1*w11_2 + y2_1*w12_2 + y3_1*w13_2...
        %     +y4_1*w14_2 + y5_1*w15_2 + y6_1*w16_2 + b2_1;
                v1_2 = y1_1*w11_2 + y2_1*w12_2 + y3_1*w13_2...
            +y4_1*w14_2 + y5_1*w15_2 + b2_1;

        y1_2 = v1_2;
        apx_y(1,i) = y1_2;
end

plot(x,d,'bo', apx_x,apx_y,'r--')
title("Funkcijos aproksimacija")
ylabel('y')
xlabel('x')

%% Papildoma uzduotis
% f(x,y) = x^2 + y^2
clc
clear
close all

[x,y] = meshgrid(0.1:1/4:1,0.1:1/4:1);
z = x.^2 + y.^2;

% Parametrai
n = 0.15;
% Pirmas sluoksnis
w11_1 = rand(1);
w12_1 = rand(1);
w21_1 = rand(1);
w22_1 = rand(1);
w31_1 = rand(1);
w32_1 = rand(1);
w41_1 = rand(1);
w42_1 = rand(1);

b1_1 = rand(1);
b2_1 = rand(1);
b3_1 = rand(1);
b4_1 = rand(1);


% Antras sluoksnis
w11_2 = rand(1);
w12_2 = rand(1);
w13_2 = rand(1);
w14_2 = rand(1);
b1_2 = rand(1);

% Mokymosi epochos
for index=1:5000

    % Neuronų atsako ir aktyvavimo funkcijų skaičiavimai, bei svorių
    % atnaujinimas su kiekvienų x(i) pavyzdžiu.
    for i = 1:height(x)
        for j=1:length(x)
            % Pirmo sluoksnio atsako ir aktyvavimo funkcijos
            v1_1 = w11_1*x(i,j) + w12_1*y(i,j) + b1_1;
            y1_1 = tanh(v1_1);
            
            v2_1 = w21_1*x(i,j) + w22_1*y(i,j) + b2_1;
            y2_1 = tanh(v2_1);
            
            v3_1 = w31_1*x(i,j) + w32_1*y(i,j) + b3_1;
            y3_1 = tanh(v3_1);
    
            v4_1 = w41_1*x(i,j) + w42_1*y(i,j) + b4_1;
            y4_1 = tanh(v4_1);
    
            % Antro sluoksnio atsako ir aktyvavimo funkcijos
            v1_2 = y1_1*w11_2 + y2_1*w12_2 + y3_1*w13_2...
                +y4_1*w14_2 + b2_1;
            y1_2 = v1_2;
    
            % Klaidos skaiciavimas
            e = z(i,j) - y1_2;

        % Mokymas
        % Antro sluoksnio
            % delta = e*f'(v). Kadangi f(v) = v, tai f'(v) = 1.
            delta1_2 = e;   
            w11_2 = w11_2 + n*delta1_2*y1_1;
            w12_2 = w12_2 + n*delta1_2*y2_1;
            w13_2 = w13_2 + n*delta1_2*y3_1;
            w14_2 = w14_2 + n*delta1_2*y4_1;
            b1_2 = b1_2 + n*delta1_2;

        % Pirmo sluoksnio
            % delta = 
            delta1_1 = (1 - tanh(v1_1)^2)*delta1_2*w11_2;
            delta2_1 = (1 - tanh(v2_1)^2)*delta1_2*w12_2;
            delta3_1 = (1 - tanh(v3_1)^2)*delta1_2*w13_2;
            delta4_1 = (1 - tanh(v4_1)^2)*delta1_2*w14_2;

            w11_1 = w11_1 + n*delta1_1*x(i,j);
            w12_1 = w12_1 + n*delta1_1*y(i,j);
            b1_1 = b1_1 + n*delta1_1;

            w21_1 = w21_1 + n*delta2_1*x(i,j);
            w22_1 = w22_1 + n*delta2_1*y(i,j);
            b2_1 = b2_1 + n*delta2_1;

            w31_1 = w31_1 + n*delta3_1*x(i,j);
            w32_1 = w32_1 + n*delta3_1*y(i,j);
            b3_1 = b3_1 + n*delta3_1;

            w41_1 = w41_1 + n*delta4_1*x(i,j);
            w42_1 = w42_1 + n*delta4_1*y(i,j);
            b4_1 = b4_1 + n*delta4_1;
        end
    end
end

% Paviršiaus aproksimavimas su naujomis x ir y vertėmis
[apx_x apx_y] = meshgrid(0.1:1/22:1,0.1:1/22:1);
apx_z =zeros(height(apx_x), length(apx_x));
for i=1:height(apx_x)
    for j = 1:length(apx_x)
        % Pirmo sluoksnio atsako ir aktyvavimo funkcijos
        v1_1 = w11_1*apx_x(i,j) + w12_1*apx_y(i,j) + b1_1;
        y1_1 = tanh(v1_1);
        
        v2_1 = w21_1*apx_x(i,j) + w22_1*apx_y(i,j) + b2_1;
        y2_1 = tanh(v2_1);
        
        v3_1 = w31_1*apx_x(i,j) + w32_1*apx_y(i,j) + b3_1;
        y3_1 = tanh(v3_1);

        v4_1 = w41_1*apx_x(i,j) + w42_1*apx_y(i,j) + b4_1;
        y4_1 = tanh(v4_1);

        % Antro sluoksnio atsako ir aktyvavimo funkcijos
        v1_2 = y1_1*w11_2 + y2_1*w12_2 + y3_1*w13_2...
            +y4_1*w14_2 + b2_1;
        y1_2 = v1_2;
        apx_z(i,j) = y1_2;
    end
end

surf(x, y, z,'FaceColor','none')
hold on
surf(apx_x, apx_y, apx_z, 'FaceColor','none','LineStyle','--','EdgeColor','r')
% legend('Realios paviršiaus funkcijos taškai','Aproksimuotos paviršiaus funkcijos taškai')
title('f(x,y)=x^{2} + y^{2}')
xlabel('x')
ylabel('y')
zlabel('z')

