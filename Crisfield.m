clear
clc
close all
format shortG
%------------------------------------------------------------------------
% Método do comprimento de arco para resolução de treliça não linear
%------------------------------------------------------------------------
% Geovane dos Santos Gomes - junho 2021
%
% Referências:
% 1-Notas de aula prof. Juan Avila - Fundamentos MecSol II - Posmec UFABC
% 2-R. Borst, M. Crisfield; J. Remmers - Nonlinear Finite Element Analysis
%   of Solids and Structures
% 3-Zheng Li - 1D Arch-length method with Modified newton raphson iteration
%------------------------------------------------------------------------
% Carregamento da solução analítica
%------------------------------------------------------------------------
analitico=load('ResultadoAnalitico.txt');
x = analitico(:,1);
y = analitico(:,2);
%------------------------------------------------------------------------
% Dados iniciais
%------------------------------------------------------------------------
l=0.1; %comprimento de arco inicial
delta_l=0.05; %incremento do comprimento de arco
L0=1; %comprimento inicial da barra da treliça
A=1; %área da seção transversal
E=1; %módulo de elasticidade
beta0=pi/3; %ângulo inicial
beta=beta0; %ângulo da barra
L=L0; %comprimento da barra
N=0; %força axial na barra
imax=500; %número máximo de iterações
tol=10e-6; %tolerância
Residual=1; %carga residual inicial (valor alto o suficiente para entrar no loop de iteração
delta_lambda_ini=0.08; %incremento do fator de carga inicial
i=1; %indicador de iterações
q=2.321; %módulo da carga
lambda=0; %fator de carga inicial
u=0; %deslocamento inicial
%------------------------------------------------------------------------
% Configuração da plotagem
%------------------------------------------------------------------------
% figure(1)
% h1=plot(x,y,'b','LineWidth',1.5);
% hold on
% grid on
% axis equal
%------------------------------------------------------------------------
% início do loop de iteração
%------------------------------------------------------------------------
while (i<=imax && abs(Residual)>=tol)
    %
    if i==1
        %----------------------------------------------------------------
        % Cálculo da força residual na configuração inicial
        %----------------------------------------------------------------
        delta_lambda(i)=delta_lambda_ini;
        F_int=2*N*sin(beta);
        F_ext=delta_lambda(i)*q;
        Residual=F_ext-F_int;
        %----------------------------------------------------------------
        % Matriz de rigidez tangente na configuração inicial
        %----------------------------------------------------------------
        k=2*N*(cos(beta))^3/(L0*cos(beta0))+2*A*E/L0*(sin(beta))^2;
        %----------------------------------------------------------------
        % Atualização do fator de carga e do deslocamento
        %----------------------------------------------------------------
        delta_u(i+1)=Residual/k;
        lambda(i+1)=lambda(i)+delta_lambda(i);
        u(i+1)=u(i)+delta_u(i+1);
        delta_lambda(i+1)=delta_lambda(i);
    else
        %----------------------------------------------------------------
        % Atualização da matriz de rigidez tangente na configuração atual
        %----------------------------------------------------------------
        L=sqrt(u(i)^2-2*L0*sin(beta0)*u(i)+L0^2);
        N=A*E/L0*(L0-L);
        beta=asin((L0*sin(beta0)-u(i))/L);
        k=2*N*(cos(beta))^3/(L0*cos(beta0))+2*A*E/L0*(sin(beta))^2;
        %----------------------------------------------------------------
        % Cálculo da força residual na configuração atual
        %----------------------------------------------------------------
        F_int=2*N*sin(beta);
        F_ext=lambda(i)*q;
        Residual=F_ext-F_int;
        %----------------------------------------------------------------
        % Atualização do fator de carga e do deslocamento
        %----------------------------------------------------------------
        var_u1=Residual./k;
        var_u2=q./k;
        a1=var_u2'*var_u2;
        a2=2*var_u2'*(delta_u(i)+var_u1);
        a3=(delta_u(i)+var_u1)'*(delta_u(i)+var_u1)-l^2+(delta_lambda(i)^2);
        %----------------------------------------------------------------
        % Resolução da equação quadrática
        % a1 x**2 + a2 x + a3
        %----------------------------------------------------------------
        if (a1==0 && a2==0)
            disp('Sem raízes reais')
        elseif (a1==0 && a2~=0)
            delta_lambda=-a3./a2;
        else
            disc=a2*a2-4*a1*a3;
            if disc<0
                disp('Sem raízes reais')
            end
            lambda1=(-a2+sqrt(disc))/(2*a1);
            lambda2=(-a2-sqrt(disc))/(2*a1);
            if lambda2>=lambda1
                var_lambda=lambda2;
            else
                var_lambda=lambda1;
            end
        end
        %
        delta_lambda(i+1)=delta_lambda(i)+var_lambda;
        delta_u(i+1)=delta_u(i)+var_u1+var_lambda*var_u2;
        %
        lambda(i+1)=0+delta_lambda(i+1);
        u(i+1)=0+delta_u(i+1);
        l=l+delta_l;
        %
        if lambda(i+1)>1;
            break
        end
        %
    end
    resposta(i,1)=u(i); %vetor deslocamento
    resposta(i,2)=lambda(i)*q-Residual; %vetor de carga
    %----------------------------------------------------------------
    % Plotagem do "path-folowing"
    %----------------------------------------------------------------
%     if i>2
%         figure(1)
%         h2=plot(resposta(i,1),resposta(i,2),'O','Color','k','Markersize',10);
%         title(['Incrementos = ',num2str(i),' ','Residual = ',num2str(Residual, '%10.5e\n')])
%         legend([h1 h2],'Solução analítica','Solução aproximada','Location','northwest')
%         xlabel('Deslocamento u')
%         ylabel('Carga \lambda*q')
%         pause(0.1)
%     end
    %
    i=i+1;
    %
end
%----------------------------------------------------------------
% Plotagem do deslocamento do nó central da treliça
%----------------------------------------------------------------
X0=[0 L0/2 L0];
Y0=[0 L0*sin(beta0) 0];
H=figure(2);
sgtitle('Arc-Length Method for 1 DoF Nonlinear Truss')
for j=2:i-1
    subplot(1,2,1)
    Y=[0 L0*sin(beta0)-u(j+1) 0];
    t1=plot(X0,Y0,'o-b',X0,Y,'o--k','LineWidth',1);
    line([L0/2 L0/2],[Y(2)+0.2 Y(2)+1.5],'linewidth',1.5,'Color','r');
    line([L0/2 L0/2],[Y(2)+0.2 Y(2)+0.2],'marker','V','MarkerSize',4,'linewidth',1,'Color','r')
    title('1 DoF Nonlinear Truss')
    legend([t1],'Original','Deformed','Location','northwest')
    axis([-0.05 1.05 -3 3]);
    axis square
    xlabel('x')
    ylabel('y')
    subplot(1,2,2)
    h1=plot(x,y,'b','LineWidth',1.5);
    hold on
    h2=plot(resposta(j,1),resposta(j,2),'o','Color','k');
    legend([h1 h2],'Analytical','Aproximated','Location','northwest')
    title('Path-following')
    xlabel('Displacement, u')
    ylabel('Load, \lambda*q')
    axis equal
    axis square
    pause(0.1)
    drawnow
    M=getframe(H);
    MM(j-1)=M;
end
video=VideoWriter('SimulCrisfield','MPEG-4');
video.FrameRate=6;
open(video);
writeVideo(video,MM);
close(video);
