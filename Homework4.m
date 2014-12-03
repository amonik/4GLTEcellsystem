clear all
close all
clear classes

BW=2e9;%Hz
N=2048;%Resource Elements
Wk=15000;%Hz
PTx=100;%Watts
distance=2000;%diameter meters
alpha=2.2;
P_d=(10/1000)^2.5;%Path Loss at 10 metes
Hk=1;%AWGN Channel
No=1e-8;%Watts/Hz
Nvar=No*Wk;
Rb=10e6;%Bis/Sec

CellNet=CellNetwork(BW,N,Wk,PTx,distance,alpha,P_d,Hk,No,Nvar,Rb);

%hold off
title('Eb/No with no Interference - Cells with 1Km Radius')
mesh(CellNet.target_matrix_x,CellNet.target_matrix_y,CellNet.G_snr)
colorbar
hold off

figure()
title('Eb/No with no Interference - Cells with 1Km Radius')
mesh(CellNet.target_matrix_x,CellNet.target_matrix_y,CellNet.G_snr)
colorbar
xlabel('X Label meters')
ylabel('Y Label meters')
zlabel('Eb/No with no Interference')

%History of users
figure()
hist(reshape(CellNet.users, prod(size(CellNet.users)), 1),101)
title('Histogram of user location within cell')
xlabel('X Label Meters')
ylabel('Number of Users')

G2=round((CellNet.G_snr*100))/100;
figure()
hist(reshape(G2, prod(size(G2)), 1),101)
title('Histogram of G-Factor - Eb/No with no Interference')
xlabel('G-Factor - Eb/No with no Interference')
ylabel('Number of Observations')

[vect1 vect2]=hist(reshape(G2, prod(size(G2)), 1),unique(G2));

for p=1:length(unique(G2))
    probability(p)=sum(vect1(1:p))/sum(vect1);
end

figure()
plot(vect2,probability)
title('Probability Density Function - P(X_GdB <= GdB)')
xlabel('G-Factor - Eb/No with no Interference')
ylabel('Probability')

%Section 2 - Eb/No with Interferernce
CellNet=CellNetwork(BW,N,Wk,PTx,distance,alpha,P_d,Hk,No,Nvar,Rb);

%hold off
title('Eb/No with Interference - Cells with 1Km Radius')
mesh(CellNet.target_matrix_x,CellNet.target_matrix_y,CellNet.G_snir)
colorbar
hold off

figure()
title('Eb/No with Interference - Cells with 1Km Radius')
mesh(CellNet.target_matrix_x,CellNet.target_matrix_y,CellNet.G_snir)
colorbar
xlabel('X Label meters')
ylabel('Y Label meters')
zlabel('Eb/No with Interference')

%History of users
figure()
hist(reshape(CellNet.users, prod(size(CellNet.users)), 1),101)
title('Histogram of user location within cell')
xlabel('X Label Meters')
ylabel('Number of Users')

G2=round((CellNet.G_snir*100))/100;
figure()
hist(reshape(G2, prod(size(G2)), 1),101)
title('Histogram of G-Factor - Eb/No with Interference')
xlabel('G-Factor - Eb/No with Interference')
ylabel('Number of Observations')

[vect1 vect2]=hist(reshape(G2, prod(size(G2)), 1),unique(G2));

for p=1:length(unique(G2))
    probability(p)=sum(vect1(1:p))/sum(vect1);
end

figure()
plot(vect2,probability)
title('Probability Density Function - P(X_GdB <= GdB)')
xlabel('G-Factor - Eb/No with Interference')
ylabel('Probability')