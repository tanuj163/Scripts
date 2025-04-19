
clear all
close all
clc

%bond length(Angstrom). 
l=2.42;

%Basis vector
a= l*cosd(81.01/2);
b=l*sind(81.02/2);

%Layer thickness
thickness=2*b;

%Vaccum length in z-direction
z=10;
%Bravis lattice 
a1=[3*a 0   0];
a2=[0 sqrt(3)*a  0];
a3=[0 0 thickness];
%position of basis atoms
Basis=[
    1/2*a2;
    1/2*a2+a3;
    1/6*a1+a3/2;
    1/2*a1;
    1/2*a1+a3;
    2/3*a1+1/2*a2+a3/2
    ];
BasisType=[
    1
    1
    2
    1
    1
    2
    ];
%number of atoms in the basis
NB=length(Basis(:,1));
%--------------------
% bulid the lattice
%--------------------
% build the lattice structure

s1 = input('Number of unit cell in a_1 direction:  ','s');
s2 = input('Number of unit cell in a_2 direction:  ','s');
N1=str2num(s1);
N2=str2num(s2);

NumAtom=N1*N2*NB;
coord=zeros(NumAtom, 3);
ticker=1;
for i=0:(N1-1)
    for j=0:(N2-1)
        for k=1:NB
            coord(ticker+k-1,1)=ticker+k-1;
            coord(ticker+k-1,2)=BasisType(k);
            coord(ticker+k-1,3:5)=i*a1+j*a2+Basis(k,:);
        end
        ticker=ticker+NB;
    end
end
        
plot(coord(:,3),coord(:,4),'o');

%Shifting all the atoms in z-direction, for the vaccum space
coord(:,5)=coord(:,5)+ z;
axis equal;
xlabel('x (Angstrom)');
ylabel('y (Angstrom)');

%----------------------
% print to a data file
%----------------------
fid=fopen('WS2.xyz','w');
fprintf(fid,'%i \n',NumAtom);
fprintf(fid,'Tungsten Sulphide lattice %i x %i \n', N1, N2);
for i=1:NumAtom
    fprintf(fid,'%i %12.8f %12.8f %12.8f\n',coord(i,2:5));
end
fclose(fid);
%--------------------------------------
% print to a lammp structure data file
%--------------------------------------
fid=fopen('WS2.u','w');
fprintf(fid,'#Position data for Tungsten Sulphide %i x %i \n', N1, N2);

fprintf(fid,'%i atoms\n',NumAtom);
fprintf(fid,'%i atom types\n\n',2);

fprintf(fid,'%12.8f %12.8f  xlo xhi\n',0, N1*norm(a1));
fprintf(fid,'%12.8f %12.8f  ylo yhi\n',0, N2*norm(a2));
fprintf(fid,'%12.8f %12.8f  zlo zhi\n\n',0, 2*z+thickness);

fprintf(fid,'Masses\n\n');
fprintf(fid,'1 32.0600  \n');
fprintf(fid,'2 183.8400 \n\n');

fprintf(fid,'Atoms\n\n');

for i=1:NumAtom
    fprintf(fid,'%i %i %12.8f %12.8f %12.8f\n',coord(i,:));
end
fclose(fid);


        