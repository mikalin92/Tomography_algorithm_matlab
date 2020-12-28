% Mika L
%reconstructs picture from tomografical experiement
%notes added for change to 3 dimensions (not tested)

load('in_out_coordinate_data.mat')
load('radioactivity_count_data.mat')

%grid size
%if 3D is wanted this should have 3 numbers
gsize=[8 8]';

%put coordinate data into a matrix
MitM=[
Mit1;%0 directly
Mit2;%on straight 90 angle  
Mit3;
Mit4;%on angle 45
Mit5;
Mit6;%ending,-45
Mit7;
Mit8;%"hand fan"1
Mit9;
Mit10;%"hand fan"2
Mit11;
Mit12;%"hand fan"3


];


rdivec=[
rdi1;
rdi2;
rdi3;
rdi4;
rdi5;
rdi6;
rdi7;
rdi8;
rdi9;
rdi10;
rdi11;
rdi12;


];

%solution is based to taking the logarithm:
logvec=-log(rdivec/max(rdivec));

msize=size(MitM);

%init pic matrix
%3D should have *gsize(3) added
KuvM=zeros(msize(1),gsize(1)*gsize(2));



% 4 nested loops, inside nCr(4,2)=6 condition-checks 
%with 8 conditions
%3d-> 15 checks to be

for k=1:msize(1)

A=[MitM(k,1) MitM(k,2)]';
B=[MitM(k,3) MitM(k,4)]';


M=zeros(gsize(1),gsize(2));

%loop every pixel
%if 3D, for loop 1:gsize(3) to be added
for m=1:gsize(1)
for n=1:gsize(2)

    
p=zeros(2,4);%if 3D, zeros(3,6)

%if 3d add more lines here:
p(:,1)=A+((m-A(1))/(B(1)-A(1)))*(B-A);%starting point grid cell row 
p(:,2)=A+((m-1-A(1))/(B(1)-A(1)))*(B-A);%starting point grid cell col 
p(:,3)=A+((n-A(2))/(B(2)-A(2)))*(B-A);%ending point grid cell row 
p(:,4)=A+((n-1-A(2))/(B(2)-A(2)))*(B-A);%ending point grid cell col 


%check every nCr possibility, if radiation hits a grid, calculate 
%geometrical path length in that grid cell
%4 sides and in+out so nCr=6
%if 3D is wanted this, with other changes, should be 6 sides and in+out
%nCr=15
%run thru 1: 6
combM=combnk(1:4,2);
combsz=size(combM);
for idex=1:combsz(1)

    %if 3d, more conditions to be nested
if m-1<=p(1,combM(idex,1)) && p(1)<=m && n-1<=p(2,combM(idex,1)) && p(2,combM(idex,1))<=n && m-1<=p(1,combM(idex,2)) && p(1,combM(idex,2))<=m && n-1<=p(2,combM(idex,2)) && p(2,combM(idex,2))<=n
    if sqrt((p(1,combM(idex,1))-p(1,combM(idex,2)))^2+(p(2,combM(idex,1))-p(2,combM(idex,2)))^2)>0
    M(m,n)=sqrt((p(1,combM(idex,1))-p(1,combM(idex,2)))^2+(p(2,combM(idex,1))-p(2,combM(idex,2)))^2);
    end
end

end


end
end

%for 3d, add *gsize(3)
Mres=reshape(M,[gsize(1)*gsize(2),1]);
KuvM(k,:)=Mres;

end



idim=size(KuvM'*KuvM);
I=eye(idim(1));


%%%%Andrei Tikhonov's regularization
L0=I;
L1=diff(I);
L2=diff(diff(I));

tikho0=L0'*L0;
tikho1=L1'*L1;
tikho2=L2'*L2;


alpha=10^(-1);
tikho=tikho2;


picvec=(KuvM'*KuvM+alpha*tikho)\KuvM'*logvec;
%matlab operand x\y means inverse(x)*y

%if 3d, use another matlab tool to show it
pic=reshape(picvec,gsize(1),gsize(2));
imagesc(pic);
title('Reconstructed tomographical image')
axis image


