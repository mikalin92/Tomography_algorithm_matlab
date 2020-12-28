% Mika L
%reconstructs picture from tomografical experiement
load('in_out_coordinate_data.mat')
load('radioactivity_count_data.mat')

%grid size
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
KuvM=zeros(msize(1),gsize(1)*gsize(2));



% 3 nested loops, inside nCr(4,2)=6 condition-checks 
%with 8 conditions
for k=1:msize(1)

A=[MitM(k,1) MitM(k,2)]';
B=[MitM(k,3) MitM(k,4)]';


M=zeros(gsize(1),gsize(2));

%loop every pixel
for m=1:gsize(1)
for n=1:gsize(2)
    
pa=A+((m-A(1))/(B(1)-A(1)))*(B-A);%starting point grid cell row 
pb=A+((m-1-A(1))/(B(1)-A(1)))*(B-A);%starting point grid cell col 
pc=A+((n-A(2))/(B(2)-A(2)))*(B-A);%ending point grid cell row 
pd=A+((n-1-A(2))/(B(2)-A(2)))*(B-A);%ending point grid cell col 


%check every nCr possibility, if radiation hits a grid, calculate 
%geometrical path length in that grid cell
%4 sides and in+out so nCr=6
if m-1<=pa(1) && pa(1)<=m && n-1<=pa(2) && pa(2)<=n && m-1<=pb(1) && pb(1)<=m && n-1<=pb(2) && pb(2)<=n
    if sqrt((pa(1)-pb(1))^2+(pa(2)-pb(2))^2)>0
    M(m,n)=sqrt((pa(1)-pb(1))^2+(pa(2)-pb(2))^2);
    end
end

if m-1<=pa(1) && pa(1)<=m && n-1<=pa(2) && pa(2)<=n && m-1<=pc(1) && pc(1)<=m && n-1<=pc(2) && pc(2)<=n
    if sqrt((pa(1)-pc(1))^2+(pa(2)-pc(2))^2)>0
    M(m,n)=sqrt((pa(1)-pc(1))^2+(pa(2)-pc(2))^2);
    end
end


if m-1<=pa(1) && pa(1)<=m && n-1<=pa(2) && pa(2)<=n && m-1<=pd(1) && pd(1)<=m && n-1<=pd(2) && pd(2)<=n
    if sqrt((pa(1)-pd(1))^2+(pa(2)-pd(2))^2)>0
    M(m,n)=sqrt((pa(1)-pd(1))^2+(pa(2)-pd(2))^2);
    end
end

if m-1<=pb(1) && pb(1)<=m && n-1<=pb(2) && pb(2)<=n && m-1<=pc(1) && pc(1)<=m && n-1<=pc(2) && pc(2)<=n
    if sqrt((pb(1)-pc(1))^2+(pb(2)-pc(2))^2)>0
    M(m,n)=sqrt((pb(1)-pc(1))^2+(pb(2)-pc(2))^2);
    end
end


if m-1<=pb(1) && pb(1)<=m && n-1<=pb(2) && pb(2)<=n && m-1<=pd(1) && pd(1)<=m && n-1<=pd(2) && pd(2)<=n
    
    if sqrt((pb(1)-pd(1))^2+(pb(2)-pd(2))^2) >0
    M(m,n)=sqrt((pb(1)-pd(1))^2+(pb(2)-pd(2))^2);   
    end
    
end

if m-1<=pc(1) && pc(1)<=m && n-1<=pc(2) && pc(2)<=n && m-1<=pd(1) && pd(1)<=m && n-1<=pd(2) && pd(2)<=n
    
    if sqrt((pc(1)-pd(1))^2+(pc(2)-pd(2))^2)>0
    M(m,n)=sqrt((pc(1)-pd(1))^2+(pc(2)-pd(2))^2);
    end
    
end


end
end

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
pic=reshape(picvec,gsize(1),gsize(2));
imagesc(pic);
title('Reconstructed tomographical image')
axis image


