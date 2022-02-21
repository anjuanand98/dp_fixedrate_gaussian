%Xnopt stores quantized regions
%Yn stores quantized points

function main
clc;
close all;
clear all;
%initializing x
num=300;%number of elements initialized in x range
Nval=[2 4 8 16 32];
%range of x: [a,b]
a=-3;
b=3;
%x parameters
mu=0;
sigma_xsq=1;
%samples of x
delta=(b-a)/num;
xval=[a+delta/2:delta:b-delta/2];
%pmf of samples of x
p=zeros(length(xval),1);
f1=@(xv) ((1/sqrt(2*pi*sigma_xsq))*exp(-(xv-mu).^2/(2*sigma_xsq)));
scal=1/integral(f1,a,b,'ArrayValued',true);
for i=1:length(p)
    p(i)=scal*integral(f1,xval(i)-delta/2,xval(i)+delta/2,'ArrayValued',true);
end
for N=Nval%number of quantization levels
[Xnopt]=nonstrategic_quantization_hloop(xval,N,p);%returns quantization decision levels
Yn=decoder(Xnopt,xval,p);%returns quantization representative levels
%encoder distortion
endistortion=encoderdistortion(Xnopt,Yn,xval,p);
%decoder distortion
dedistortion=decoderdistortion(Xnopt,Yn,xval,p);
%storing data for N=2 4 8 16 32
switch N
    case 2
        Xnopt2=Xnopt;
        Yn2=Yn;
        endistortion2=endistortion;
        dedistortion2=dedistortion;
    case 4
        Xnopt4=Xnopt;
        Yn4=Yn;
        endistortion4=endistortion;
        dedistortion4=dedistortion;
    case 8
        Xnopt8=Xnopt;
        Yn8=Yn;
        endistortion8=endistortion;
        dedistortion8=dedistortion;
    case 16
        Xnopt16=Xnopt;
        Yn16=Yn;
        endistortion16=endistortion;
        dedistortion16=dedistortion;
    case 32
        Xnopt32=Xnopt;
        Yn32=Yn;
        endistortion32=endistortion;
        dedistortion32=dedistortion;
end
end
% % %dataxcubedgaussianfixed.mat holds all the output data
save('dataxcubedgaussianfixed.mat','Nval','Xnopt2','Yn2','endistortion2','dedistortion2','Xnopt4','Yn4','endistortion4','dedistortion4','Xnopt8','Yn8','endistortion8','dedistortion8','Xnopt16','Yn16','endistortion16','dedistortion16','Xnopt32','Yn32','endistortion32','dedistortion32');
%plotting encoder and decoder distortions
f=figure;
plot(Nval,[endistortion2, endistortion4, endistortion8, endistortion16, endistortion32],'*-');
hold on;
plot(Nval,[dedistortion2, dedistortion4, dedistortion8, dedistortion16, dedistortion32],'o-');
hold off;
legend({'encoder distortion','decoder distortion'},'FontSize',14);
xlabel('N','FontSize',14);
ylabel('distortion','FontSize',14);
title('fixed rate xcubed gaussian','FontSize',14);
saveas(f,'fixedrate_xcubed_gaussian.fig');
saveas(f,'fixedrate_xcubed_uniform.png');
saveas(f,'fixedrate_xcubed_gaussian.m');

function [y]=decoder(x,xval,p)
%inputs: x - quantization decision levels, xval - samples of x, p - pmf of xval
%outputs: y - quantization representative levels
N=length(x)-1;
y=zeros(1,N);
for i=1:N%iterate over each region
    in1=find(x(i)==xval);
    if in1~=1%non-overlapping regions (]
        in1=in1+1;
    end
    in2=find(x(i+1)==xval);
    y(i)=(xval(in1:in2)*p(in1:in2))/(ones(1,length(p(in1:in2)))*p(in1:in2));
end

function endistortion=encoderdistortion(x,y,xval,p)
%inputs: x - quantization decision levels, y - quantization representative
%levels, xval - samples of x, p - pmf of xval
%output: endistortion - encoder distortion
N=length(y);
endistortion=0;
    for n=1:N%iterate over each region
        in1=find(x(n)==xval);
        if in1~=1%non-overlapping regions (]
            in1=in1+1;
        end
        in2=find(x(n+1)==xval);
        endistortion=endistortion+((xval(in1:in2).^3-y(n)).^2)*p(in1:in2);  
    end
    
function dedistortion=decoderdistortion(x,y,xval,p)
%inputs: x - quantization decision levels, y - quantization representative
%levels, xval - samples of x, p - pmf of xval
%output: dedistortion - decoder distortion
N=length(y);
dedistortion=0;
for n=1:N%iterate over each region
    in1=find(x(n)==xval);
    if in1~=1%non-overlapping regions (]
        in1=in1+1;
    end
    in2=find(x(n+1)==xval);
    dedistortion=dedistortion+((xval(in1:in2)-y(n)).^2)*p(in1:in2);  
end 

function [Xnopt]=nonstrategic_quantization_hloop(xval,N,p)
%inputs: xval - samples of x, N - number of quantization levels, p - pmf of xval
%output: Xnopt - quantization decision levels
%initializing x optimal end points
Xnopt=10^7*ones(N+1,1);
Xnopt(1)=xval(1);
Xnopt(end)=xval(end);
%initializing [alpha, beta] pairs for x
D1val=zeros((length(xval)-1)*length(xval)/2,2);%stores all possible intervals [alpha,beta]
k=1;%iteration over all possible (alpha,beta) pairs
for i=1:length(xval)-1%for each alpha value possible
    for j=i+1:length(xval)%for each beta value possible, given the alpha value
        D1val(k,1:2)=[xval(i) xval(j)];%(alpha,beta)
        k=k+1;
    end
end
Dnxoxval=D1val(1:length(xval)-1,2);%(xo,x) - possible x values
%initializing D1 for each [alpha, beta] pair in D1val
D1=zeros(size(D1val,1),1);%D1([alpha,beta])=minimum over y integral(alpha to beta (f(x,x-y)p(x)dx))
%calculating D1
for i=1:size(D1val,1)%iterate over possible (alpha,beta) pairs
    q=decoder([D1val(i,1);D1val(i,2)],xval,p);
    D1(i)=encoderdistortion([D1val(i,1);D1val(i,2)],q,xval,p);%distortion D1
end
Dn=10^7*ones(length(Dnxoxval),N-1);%stores Dn values for each n (2:N), for each x in (xo,x)
Xn=10^7*ones(length(Dnxoxval),N-1);%stores X_(n-1) values for each n (1:N-1), for each x in (xo,x)
for n=2:N
    for i=n+1:length(xval)%Dn(xo,x(i)) - there should be sufficient levels between xo and x(i), so i starts from n, e.g., n=2 level quantization implies you need atleast three points in x 
        [minval xboundary]=nonstr_recursion_1(xval,i,n,Xn,Dn,Dnxoxval,D1val,D1);%finding Dn and X_(n-1) values for each (xo,x) for a given n level 
        %indexed such that x(i) location in Dnxoxval gives corresponding Dn and Xn values
        Dn(find((Dnxoxval==xval(i))),n-1)=minval;
        Xn(find((Dnxoxval==xval(i))),n-1)=xboundary;
    end
end
for n=N:-1:2%backward iteration to find Xnopt
    Xnopt(n)=Xn(find((Dnxoxval==Xnopt(n+1))),n-1);%X_(n-1) optimal=X_(n-1)(xo,X_n opt)
end

function [minval xboundary]=nonstr_recursion_1(xval,xind,n,Xn,Dn,Dnxoxval,D1val,D1)%returns distortion and corresponding value of x and theta
    if n==1%n=1 values are already computed in D1
        minval=D1(find(D1val(:,1)==xval(1) & D1val(:,2)==xval(xind)));
        xboundary=xval(xind);
        return;
    end
    alpharange=[n:xind-1];%indices of alpha values possible in x range (xo<alpha<x)
    arr=zeros(length(alpharange),1);%for each alpha value possible, finding minimum (D_(n)(xo,x)=min over alpha, ao<alpha<x D_(n-1)(xo,alpha)+D1(alpha,x))
    %D_n(xo<alpha<x, thetao<beta<theta)
    for arrx=1:length(alpharange)
        if n-2>=1 && Xn(find((Dnxoxval==xval(alpharange(arrx)))),n-2)~=10^7
            xboundary1=Xn(find((Dnxoxval==xval(alpharange(arrx)))),n-2);
            minval1=Dn(find((Dnxoxval==xval(alpharange(arrx)))),n-2);
        else
        [minval1 xboundary1]=nonstr_recursion_1(xval,alpharange(arrx),n-1,Xn,Dn,Dnxoxval,D1val,D1);%D_(n-1)(xo,alpha) 
        end
        %D_n(xo,x)=min over alpha, xo<alpha<x (D_(n-1)(xo,alpha)+D1(alpha,x))
        arr(arrx)=minval1+D1(find(D1val(:,1)==xval(alpharange(arrx)) & D1val(:,2)==xval(xind)));
    end
    minval=min(arr);%returning minimum distortion value
    in=find(arr==minval);
    xboundary=xval(alpharange(in(1)));