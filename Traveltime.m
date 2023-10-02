clear all;
close all;
clc;
tic;

L=7000; % length
D=2500; % depth
d=2000.0; % depth of soil layer from the bottom
Hr=50; % regional depth
HL=50; % local depth
P=365; % period
Kx1=2; % conductivity of bottom layer along x
Kz1=2; % conductivity of bottom layer along z
Kx2=2; % conductivity of top layer along x
Kz2=2; % conductivity of top layer along z
Ss1=.0001; % specific storage of bottom layer along
Ss2=.0001; % specific storage of top layer along
Dx1=Kx1/Ss1;
Dx2=Kx2/Ss2;
bta1=(Kz1/Kx1)^.5;
bta2=(Kz2/Kx2)^.5;
An1=(Kz1/Kx1)^.5;
An2=(Kz2/Kx2)^.5;
M=7; % No of local undulations

ne1=.25; % porosity of layer 1
ne2=.25; % porosity of layer 2

mloop=19;
mloop1=19;
nloop=99999;
nloop1=19;
tloop=210000;
xa=[6900];
dz=.0001;

%steady-state solution
for m=1:mloop

    Nm=(m-1)*pi/L;

A1=[1 -1 0 0];
aa=Nm*d/An1;
ab=Nm*d/An2;
A2=[exp(aa) exp(-aa)  -exp(ab) -exp(-ab)];
ac=Kz1/An1;
ad=Kz2/An2;
A3=[ac*exp(aa) -ac*exp(-aa)  -ad*exp(ab) ad*exp(-ab)];
ae=Nm*D/An2;
A4=[0 0 exp(ae) exp(-ae)];
A=[A1 
    A2
    A3
    A4];
B=[0 0 0 1]';
C=A\B;
C1(m)=C(1);
C2(m)=C(2);
C3(m)=C(3);
C4(m)=C(4);

if m==1
    Ams(m)=D+Hr;
elseif m==2
    Ams(m)=-Hr;
else
    Ams(m)=0;
end
end

for m=1:mloop
    %%for sin
if m==1
    Am0(m)=D+Hr;
elseif  m==2
    Am0(m)=-Hr;
elseif m==M+1
    Am0(m)=0;
else
    Am0(m)=0;
end

end

% transient solution

for m=1:mloop1
    
    Nm=(m-1)*pi/L;

    s=0;
    for n=1:nloop
        
        x1=(Nm*Nm*Dx1)^.5;
        x2=(Nm*Nm*Dx2)^.5;
        if x1>x2
            x3=x2;
        else
            x3=x1;
        end
        y1=x3+((n-1)*dz);
        y2=x3+(n*dz);
        
        Kz=Kz1/Kz2;
        bta=bta1/bta2;
        
        aa=y1^2/Dx1;
        r11=(aa-Nm^2)^.5;
        ba=y2^2/Dx1;
        r12=(ba-Nm^2)^.5;
        
        bb=y1^2/Dx2;
        r21=(bb-Nm^2)^.5;
        bc=y2^2/Dx2;
        r22=(bc-Nm^2)^.5;
        
        ab=sin(r11*d/bta1);
        ac=sin(r21*(D-d)/bta2);
        g1=Kz*r11*ab*ac;
        ad=cos(r11*d/bta1);
        ae=cos(r21*(D-d)/bta2);
        g2=bta*r21*ad*ae;
        z1=g1-g2;
        
        af=sin(r12*d/bta1);
        ag=sin(r22*(D-d)/bta2);
        g3=Kz*r12*af*ag;
        ah=cos(r12*d/bta1);
        ai=cos(r22*(D-d)/bta2);
        g4=bta*r22*ah*ai;
        z2=g3-g4;
        
        if z1>=0 && z2>=0
            aa=0;
           
        elseif z1<=0 && z2<=0
            aa=0;

        elseif z1>0 && z2<0
            
            aa=abs(z1)+abs(z2);
            ab=dz*abs(z1)/aa;
            ys=y1+ab;
            ab=1;
            s=s+ab;
            Nn(m,s)=ys;
            if s>nloop1
                break
            end
            
        else 
            aa=abs(z1)+abs(z2);
            ab=dz*abs(z1)/aa;
            ys=y1+ab;            
            ab=1;
            s=s+ab;
            Nn(m,s)=ys;
            if s>nloop1
                break
            end
        end
    end
end

for m=1:mloop1
    
    Nm=(m-1)*pi/L;
    for n=1:nloop1+1
    
    a1=(Nn(m,n)^2)/Dx1;
    a2=(Nn(m,n)^2)/Dx2;
    r1=(a1-Nm^2)^.5;
    r2=(a2-Nm^2)^.5;

    ca=Nm/An1;
    cb=Nm/An2;

    f1=@(y)(exp(ca*y).*cos(r1*y/bta1));
    w1=integral(f1,0,d);
    f2=@(y)(exp(-ca*y).*cos(r1*y/bta1));
    w2=integral(f2,0,d);
    ca=Ss1*((w1*C1(m))+(w2*C2(m)));

    f3=@(y)(exp(cb*y).*sin(r2*(D-y)/bta2));
    w3=integral(f3,d,D);
    f4=@(y)(exp(-cb*y).*sin(r2*(D-y)/bta2));
    w4=integral(f4,d,D);

    cb=(w3*C3(m))+(w4*C4(m));
    cc=cos(r1*d/bta1)/sin(r2*(D-d)/bta2);
    cd=Ss2*cb*cc;

    f5=@(y)(cos(r1*y/bta1).*cos(r1*y/bta1));
    w5=integral(f5,0,d);
    
    f6=@(y)(sin(r2*(D-y)/bta2).*sin(r2*(D-y)/bta2));
    w6=integral(f6,d,D);

    ce=Ss1*w5;
    cf=Ss2*w6*cc^2;
    cg=Ams(m)-Am0(m);

    Emno(m,n)=cg*(ca+cd)/(ce+cf);
    Gmn(m,n)=(ca+cd)/(ce+cf);
    end
end


% travel time and pathline calculation

for v=1:length(xa)
z=D;
x=xa(v);
X(1)=x;
Z(1)=z;
dt=.4;
for u=1:tloop
    t=(u-1)*dt;

    alp=sin(pi*t/P)^2;

    for m=1:mloop
if m==1
    Amt(m)=D+Hr+(alp*HL);
elseif  m==2
    Amt(m)=-Hr;
elseif m==M+1
    Amt(m)=-alp*HL;
else
    Amt(m)=0;
end

end

    for m=1:mloop1
    
    Nm=(m-1)*pi/L;
    for n=1:nloop1+1
    
        dd=Nn(m,n);
        ea=pi/P;
        eb=dd^2*sin(2*ea*t);
        ec=2*ea*cos(2*ea*t);
        ed=2*ea*exp(-t*dd^2);
        ee=eb-ec+ed;
        ef=(2*ea)^2+(dd^4);

        if m==1
            fa=HL*ea*ee/ef;
        elseif m==M+1
            fa=-HL*ea*ee/ef;
        else
            fa=0;
        end

        eg=Emno(m,n)*exp(-t*dd^2);
        eh=fa*Gmn(m,n);

        Emnt(m,n)=eg-eh;
    end
end

    if z>=d

sb=0;
for m=1:mloop1
    
    Nm=(m-1)*pi/L;
    for n=1:nloop1+1
    
    a1=(Nn(m,n)^2)/Dx1;
    a2=(Nn(m,n)^2)/Dx2;
    r1=(a1-Nm^2)^.5;
    r2=(a2-Nm^2)^.5;

    bc=cos(r1*d/bta1)/sin(r2*(D-d)/bta2);
    bd=-bc*r2*cos(r2*(D-z)/bta2);
    be=Emnt(m,n)*bd.*cos(Nm*x)/bta2;
    sb=sb+be;
    end
end

s2=0;
for m=1:mloop

    Nm=(m-1)*pi/L;
    Z2=C3(m)*exp(Nm*z/An2)-C4(m)*exp(-Nm*z/An2);
    bb=Amt(m)*cos(Nm*x).*Z2*Nm/An2;
    s2=s2+bb;
end

Vz2=-Kz2*(s2+sb);
dz=Vz2*dt/ne2;
z=z+dz;
Z(u+1)=z;

sb=0;
for m=1:mloop1
    
    Nm=(m-1)*pi/L;
    for n=1:nloop1+1
    
    a1=(Nn(m,n)^2)/Dx1;
    a2=(Nn(m,n)^2)/Dx2;
    r1=(a1-Nm^2)^.5;
    r2=(a2-Nm^2)^.5;

    bc=cos(r1*d/bta1)/sin(r2*(D-d)/bta2);
    bd=-Nm*bc*sin(r2*(D-z)/bta2);
    be=Emnt(m,n)*bd.*sin(Nm*x);
    sb=sb+be;
    end
end

s2=0;
for m=1:mloop

    Nm=(m-1)*pi/L;
    Z2=C3(m)*exp(Nm*z/An2)+C4(m)*exp(-Nm*z/An2);
    bb=-Amt(m)*Nm*sin(Nm*x).*Z2;
    s2=s2+bb;
end

Vx2=-Kx2*(s2+sb);
dx=Vx2*dt/ne2;
x=x+dx;
X(u+1)=x;

    else    

sa=0;
for m=1:mloop1
    
    Nm=(m-1)*pi/L;
    for n=1:nloop1+1
    
    a1=(Nn(m,n)^2)/Dx1;
    r1=(a1-Nm^2)^.5;

    ba=-r1*sin(r1*z/bta1);
    bb=Emnt(m,n)*ba.*cos(Nm*x)/bta1;
    sa=sa+bb;
    end
end

s1=0;
for m=1:mloop

    Nm=(m-1)*pi/L;
    Z1=C1(m)*exp(Nm*z/An1)-C2(m)*exp(-Nm*z/An1);

    ba=Amt(m)*cos(Nm*x).*Z1*Nm/An1;
    s1=s1+ba;
end

Vz1=-Kz1*(s1+sa);
dz=Vz1*dt/ne1;
z=z+dz;
Z(u+1)=z;

sa=0;
for m=1:mloop1
    
    Nm=(m-1)*pi/L;
    for n=1:nloop1+1
    
    a1=(Nn(m,n)^2)/Dx1;
    r1=(a1-Nm^2)^.5;

    ba=-Nm*cos(r1*z/bta1);
    bb=Emnt(m,n)*ba.*sin(Nm*x);
    sa=sa+bb;
    end
end

s1=0;
for m=1:mloop

    Nm=(m-1)*pi/L;
    Z1=C1(m)*exp(Nm*z/An1)+C2(m)*exp(-Nm*z/An1);

    ba=-Amt(m)*Nm*sin(Nm*x).*Z1;
    s1=s1+ba;
end

Vx1=-Kx1*(s1+sa);
dx=Vx1*dt/ne1;
x=x+dx;
X(u+1)=x;
    end

    if z>D
        break
    end
end

Tt(v)=t/365; % travel time
hold on
plot(X,Z)
end
toc;
disp(Tt)

 