function shootthemoon(v0,theta0,phi,tmax)
    c=[0;0.5;0.5;1];
    a=[0 0 0 0;0.5 0 0 0;0 0.5 0 0;0 0 1 0];
    w=[1/6 1/3 1/3 1/6];
    Me=6e24;
    Mm=7.35e22;
    G=6.67e-11*1e-9*3600^2;
    d=384000;
    re=Mm/(Me+Mm)*d;
    rm=Me/(Me+Mm)*d;
    T=27.322*24;
    omega=2*pi/T;
    phim=@(t) omega*t+phi;
    phie=@(t) omega*t+phi+pi;
    xe=@(t) re*cos(phie(t));
    ye=@(t) re*sin(phie(t));
    xm=@(t) rm*cos(phim(t));
    ym=@(t) rm*sin(phim(t));
    de=@(t,s) sqrt((s(1)-xe(t))^2+(s(3)-ye(t))^2);
    dm=@(t,s) sqrt((s(1)-xm(t))^2+(s(3)-ym(t))^2);
    Fx=@(t,s) -G*Me/(de(t,s)^3)*(s(1)-xe(t))-G*Mm/(dm(t,s)^3)*(s(1)-xm(t));
    Fy=@(t,s) -G*Me/(de(t,s)^3)*(s(3)-ye(t))-G*Mm/(dm(t,s)^3)*(s(3)-ym(t));
    F=@(t,s) [s(2);Fx(t,s); s(4);Fy(t,s)];
    s(:,1)=[6800;v0*cosd(theta0);0;v0*sind(theta0)];
    t(1)=0;
    dt=0.1;
    i=1;
    while t(i)<tmax
        k=zeros(length(s(:,1)),length(c));
            for j=1:length(c)
                k(:,j)=dt*F(t(i)+dt*c(j),s(:,i)+k*a(j,:)');
            end
            s(:,i+1)=s(:,i)+k*w';
            t(i+1)=t(i)+dt;
            i=i+1;
    end
    thetae=0:360;
    xep=xe(t);
    yep=ye(t);
    xmp=xm(t);
    ymp=ym(t);
    Ex=6400*cosd(thetae);
    Ey=6400*sind(thetae);
    Mx=3600*cosd(thetae);
    My=3600*sind(thetae);
    
    %plot(s(1,:),s(3,:));
    %hold on;
    %plot(xep,yep,'r');
    plot(Ex,Ey,'r');
    plot(Mm,My,'r');
    %plot(xmp,ymp,'g');
    for i=1:length(t)
        clf
        plot(s(1,i),s(3,i),'r:o');
        axis([-d d -d d]);
        hold on;
        plot(xep(i),yep(i),'bo','markersize',20,'markerfacecolor',[1 0 0]);
        plot(xmp(i),ymp(i),'bo','markersize',10,'markerfacecolor',[0 1 0]);
        
        plot(s(1,1:i),s(3,1:i),'r');
        plot(xep(1:i),yep(1:i),'b');
        plot(xmp(1:i),ymp(1:i),'b');
        pause(0.001);
        hold off;
        
end

   