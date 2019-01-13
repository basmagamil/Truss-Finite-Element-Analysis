function [ratiof,stressf,ratiob,stressb] = analyze_structure(h,Ab)

%take input as [element no, node no1, node no2, theta with positive x]
enn=[1 1 2 0;
    2 2 3 0;
    3 3 4 0;
    4 4 5 0;
    5 5 6 0;
    6 7 1 120;
    7 7 2 60;
    8 8 2 120;
    9 8 3 60;
    10 9 3 120;
    11 9 4 60;
    12 10 4 120;
    13 10 5 60;
    14 11 5 120;
    15 11 6 60;
    16 7 8 0;
    17 8 9 0;
    18 9 10 0;
    19 10 11 0];
ennframe=enn(1:5,:);
ennbar=enn(6:19,:);

%indexing nodes
nodesframe=1:6;
nodesbar=7:11;

%indexing dof
dofframe=1:1:length(nodesframe)*3;
dofframe=reshape(dofframe,3,[])';
dofbar=1+length(nodesframe)*3:1:length(nodesframe)*3+length(nodesbar)*2;
dofbar=reshape(dofbar,2,[])';

%ndoftotal=[node no, u, v, theta]
n3dofframe=[nodesframe' dofframe];
n2dofbar=[nodesbar' dofbar];
n2dofallbar=[n2dofbar zeros(length(n2dofbar),1)];
ndoftotal=[n3dofframe; n2dofallbar];

%connectframe=[element number, node1 dof, node2 dof];
connectframe=[];
for m=1:length(ennframe)
    connectframe=[connectframe; m ndoftotal(enn(m,2),2:4) ndoftotal(enn(m,3),2:4)];
end
connectbar=[];
for m=ennbar(1,1):ennbar(length(ennbar),1)
    connectbar=[connectbar; m ndoftotal(enn(m,2),2:3) ndoftotal(enn(m,3),2:3)];
end

%givens
L=2;
b=3.5;
Ef=41e09;
Eb=207e09;
wd=900*9.8;
%calculate
I=(1/12)*b*(h^3); 
Af=b*h;

%Building k local & global frame & truss for all elements
for i=1:5
    klf=[(Ef*Af)/L 0 0 -(Ef*Af)/L 0 0;
        0 (12*Ef*I)/L^3 (6*Ef*I)/L^2 0 -(12*Ef*I)/L^3 (6*Ef*I)/L^2;
        0 (6*Ef*I)/L^2 (4*Ef*I)/L 0 -(6*Ef*I)/L^2 (2*Ef*I)/L;
        -(Ef*Af)/L 0 0 (Ef*Af)/L 0 0;
        0 -(12*Ef*I)/L^3 -(6*Ef*I)/L^2 0 (12*Ef*I)/L^3 -(6*Ef*I)/L^2;
        0 (6*Ef*I)/L^2 (2*Ef*I)/L 0 -(6*Ef*I)/L^2 (4*Ef*I)/L];
    kgf=klf;
    kgftotal(:,:,i)=kgf;
end

for i=6:19
    klb=((Eb*Ab)/L).*[1 0 -1 0;
        0 0 0 0;
        -1 0 1 0;
        0 0 0 0];
    R=[cosd(enn(i,4)) sind(enn(i,4)); -sind(enn(i,4)) cosd(enn(i,4))];
    Tb=[R zeros(2,2); zeros(2,2) R];
    kgb=transpose(Tb)*klb*Tb;
    kgbtotal(:,:,i)=kgb;
end

%%%preparing for kglobal
%frame
rowf=[];
colf=[];
kgfelements=[];
for m=1:5
    for i=2:7
        for j=2:7
            rowf=[rowf; connectframe(m,i)];
            colf=[colf; connectframe(m,j)];
            kgfelements= [kgfelements; kgftotal(i-1,j-1,m)];
        end
    end
end
%bar
rowb=[];
colb=[];
kgbelements=[];
for m=6:19
    for i=2:5
        for j=2:5
            rowb=[rowb; connectbar(m-5,i)];
            colb=[colb; connectbar(m-5,j)];
            kgbelements= [kgbelements; kgbtotal(i-1,j-1,m)];
        end
    end
end
%sum & sparse
row=[rowf; rowb];
col=[colf; colb];
kgelements=[kgfelements; kgbelements];
kgtotal=sparse(row,col,kgelements);
kg=full(kgtotal);
kg28=kg;

%%LOAD
%local global load
for i=1:5
    Lgf=[-wd*L/2; -(wd*L^2)/12; -wd*L/2; (wd*L^2)/12];
    Lftotal(:,:,i)=Lgf;
end

connectLf=[connectframe(:,1), connectframe(:,3), connectframe(:,4), connectframe(:,6), connectframe(:,7)];

%load sparse
%frame
rowLf=[];
Lfelements=[];
for m=1:5
    for i=2:5
        rowLf=[rowLf; connectLf(m,i)];
        Lfelements= [Lfelements; Lftotal(i-1,1,m)];
    end
end
%set bars to zero
rowLb=[];
Lbelements=[];
for m=6:19
    for i=2:5
        rowLb=[rowLb; connectbar(m-5,i)];
        Lbelements= [Lbelements; 0];
    end
end
%sum & sparse
rowL=[rowLf; rowLb];
Lgelements=[Lfelements; Lbelements];
Lgtotal=sparse(rowL,ones(length(rowL),1),Lgelements);
Lg=full(Lgtotal);
Lg28=Lg;

%%%Deleting rows of delta corresponding supports
kg(1,:)=[]; kg(:,1)=[];
kg((2-1),:)=[]; kg(:,(2-1))=[];
kg((17-2),:)=[]; kg(:,(17-2))=[];

Lg(1,:)=[]; 
Lg((2-1),:)=[]; 
Lg((17-2),:)=[]; 

%calc deltas
delta25=kg\Lg;
delta28=[0;0;delta25(1:17-2-1,:);0;delta25(17-2:25)];

%%%internal forces
for i=1:5
    ff(:,:,i)=kgftotal(:,:,i)*[delta28(connectframe(i,2:7),1)];
    M(i,1)=max(abs(ff(3,:,i)),abs(ff(6,:,i)));
    stressf(i,1)=(abs(ff(1,:,i))/Af)+(M(i,1)*(h/2)/I);
end
for i=6:19
    fb(:,:,i)=kgbtotal(:,:,i)*[delta28(connectbar(i-5,2:5),1)];
    R=[cosd(enn(i,4)) sind(enn(i,4)); -sind(enn(i,4)) cosd(enn(i,4))];
    Tb=[R zeros(2,2); zeros(2,2) R];
    flb(:,:,i)=Tb*fb(:,:,i);
    stressb(i-5,1)=abs(flb(3,:,i))/Ab;
end
stressfallow=130e06;
stressballow=220e06;
ratiof=stressf./stressfallow;
ratiob=stressb./stressballow;

end

