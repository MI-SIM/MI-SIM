function dP=four_mod(time, init, S1in, D, Y1, kdec1, Y2, kdec2, Y3, kdec3, km1, Ks1, km2, Ks2, km3, Ks3, Ks3c, KI2,gamma0,gamma1,gamma2,S2in,S3in,h,h1,tt,motif,flag,cno)
%four_mod - calls motif function to be solved with the ODE solver
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: SOFTWARE HOSTED SITE
% September 2015; Last revision: 04-Jan-2016


%Update solver progress according to routine run

if isempty(flag)
    handles.timestamp=h;
    handles.progress=h1;
    c = clock; % [year month day hour minute seconds]
    TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
    total_tim=(time+1)*100/(tt+1);
    
    set(handles.timestamp,'String',['Time: ',TSMP])
    set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
    
elseif flag==3
    handles.timestamp=h;
    handles.progress=h1;
    c = clock; % [year month day hour minute seconds]
    TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
    tm=cno+time;
    total_tim=(tm+1)*100/(tt+1);
    
    set(handles.timestamp,'String',['Time: ',TSMP])
    set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
    
else
    handles.timestamp=h;
    c = clock; % [year month day hour minute seconds]
    TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
    set(handles.timestamp,'String',['Time: ',TSMP])
end
%allow for object to update
pause(0.000000000000000001)

%Take the substrate and biomass values at each time point
outsiz=length(init);
S1=init(1);
X1=init(2);
S2=S1; 
I2=1; %Default inhibition term

%Set initial conditions dependent on motif chosen
if outsiz==3
    X2=init(3);
elseif outsiz==4
    S2=init(3);
    X2=init(4);
elseif outsiz==5
    S2=init(3);
    X2=init(4);
    S3=init(5);
elseif outsiz==6
    S2=init(3);
    X2=init(4);
    S3=init(5);
    X3=init(6);
end

%Both growth functions take Monod form and the hydrogen inhibition function
%takes the non-competitive form.

%Define Monod growth functions and inhibition function
if outsiz ==3 || outsiz==4 || outsiz==5
    f1=(km1*S1/(Ks1+S1));
elseif outsiz==6
    f1=(km1*S1/(Ks1+S1))*(S3/(Ks3c+S3));
end
f2=(km2*S2/(Ks2+S2));
if outsiz==5 || outsiz==6
    f3=(km3*S3/(Ks3+S3));
end

%Generate system of ODEs for chosen motif
switch motif
    
    case 'fc'
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=[]; eq6=[];
    case 'sc'
        eq1=D*(S1in-S1)-f1*X1-f1*X2;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1; 
        eq3=-D*X2+Y2*f1*X2-kdec2*X2;
        eq4=[]; eq5=[]; eq6=[];
    case 'fci'
        KI3=1e-3; %Default, change these
        I3=1/(1+S3/KI3);
        eq1=D*(S1in-S1)-f1*X1*I3;
        eq2=-D*X1+Y1*f1*X1*I3-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1*I3-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=-D*S3+f2*X2; 
        eq6=[];
    case 'ni'
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=D*(S2in-S2)-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=[]; eq6=[];
    case 'syn'
        I2=1/(1+S2/KI2);
        eq1=D*(S1in-S1)-f1*X1*I2;
        eq2=-D*X1+Y1*f1*X1*I2-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1*I2-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=[]; eq6=[];
    case 'pi'
        I4 = 1/(1+KI2/S2);
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1;
        eq4=-D*X2+Y3*f3*X2*I4-kdec2*X2;
        eq5=D*(S3in-S3)-f3*X2*I4; 
        eq6=[];
    case 'ths'
        I2=1/(1+S3/KI2);
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=D*(S2in-S2)+gamma0*(1-Y1)*f1*X1-f2*X2*I2;
        eq4=-D*X2+Y2*f2*X2*I2-kdec2*X2;
        eq5=D*(S3in-S3)+gamma1*(1-Y2)*f2*X2*I2-f3*X3-gamma2*f1*X1; 
        eq6=-D*X3+Y3*f3*X3-kdec3*X3;
end

%Concatenate equations as output
dP=[eq1;eq2;eq3;eq4;eq5;eq6];


