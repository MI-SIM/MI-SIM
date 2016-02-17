function var=set_variables(arg)
%set_variables - set the variables and parameters symbolically to be used in the model
%to analysis
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: SOFTWARE HOSTED SITE
% September 2015; Last revision: 04-Jan-2016

syms S1 X1 S2 X2 S3 X3 km1 Ks1 Y1 kdec1 km2 Ks2 Y2 kdec2 km3 Ks3 Ks3c Y3 kdec3 KI2 D S1in S2in S3in Ks1a Ks2a Ks3a time gamma0 gamma1 gamma2 S3in 

if arg==1 %Monod
        var=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time Ks1 Ks2 Ks3 Ks3c gamma0 gamma1 gamma2];
elseif arg==2 %Contois...needs amending
        var=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time Ks1a Ks2a Ks3a gamma0];
end

%For further growth functions
% if arg==3
%         var=[kmp Yp kdecp kmH YH kdecH KIH D Spin time kp kH];
% elseif arg==4
%         var=[kmp Yp kdecp kmH YH kdecH KIH D Spin time Ksp KsH N];
% elseif arg==5
%         var=[kmp Yp kdecp kmH YH kdecH KIH D Spin time Ksp KsH Ks0p Ks0H SHin];
% end