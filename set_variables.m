function [var,var_str]=set_variables(growth)
%set_variables - set the variables and parameters symbolically to be used in the model
%to analysis
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: https://github.com/MI-SIM/MI-SIM
% September 2015; Last revision: 23-Feb-2016

syms S1 X1 S2 X2 S3 X3 km1 Ks1 Y1 kdec1 km2 Ks2 Y2 kdec2 km3 Ks3 Ks3c Y3 kdec3 KI2 D S1in S2in S3in Ks1a Ks2a Ks3a time gamma0 gamma1 gamma2 S3in 

switch growth
    case 'Monod'
        var=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time Ks1 Ks2 Ks3 Ks3c gamma0 gamma1 gamma2];
        var_str={'km1', 'Y1', 'kdec1', 'km2', 'Y2', 'kdec2', 'km3','Y3','kdec3','KI2','D','S1in','S2in','S3in','time','Ks1','Ks2','Ks3','Ks3c','gamma0','gamma1','gamma2'};

    case 'Contois'    
        var=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time Ks1a Ks2a Ks3a gamma0]; %TBC
        
    case 'Tessier'
        %TBC
    case 'Moser'
        %TBC
    case 'Logistic'
        %TBC
    case 'Andrews'
        %TBC
    case 'Thermodynamic'
        %TBC
end

