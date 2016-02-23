%growth_functions - generate the system growth functions
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: SOFTWARE HOSTED SITE
% February 2015; Last revision: 19-Feb-2016

%Define growth functions according to selection growth model
adt=1;
switch growth
    case 'Monod'         %Monod growth function
        
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
        f3tx = '';
        
        switch motif
            
            case 'fc'
                
                I2tx = '';
                
            case 'sc'
                
                I2tx = '';
                
            case 'fci'
                
                I2tx = '';
                
            case 'ni'
                
                I2tx = '';
                
            case 'syn'
                
                I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                
            case 'pi'
                
                I2tx = '$I_4 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                
            case 'ths'
                if dse==1
                    adt=S3/(Ks3c+S3);
                end
                f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}\frac{S_3}{K_{S,3c} + S_3}$';
                f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3} + S_3}$';
                I2tx = '$I_2 = \frac{1}{1+\frac{S_3}{K_{i,2}}}$';
        end
        if dse==1
            f1=(km1*S1/(Ks1+S1))*adt;
            f2=(km2*S2/(Ks2+S2));
            f3=(km3*S3/(Ks3+S3));
        end
        
        handles.fnctext = strvcat(f1tx,f2tx,f3tx,I2tx);
        
    case 'Contois'       %Contois growth function
        
        
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{C,1}*X_1 + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{C,2}*X_2 + S_2}$';
        f3tx = '';
        
        switch motif
            
            case 'fc'
                
                I2tx = '';
                
            case 'sc'
                
                I2tx = '';
                
            case 'fci'
                
                I2tx = '';
                
            case 'ni'
                
                I2tx = '';
                
            case 'syn'
                
                I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                
            case 'pi'
                
                I2tx = '$I_4 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                
            case 'ths'
                if dse==1
                    adt=S3/(Ks3c*X3+S3);
                end
                f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{C,1}*X_1 + S_1}\frac{S_3}{K_{C,3c}*X_3 + S_3}$';
                f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{C,3}*X_3 + S_3}$';
                I2tx = '$I_3 = \frac{1}{1+\frac{S_3}{K_{i,3}}}$';
        end
        
        if dse==1
            f1=(km1*S1/(Ks1a*X1+S1))*adt;
            f2=(km2*S2/(Ks2a*X2+S2));
            f3=(km3*S3/(Ks3a*X3+S3));
        end
    case 'Tessier'      %Tessier growth function
        
    case 'Moser'      %Moser growth function
        
    case 'Logistic'      %Logistic growth function
        
    case 'Andrews'      %Andrews growth function
        
    case 'Thermodynamic'    %Thermodynamic function
        
end