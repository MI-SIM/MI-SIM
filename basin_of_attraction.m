%basin_of_attraction - Plot basin of attraction
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: SOFTWARE HOSTED SITE
% September 2015; Last revision: 05-Jan-2016

%Get the values entered in the GUI

km1 = str2double(get(handles.km1_in,'String'));
Y1 = str2double(get(handles.y1_in,'String'));
kdec1 = str2double(get(handles.kdec1_in,'String'));
km2 = str2double(get(handles.km2_in,'String'));
Y2 = str2double(get(handles.y2_in,'String'));
kdec2 = str2double(get(handles.kdec2_in,'String'));
km3 = str2double(get(handles.km3_in,'String'));
Y3 = str2double(get(handles.y3_in,'String'));
kdec3 = str2double(get(handles.kdec3_in,'String'));
S1in = str2double(get(handles.s1in_in,'String'));
S2in = str2double(get(handles.s2in_in,'String'));
S3in = str2double(get(handles.s3in_in,'String'));
KI2 = str2double(get(handles.ki2_in,'String'));
Ks1 = str2double(get(handles.ks1_in,'String'));
Ks2 = str2double(get(handles.ks2_in,'String'));
Ks3 = str2double(get(handles.ks3_in,'String'));
Ks3c = str2double(get(handles.ks32_in,'String'));
Ks1a = str2double(get(handles.ks1a_in,'String'));
Ks2a = str2double(get(handles.ks2a_in,'String'));
Ks3a = str2double(get(handles.ks3a_in,'String'));
gamma0 = str2double(get(handles.gamma0,'String'));
gamma1 = str2double(get(handles.gamma1,'String'));
gamma2 = str2double(get(handles.gamma2,'String'));
D = str2double(get(handles.d_in,'String'));
time1 = str2double(get(handles.time_in,'String'));
S1_init=str2double(get(handles.s1_init,'String'));
X1_init=str2double(get(handles.x1_init,'String'));
S2_init=str2double(get(handles.s2_init,'String'));
X2_init=str2double(get(handles.x2_init,'String'));
S3_init=str2double(get(handles.s3_init,'String'));
X3_init=str2double(get(handles.x3_init,'String'));

%Get the values of the fixed points
[eqs, f1, f2, f3, I2, I3, I4]=define_system_equations(motif,growth);
if growth~=1 && growth~=2
    warning('Please select either Monod or Contois')
    return
else
end

%Plug in the values for numeric fixed points
variables=set_variables(growth);
if growth==1
    parameters=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time1 Ks1 Ks2 Ks3 Ks3c gamma0 gamma1 gamma2]; %Add KI3
elseif growth==2
    parameters=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in s2in S3in time1 Ks1a Ks2a Ks3a Ks3c gamma0 gamma1 gamma2];
end
%Check number of equations
noeq=length(eqs);

for k=1:noeq
    eval(['eq',num2str(k),'_numerical=subs(eqs(',num2str(k),'),variables,parameters);']);
end
%Solve system of ODEs
S1_fix_num=[]; X1_fix_num=[]; S2_fix_num=[]; X2_fix_num=[]; S3_fix_num=[]; X3_fix_num=[];
syms S1 X1 S2 X2 S3 X3
if noeq==3
    [S1_fix_num,X1_fix_num,X2_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0]);
elseif noeq==4
    [S1_fix_num,S2_fix_num,X1_fix_num,X2_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0]);
elseif noeq==5
    [S1_fix_num,X1_fix_num,S2_fix_num,X2_fix_num,S3_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0]);
elseif noeq==6
    [S1_fix_num,X1_fix_num,S2_fix_num,X2_fix_num,S3_fix_num,X3_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0 eq6_numerical],S1,X1,S2,X2,S3,X3);
end

fixed_numerical=[S1_fix_num,X1_fix_num,S2_fix_num,X2_fix_num,S3_fix_num,X3_fix_num];
fixed_numerical=double(fixed_numerical);
%Remove invalid fixed points (i.e. those with negative values)
for i=1:length(fixed_numerical(1,:))
    fixed_numerical(fixed_numerical(:,i)<0,:) = [];
end
format shortEng
format compact
fixed_numerical1=fixed_numerical;

%Define parameter sweep
minp1a = get(handles.min_p1,'String');
minp1 = str2double(minp1a);
maxp1a = get(handles.max_p1,'String');
maxp1 = str2double(maxp1a);
stepp1a = get(handles.step_p1,'String');
stepp1 = str2double(stepp1a);
low_lim=minp1;
up_lim=maxp1;

%Vary the initial condition
x=linspace(low_lim,up_lim,stepp1);
y=linspace(low_lim,up_lim,stepp1);

%Label the axis
xvarls = get(handles.simparam1,'String'); yvarls=get(handles.simparam2,'String');
xvarop = get(handles.simparam1,'Value'); yvarop = get(handles.simparam2,'Value');

xvar = xvarls(xvarop,:); yvar = yvarls(yvarop,:);
axes(handles.trajectoryplot); cla
xlabel(xvar,'FontSize',18)
ylabel(yvar,'FontSize',18)

%Set up progress bar
total_steps = length(x)*length(y);
no_steps=zeros(length(x),length(y));

set(handles.func_prog,'String','Running: Basin of Attraction','ForegroundColor','r')

solver=char(strtrim(solver));
h=handles.timestamp; h1=[]; tt=[];

z_boa=zeros(length(x),length(y));
initial=zeros(length(x),length(y),2);  %set up a matrix of zeros
for i=1:length(x)
    for j=1:length(y)
        
        initial(i,j,1)=x(i);    %fill in the matrix with values
        initial(i,j,2)=y(j);    %fill in the matrix with values
        
        % varies two of the initial conditions and lets the user set the others
        if noeq==3
            if strcmp(xvar,'S1') && strcmp(yvar,'S1')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S1')
                S1i=initial(i,j,1);
            elseif strcmp(yvar,'S1')
                S1i=initial(i,j,2);
            else
                S1i=S1_init;
            end
            
            if strcmp(xvar,'X1') && strcmp(yvar,'X1')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X1')
                X1i=initial(i,j,1);
            elseif strcmp(yvar,'X1')
                X1i=initial(i,j,2);
            else
                X1i=X1_init;
            end
            
            if strcmp(xvar,'X2') && strcmp(yvar,'X2')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X2')
                X2i=initial(i,j,1);
            elseif strcmp(yvar,'X2')
                X2i=initial(i,j,2);
            else
                X2i=X2_init;
            end
            
            init1=[S1i, X1i, X2i];
            
        elseif noeq==4
            
            if strcmp(xvar,'S1') && strcmp(yvar,'S1')
                msgbox('Select two distinct variables','Error','error')
                %return
                
            elseif strcmp(xvar,'S1')
                S1i=initial(i,j,1);
            elseif strcmp(yvar,'S1')
                S1i=initial(i,j,2);
            else
                S1i=S1_init;
            end
            
            if strcmp(xvar,'S2') && strcmp(yvar,'S2')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S2')
                S2i=initial(i,j,1);
            elseif strcmp(yvar,'S2')
                S2i=initial(i,j,2);
            else
                S2i=S2_init;
            end
            
            if strcmp(xvar,'X1') && strcmp(yvar,'X1')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X1')
                X1i=initial(i,j,1);
            elseif strcmp(yvar,'X1')
                X1i=initial(i,j,2);
            else
                X1i=X1_init;
            end
            
            if strcmp(xvar,'X2') && strcmp(yvar,'X2')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X2')
                X2i=initial(i,j,1);
            elseif strcmp(yvar,'X2')
                X2i=initial(i,j,2);
            else
                X2i=X2_init;
            end
            
            init1=[S1i, X1i, S2i, X2i];
            
            
            
        elseif noeq==5
            
            if strcmp(xvar,'S1') && strcmp(yvar,'S1')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S1')
                S1i=initial(i,j,1);
            elseif strcmp(yvar,'S1')
                S1i=initial(i,j,2);
            else
                S1i=S1_init;
            end
            
            if strcmp(xvar,'S2') && strcmp(yvar,'S2')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S2')
                S2i=initial(i,j,1);
            elseif strcmp(yvar,'S2')
                S2i=initial(i,j,2);
            else
                S2i=S2_init;
            end
            
            if strcmp(xvar,'X1') && strcmp(yvar,'X1')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X1')
                X1i=initial(i,j,1);
            elseif strcmp(yvar,'X1')
                X1i=initial(i,j,2);
            else
                X1i=X1_init;
            end
            
            if strcmp(xvar,'X2') && strcmp(yvar,'X2')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X2')
                X2i=initial(i,j,1);
            elseif strcmp(yvar,'X2')
                X2i=initial(i,j,2);
            else
                X2i=X2_init;
            end
            
            if strcmp(xvar,'S3') && strcmp(yvar,'S3')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S3')
                S3i=initial(i,j,1);
            elseif strcmp(yvar,'S3')
                S3i=initial(i,j,2);
            else
                S3i=S3_init;
            end
            
            
            
            init1=[S1i, X1i, S2i, X2i, S3i];
        elseif noeq==6
            
            if strcmp(xvar,'S1') && strcmp(yvar,'S1')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S1')
                S1i=initial(i,j,1);
            elseif strcmp(yvar,'S1')
                S1i=initial(i,j,2);
            else
                S1i=S1_init;
            end
            
            if strcmp(xvar,'S2') && strcmp(yvar,'S2')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S2')
                S2i=initial(i,j,1);
            elseif strcmp(yvar,'S2')
                S2i=initial(i,j,2);
            else
                S2i=S2_init;
            end
            
            if strcmp(xvar,'X1') && strcmp(yvar,'X1')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X1')
                X1i=initial(i,j,1);
            elseif strcmp(yvar,'X1')
                X1i=initial(i,j,2);
            else
                X1i=X1_init;
            end
            
            if strcmp(xvar,'X2') && strcmp(yvar,'X2')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X2')
                X2i=initial(i,j,1);
            elseif strcmp(yvar,'X2')
                X2i=initial(i,j,2);
            else
                X2i=X2_init;
            end
            
            if strcmp(xvar,'S3') && strcmp(yvar,'S3')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S3')
                S3i=initial(i,j,1);
            elseif strcmp(yvar,'S3')
                S3i=initial(i,j,2);
            else
                S3i=S3_init;
            end
            
            if strcmp(xvar,'X3') && strcmp(yvar,'X3')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X3')
                X3i=initial(i,j,1);
            elseif strcmp(yvar,'X3')
                X3i=initial(i,j,2);
            else
                X3i=X3_init;
            end
            
            init1=[S1i, X1i, S2i, X2i, S3i, X3i];
        end
        
        %stops the ode solver once we have reached a steady state
        %options1=odeset('Events',@eventfun);        %eventfun(tout,yout,fixed_numerical)
        %Options for ODE solver
        options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[]);
        flag = 1;
        cno=[];
        warning off
        if growth==1
            
         eval(['[tout,yout]=',solver,'(@model_gen, [0:0.01:time1], init1, options, S1in, D, Y1, kdec1, Y2, kdec2, Y3, kdec3, km1, Ks1, km2, Ks2, km3, Ks3, Ks3c, KI2,gamma0, gamma1, gamma2, S2in, S3in,h,h1,tt,motif,flag,cno);']);
                    
        elseif growth==2
            %[tout,yout]=ode23s(@four_mod2, 0:0.01:time1, init1, options, Spin1, D1, Yp1, kdecp1, YH1, kdecH1, kmp1, kmH1, Ksxp1, KsxH1, KIH1);
        end
        warning on
        
        iter=j+length(y)*(i-1);
        total_tim=(iter)*100/(length(x)*length(y));
        set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
        %pause(0.000000000000000001)
        
        %choose error
        error=str2double(get(handles.error_val,'String'));
        if error<0 || isnumeric(error)==0
            set(handles.error_val,'String',1e-4);
            error=str2double(get(handles.error_val,'String'));
            msgbox('The value entered for Error was not valid so the default value wa reset','Error','error')
        else
        end
        
        %Get basin of attraction
        for k=1:size(fixed_numerical,1)
            if norm(yout(end,:)-fixed_numerical(k,:))<error
                z_boa(i,j)=k-1;
            end
        end

    end
end
%plot the basin of attraction
axis([low_lim up_lim low_lim up_lim]);
colormap(hot); warning off
if mean(std(z_boa'))==0
else
    contourf(x,y,z_boa','linestyle','none','levelstep',1)
end

hold on

for k=1:size(fixed_numerical,1)
    %Label plot
    eval(['zons',num2str(k-1),'=zeros(length(x)*length(y),1);']);
    
    %Find indices
    eval(['zzer',num2str(k-1),'=find(z_boa==',num2str(k-1),');']);
    eval(['zons',num2str(k-1),'(zzer',num2str(k-1),')=1;']);
    
    %Reshape
    eval(['rM',num2str(k-1),'=reshape(zons',num2str(k-1),',length(x),length(y));']);
    
    xl=0:length(x)-1;
    yl=0:length(y)-1;
    
    [Xl Yl]=meshgrid(xl,yl);
    
    eval(['cX',num2str(k-1),'=mean(Xl(rM',num2str(k-1),'==1));']);
    eval(['cY',num2str(k-1),'=mean(Yl(rM',num2str(k-1),'==1));']);
    
    set(handles.func_prog,'String','Completed: Basin of Attraction','ForegroundColor',[0 0.6 1])
    
    eval(['ccX',num2str(k-1),'=cX',num2str(k-1),'*max(x)/length(x);']);
    eval(['ccY',num2str(k-1),'=cY',num2str(k-1),'*max(y)/length(y);']);
    
    eval(['text(ccX',num2str(k-1),',ccY',num2str(k-1),',''FP',num2str(k),''',''fontsize'',14,''color'',''w'',''backgroundcolor'',''k'',''fontweight'',''bold'');']);
    
end

%Save Report
newfig_l=figure;
axesobject_l=copyobj(handles.trajectoryplot,newfig_l);
title('Basin of Attraction')
export_fig temp_fig/ExpandReport -pdf -r600 -append
close(newfig_l)
