%run_script_gui - m-file for running the selected analysis routine [via the
%RUN button]
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: SOFTWARE HOSTED SITE
% September 2015; Last revision: 04-Jan-2016

%% Preamble
%Enable reporting
set(handles.report_menu,'Enable','on');
set(handles.reportpanel,'Visible','on'); set(handles.email_add,'Enable','on'); set(handles.server_add,'Enable','on')
set(handles.email_txt,'Enable','on'); set(handles.serv_txt,'Enable','on');

%Get value indicating if [jacobian sparsity matrix; normalised eigenvalues;
%plot overlay] will be used in the solver
jacobanaly=get(handles.jacobian_but,'Value');
ovra=get(handles.normeig,'Value');
ovro=get(handles.overlay,'Value');


if handles.mpltopt==0
    %Clear figure axes and text information
    if ovra==0
        if isfield(handles,'solutionplot')
            cla(handles.solutionplot)
        elseif isfield(handles,'hsub')
            for k=1:length(handles.hsub)
                cla(handles.hsub(k))
            end
        end
    end
    
    if ovro==0
        cla(handles.trajectoryplot)
    end
    if exist('handles.htxt')
        delete(handles.htxt)
    end
    if exist('handles.stbtxt')
        delete(handles.stbtxt)
    end
    %Define system equations
    [eqs, f1, f2, f3, I2, I3, I4]=define_system_equations(motif,growth);
    
    %Check that an existing modelled growth function is selected
    if growth~=1 && growth~=2
        warning('Please select either Monod or Contois')
        return
    else
    end
end

%% Algorithms
%Case structure for analysis routines
switch handles.simtype
    case 'single_p' %Single-point analysis
        
        %Get the symbolic variables and parameters
        variables=set_variables(growth);
        if growth==1
            parameters=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time1 Ks1 Ks2 Ks3 Ks3c gamma0 gamma1 gamma2];
        elseif growth==2
            parameters=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in s2in S3in time1 Ks1a Ks2a Ks3a Ks3c gamma0 gamma1 gamma2];
        end
        
        %Check number of equations
        noeq=length(eqs);
        
        %REPLACE HERE WITH MATLABFUNCTION
        for k=1:noeq
            eval(['eq',num2str(k),'_numerical=subs(eqs(',num2str(k),'),variables,parameters);']);
        end
        
        %% Fixed points solutions
        S1_fix_num=[]; X1_fix_num=[]; S2_fix_num=[]; X2_fix_num=[]; S3_fix_num=[]; X3_fix_num=[];
        
        if noeq==3
            [S1_fix_num,X1_fix_num,X2_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0]);
            set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
            set(handles.s2_check,'Enable','off','value',0); set(handles.x2_check,'Enable','on','value',1);
            set(handles.s3_check,'Enable','off','value',0); set(handles.x3_check,'Enable','off','value',0);
        elseif noeq==4
            [S1_fix_num,S2_fix_num,X1_fix_num,X2_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0]);
            set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
            set(handles.s2_check,'Enable','on','value',1); set(handles.x2_check,'Enable','on','value',1);
            set(handles.s3_check,'Enable','off','value',0); set(handles.x3_check,'Enable','off','value',0);
        elseif noeq==5
            [S1_fix_num,X1_fix_num,S2_fix_num,X2_fix_num,S3_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0]);
            set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
            set(handles.s2_check,'Enable','on','value',1); set(handles.x2_check,'Enable','on','value',1);
            set(handles.s3_check,'Enable','on','value',1); set(handles.x3_check,'Enable','off','value',0);
        elseif noeq==6
            syms S1 X1 S2 X2 S3 X3
            assume(S1>=0); assume(X1>=0); assume(S2>=0); assume(X2>=0); assume(S3>=0); assume(X3>=0); %Make assumptions (X,S real non-negative)
            [S1_fix_num,X1_fix_num,S2_fix_num,X2_fix_num,S3_fix_num,X3_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0 eq6_numerical==0],[S1,X1,S2,X2,S3,X3],'real',true);
            
            set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
            set(handles.s2_check,'Enable','on','value',1); set(handles.x2_check,'Enable','on','value',1);
            set(handles.s3_check,'Enable','on','value',1); set(handles.x3_check,'Enable','on','value',1);
        end
        
        fixed_numerical=[S1_fix_num,X1_fix_num,S2_fix_num,X2_fix_num,S3_fix_num,X3_fix_num];
        f = double(fixed_numerical);
        
        %Remove invalid fixed points (fp<0)
        [ii,jj]=find(f<0);
        iu = unique(ii);
        f(iu,:)=[];
        
        format short
        format compact
        fixed_numerical1=f;
        
        %Set the fixed points in the GUI        
        if noeq==3
            set(handles.fixed_points_s1,'String',double(fixed_numerical1(:,1)));
            set(handles.fixed_points_x2,'String',double(fixed_numerical1(:,2)));
            set(handles.fixed_points_x1,'String',double(fixed_numerical1(:,3)));
            S2_init=[];
        elseif noeq==4
            set(handles.fixed_points_s1,'String',double(fixed_numerical1(:,1)));
            set(handles.fixed_points_x1,'String',double(fixed_numerical1(:,2)));
            set(handles.fixed_points_s2,'String',double(fixed_numerical1(:,3)));
            set(handles.fixed_points_x2,'String',double(fixed_numerical1(:,4)));
        elseif noeq==5
            set(handles.fixed_points_s1,'String',double(fixed_numerical1(:,1)));
            set(handles.fixed_points_x1,'String',double(fixed_numerical1(:,2)));
            set(handles.fixed_points_s2,'String',double(fixed_numerical1(:,3)));
            set(handles.fixed_points_x2,'String',double(fixed_numerical1(:,4)));
            set(handles.fixed_points_s3,'String',double(fixed_numerical1(:,5)));
        elseif noeq==6
            set(handles.fixed_points_s1,'String',double(fixed_numerical1(:,1)));
            set(handles.fixed_points_x1,'String',double(fixed_numerical1(:,2)));
            set(handles.fixed_points_s2,'String',double(fixed_numerical1(:,3)));
            set(handles.fixed_points_x2,'String',double(fixed_numerical1(:,4)));
            set(handles.fixed_points_s3,'String',double(fixed_numerical1(:,5)));
            set(handles.fixed_points_x3,'String',double(fixed_numerical1(:,6)));
        end
        
        %Set labels
        lfn=length(fixed_numerical(:,1));
        if lfn<=5
            for k=1:lfn
                eval(['set(handles.fp',num2str(k),',','''Visible'',''on'')'])
            end
        else
            for k=1:lfn
                eval(['set(handles.fp',num2str(k),',','''Visible'',''on'')'])
            end
        end
        
        %% Stability
        %Shows the stability of the fixed points in the GUI
        stability_analysis
        set(handles.fixed_points_stability,'String',s');
        set(handles.fixed_points_stability_rh,'String',s2');
        
        %Options for ODE solver
        if jacobanaly==1
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[],'Jpattern',Jpat);
        else
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[]);
        end
        clear init
        %Set initial conditions
        if noeq==3
            init2=[S1_init, X1_init, X2_init];
        elseif noeq==4
            init2=[S1_init, X1_init, S2_init, X2_init];
        elseif noeq==5
            init2=[S1_init, X1_init, S2_init, X2_init, S3_init];
        elseif noeq==6
            init2=[S1_init, X1_init, S2_init, X2_init, S3_init, X3_init];
        end
        
        solver=char(strtrim(solver));
        
        %% Dynamics
        %Run solver - Dynamics
        if growth==1
            set(handles.func_prog,'String','Running: Dynamics','ForegroundColor','r')
            h=handles.timestamp;
            h1=handles.progress;
            tt=time1; flag=[]; cno=[];
            te = cputime;
            eval(['[tout,yout]=',solver,'(@model_gen, [0:0.01:time1], init2, options, S1in, D, Y1, kdec1, Y2, kdec2, Y3, kdec3, km1, Ks1, km2, Ks2, km3, Ks3, Ks3c, KI2,gamma0, gamma1, gamma2, S2in, S3in,h,h1,tt,motif,flag,cno);']);
            tfe= cputime-te;
        elseif growth==2
            eval(['[tout,yout]=ode23s(@four_mod2, [0:0.01:time1], init2, options, S1in, D, Y1, kdec1, Y2, kdec2, Y3, kdec3, km1, Ks1a, km2, Ks2a, km3, Ks3a, KI2,gamma0,S3in);']);
        end
        set(handles.func_prog,'String','Completed: Dynamics','ForegroundColor',[0 0.6 1])
        handles.tout=tout;
        handles.yout=yout;
        handles.fn=fixed_numerical1;
        %Plot results
        which_plot='sol';
        plot_results;
        set(handles.plot_button,'enable','on')
        set(handles.plot_button2,'enable','on')
        set(handles.twodphase,'Enable','on'); set(handles.threedphase,'Enable','on')
        
    case 'multiple_p' %Multiple-point analysis
        if handles.mpltopt==0 %%%CHECK
            %Plug in the values for numerical fixed points, ommitting selected
            %parameters (i.e. the parameter pair to sweep)
            
            variables = symvar(eqs);
            str1=get(handles.simparam1,'String'); str2=get(handles.simparam2,'String');
            val1=get(handles.simparam1,'Value'); val2=get(handles.simparam2,'Value');
            param1=str1(val1,:); param2=str2(val2,:);
            str1([val1,val2],:)=[];
            str1=strvcat(str1,'gamma0');
            
            %Check number of equations and additional gammas for EQS==6
            noeq=length(eqs);
            if noeq==6
                str1=strvcat(str1,'gamma1','gamma2');
            end
            
            %Match indices of symbolic parameters (variables) and
            %numerical parameters
            ind=1; str1_cell=cellstr(str1);
            for k=1:length(variables)
                if ~isempty(find(ismember(str1_cell,char(variables(k)))))
                    gind(ind)=find(ismember(str1_cell,char(variables(k))));
                    pind(ind)=k;
                    ind=ind+1;
                end
            end
                        
            variables=variables(pind);
            str1=str1(gind,:);
            
            %Evaluate parameters to numeric
            for k=1:length(str1)
                parameters(k)=eval(str1(k,:));
            end
            
            for k=1:noeq
                eval(['eq',num2str(k),'_numerical=subs(eqs(',num2str(k),'),variables,parameters);']);
            end
            
            % Set fixed points empty 
            S1_fix_num=[]; X1_fix_num=[]; S2_fix_num=[]; X2_fix_num=[]; S3_fix_num=[]; X3_fix_num=[];
            
            %% Solutions
            % Three ODEs
            if noeq==3
                %Solve ODEs
                [S1_fix_num,X1_fix_num,X2_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0],'S1','X1','X2');
                set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
                set(handles.s2_check,'Enable','off','value',0); set(handles.x2_check,'Enable','on','value',1);
                set(handles.s3_check,'Enable','off','value',0); set(handles.x3_check,'Enable','off','value',0);
                %Substitute operating parameters
                r1=str2double(get(handles.min_p1,'String')); r2=str2double(get(handles.max_p1,'String')); s1a=str2double(get(handles.step_p1,'String'));
                r3=str2double(get(handles.min_p2,'String')); r4=str2double(get(handles.max_p2,'String')); s2a=str2double(get(handles.step_p2,'String'));
                %Get parameter sweep
                range1=linspace(r1,r2,s1a);
                range2=linspace(r3,r4,s2a);
                %Create matlab functions of fixed point solutions (reduces
                %compute time)
                S1n = matlabFunction(S1_fix_num); X1n = matlabFunction(X1_fix_num); X2n = matlabFunction(X2_fix_num);
                lfp=length(S1_fix_num);
                cn=1;
                %Get numerical fixed-point solutions across parameter space
                for k=range1
                    for kk=range2
                        sss1(1:lfp,cn)=S1n(kk,k);
                        ssx1(1:lfp,cn)=X1n(kk,k);
                        ssx2(1:lfp,cn)=X2n(kk,k);
                        cn=cn+1;
                    end
                end
                S1new=sss1;
                X1new=ssx1;
                X2new=ssx2;
                        
                %% Determine steady-states for each point
                
                %Convert valid biomass (X) to numeric, invalid to NaN
                ssx1(ssx1>0)=5.5; ssx1(ssx1<0)=nan;
                ssx2(ssx2>0)=6.5; ssx2(ssx2<0)=nan;
                
                %Check for existence for remaining variable (substrate, S)
                sss1(sss1<0)=nan; sss1(sss1>0)=0;
                
                %Reshape fixed point solutions at each FP and assign
                %steady-state label (SS)
                for ks=1:size(ssx1,1);
                    ssxr1{ks}=reshape(ssx1(ks,:),s1a,s2a);
                    ssxr2{ks}=reshape(ssx2(ks,:),s1a,s2a);
                    sssr1{ks}=reshape(sss1(ks,:),s1a,s2a);
                    ssx{ks}=ssxr1{ks}+ssxr2{ks}+sssr1{ks};
                    ssx{ks}(ssx{ks}==0)=1; %SS1
                    ssx{ks}(ssx{ks}==5.5)=2; %SS2
                    ssx{ks}(ssx{ks}==6.5)=4; %SS4
                    ssx{ks}(ssx{ks}==12)=3; %SS3
                    ssx{ks}(isnan(ssx{ks}))=0; %Does not exist
                end
                handles.ssxr=ssx;
                %Contour plot of Steady-states
                axes(handles.trajectoryplot)
                fst=handles.ssxr;
                final_exist=sum(cat(3,fst{:}),3); %Sum SS to give unique region
                contourf(range1,range2,final_exist,'LineStyle','none','LevelStep',1)
                set(handles.colormp,'Visible','on');
                xlabel(param1)
                ylabel(param2)
                
                %Progress bar
                total_steps=ncols;
                
                set(handles.func_prog,'String','Running: Creating matrix of FP solutions','ForegroundColor','r')
                
                %Get Jacobian of system
                sym_vars = symvar(eqs);
                var_new=[]; var_new_char=[];
                for k=1:length(sym_vars)
                    if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'X2')
                        var_new=var_new;
                        var_new_char=var_new_char;
                    else
                        var_new=[var_new;sym_vars(k)];
                        var_new_char=strvcat(var_new_char,char(sym_vars(k)));
                    end
                end
                for kk=1:length(var_new_char)
                    eqs=subs(eqs,char(var_new(kk)),subs(var_new_char(kk,:)));
                end
                
                var_order=strvcat('S1','X1','X2');
                %Calculate Jacobian matrix
                jac_sys=jacobian(vpa(eqs),sym({'S1','X1','X2'}));
                
                %Calculate characteristic polynomial
                cpoly=charpoly(jac_sys);
                [nrows,ncols]=size(S1new);
                
                %for the progess bar
                total_steps=numel(S1new);
                no_steps=zeros(nrows,ncols);
                
                set(handles.func_prog,'String','Running: Multiple operating points','ForegroundColor','r')
                axes(handles.spmatrix)
                
                xln=repmat(1:s1a,1,s2a);
                ylm=repmat(1:s1a,s2a,1);
                yln=reshape(ylm,1,s1a*s2a);
                
                if exist('htxt')
                    delete(htxt)
                end
                
                %Create Matlab Functions of the Jacobian and characteristic
                %polynomial
                jacs=matlabFunction(jac_sys);
                cpl=matlabFunction(cpoly);
                
                %Stability calculation
                for ii=1:nrows
                    for kk=1:ncols
                        
                        %Update progress information
                        iter=kk+ncols*(ii-1);
                        total_tim=(iter)*100/(total_steps);
                        set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                        c = clock; % [year month day hour minute seconds]
                        TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                        set(handles.timestamp,'String',['Time: ',TSMP])
                        
                        %Update sparsity matrix
                        jac_P=jacs(S1new(ii,kk),X1new(ii,kk),X2new(ii,kk));
                        colormap(hot)
                        xtxt=[param1,' = ',num2str(range1(yln(kk)))];
                        ytxt=[param2,' = ',num2str(range2(xln(kk)))];
                        fpx=['Fixed Point ',num2str(ii)];
                        
                        if isnan(jac_P)==0
                            minm=min(min(double(jac_P)));
                            maxm=max(max(double(jac_P)));
                            imagesc(double(jac_P),[minm maxm])
                            text(-2.5,2.5,strvcat('Jacobian Matrix',xtxt,ytxt,fpx),'interpreter','latex','horiz','left','vert','middle','fontsize',14);
                        else
                            imagesc(ones(4,4),[0 1])
                            text(1.5,2.5,'No solution')
                            text(-2.5,2.5,strvcat('Jacobian Matrix',xtxt,ytxt,fpx),'interpreter','latex','horiz','left','vert','middle','fontsize',14);
                        end
                        set(gca,'Xticklabel',cellstr(var_order),'Yticklabel',cellstr(var_order))
                        pause(0.000000000000000001)
                        
                        %Get numerical characteristic polynomial
                        if S1new(ii,kk)>=0 && X1new(ii,kk) >= 0 && X2new(ii,kk) >= 0
                            cpoly=cpl(S1new(ii,kk),X1new(ii,kk),X2new(ii,kk));
                            
                            %Calculate eigenvalues
                            eigen=roots(cpoly);
                            
                            %Determine stability
                            if all(real(eigen<0))==1 && X1new(ii,kk)==0 && X2new(ii,kk)==0
                                stz(ii,kk)=1; sta(ii,kk)=1;
                            elseif all(real(eigen<0))==1 && X1new(ii,kk)>=0 && X2new(ii,kk)==0
                                stz(ii,kk)=2; sta(ii,kk)=1;
                            elseif all(real(eigen<0))==1 && X1new(ii,kk)>=0 && X2new(ii,kk)>=0
                                stz(ii,kk)=3; sta(ii,kk)=1;
                            elseif all(real(eigen<0))==0 && X1new(ii,kk)==0 && X2new(ii,kk)==0
                                stz(ii,kk)=1; sta(ii,kk)=0.1;
                            elseif all(real(eigen<0))==0 && X1new(ii,kk)>=0 && X2new(ii,kk)==0
                                stz(ii,kk)=2; sta(ii,kk)=0.2;
                            elseif all(real(eigen<0))==0 && X1new(ii,kk)>=0 && X2new(ii,kk)>=0
                                stz(ii,kk)=3; sta(ii,kk)=0.3;
                            end
                        else
                            stz(ii,kk)=0; sta(ii,kk)=0;
                        end
                        
                        
                    end
                end
                
                %Determine which points in SS are stable
                [ln,wd]=size(fst{1});
                [nr,nc]=size(stz);
                
                for km=1:nr
                    snr{km}=reshape(stz(km,:),ln,wd);
                    sna{km}=reshape(sta(km,:),ln,wd);
                end
                
                z1=1:ln;
                z2=1:wd;
                
                for km=1:nr
                    b=1;
                    combined_s{km}=snr{km};
                    combined_a{km}=sna{km};
                    for kk=z1
                        for kj=z2
                            sz(km,b)=combined_s{km}(kk,kj);
                            sza(km,b)=combined_a{km}(kk,kj);
                            b=b+1;
                        end
                    end
                end
                %Assign unique regions combining SS & stability information
                for kk=1:length(sz)
                    comb_s=[sz(:,kk),sza(:,kk)];
                    und=unique(comb_s,'rows');
                    unind=find(und(:,1)>0);
                    ssund=und(unind,:);
                    
                    if all(ssund(:,1)==1)
                        Js(kk)=1;
                    elseif all(ssund(:,1)==2)
                        Js(kk)=2;
                    elseif all(ssund(:,1)==3)
                        Js(kk)=3;
                    elseif any(ssund(:,1)==1) && any(ssund(:,1)==2) && ~ismember(3,ssund(:,1))
                        Js(kk)=4;
                    elseif any(ssund(:,1)==1) && any(ssund(:,1)==3) && ~ismember(2,ssund(:,1))
                        Js(kk)=5;
                    elseif any(ssund(:,1)==1) && any(ssund(:,1)==2) && any(ssund(:,1)==3)
                        Js(kk)=6;
                    elseif any(ssund(:,1)==2) && any(ssund(:,1)==3) && ~ismember(1,ssund(:,1))
                        Js(kk)=7;
                    end
                    
                    if all(ssund(:,2)==1)
                        Jsf(kk)=Js(kk);
                    elseif any(ssund(:,2)==0.1)
                        funstab=find(ssund(:,2)==0.1);
                        unss=ssund(funstab,1);
                        if length(unss)==1 || length(unss)==2
                            Jsf(kk)=Js(kk)+0.51;
                        elseif length(unss)==3
                            Jsf(kk)=0;
                        end
                    elseif any(ssund(:,2)==0.2)
                        funstab=find(ssund(:,2)==0.2);
                        unss=ssund(funstab,1);
                        if length(unss)==1 || length(unss)==2
                            Jsf(kk)=Js(kk)+0.52;
                        elseif length(unss)==3
                            Jsf(kk)=0;
                        end
                        
                    elseif any(ssund(:,2)==0.3)
                        funstab=find(ssund(:,2)==0.3);
                        unss=ssund(funstab,1);
                        if length(unss)==1 || length(unss)==2
                            Jsf(kk)=Js(kk)+0.53;
                        elseif length(unss)==3
                            Jsf(kk)=0;
                        end
                        
                    end
                    
                end
                
                Jj=reshape(Jsf,ln,wd);
                
                final_stable=Jj';
                
                %Store values for plotting %REMOVE THIS AND REPLACE WITH
                %LABELLING
                handles.stzp=final_stable; handles.ssts=floor(Js');
                handles.fe=final_exist;
                handles.range1=range1;
                handles.range2=range2; handles.siz1=s1a; handles.siz2=s2a; handles.prm1=param1; handles.prm2=param2;
                bifstr=strvcat('Select fixed-point to view...','Steady-states','Stability'); %REMOVE THIS
                set(handles.bifursp,'Visible','on','String',bifstr)
                set(handles.func_prog,'String','Completed: Multiple operating points','ForegroundColor',[0 0.6 1])
                guidata(hObject,handles);
                
            %Four equations
            elseif noeq==4
                warning off
                set(handles.progress,'String',['Progress: ',num2str(0),'%'])
                syms S1 X1 S2 X2 S3 X3 D S1in %Specify operating parameters rather than D S1in
                
                %Make assumptions
                assume(S1>=0); assume(X1>=0); assume(S2>=0); assume(X2>=0); assume(S3>=0); assume(X3>=0);
                set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
                set(handles.s2_check,'Enable','on','value',1); set(handles.x2_check,'Enable','on','value',1);
                set(handles.s3_check,'Enable','off','value',0); set(handles.x3_check,'Enable','off','value',0);
                %Substitute operating parameters
                r1=str2double(get(handles.min_p1,'String')); r2=str2double(get(handles.max_p1,'String')); s1a=str2double(get(handles.step_p1,'String'));
                r3=str2double(get(handles.min_p2,'String')); r4=str2double(get(handles.max_p2,'String')); s2a=str2double(get(handles.step_p2,'String'));
                range1=linspace(r1,r2,s1a);
                range2=linspace(r3,r4,s2a);
                
                %Solve equations at equilibrium
                ssol=get(handles.solmeth_select,'String');
                vsol=get(handles.solmeth_select,'Value');
                calculationType=char(ssol(vsol,:));
                switch calculationType
                    case 'Matlab'
                        
                        [S1_fix_num,X1_fix_num,S2_fix_num,X2_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0],'S1','X1','S2','X2','Real',true);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------------------------
                        
                        need to get rid of unrealistic fixed points (i.e. those with negative values)
                        [~,ncols]=size(sss1);
                        
                        %Determine steady-states for each point
                        
                        %Convert positives to 1, negatives to 0
                        ssx1(ssx1>0)=5.5; ssx1(ssx1<0)=nan;
                        ssx2(ssx2>0)=6.5; ssx2(ssx2<0)=nan;
                        
                        %Check for existence
                        sss1(sss1<0)=nan; sss1(sss1>0)=0;
                        sss2(sss2<0)=nan; sss2(sss2>0)=0;
                        
                        %Create handle of matrices
                        for ks=1:size(ssx1,1);
                            ssxr1{ks}=reshape(ssx1(ks,:),s1a,s2a);
                            ssxr2{ks}=reshape(ssx2(ks,:),s1a,s2a);
                            sssr1{ks}=reshape(sss1(ks,:),s1a,s2a);
                            sssr2{ks}=reshape(sss2(ks,:),s1a,s2a);
                            ssx{ks}=ssxr1{ks}+ssxr2{ks}+sssr1{ks}+sssr2{ks};
                            ssx{ks}(ssx{ks}==0)=1; %SS1
                            ssx{ks}(ssx{ks}==5.5)=2; %SS2
                            ssx{ks}(ssx{ks}==6.5)=4; %SS4
                            ssx{ks}(ssx{ks}==12)=3; %SS3
                            ssx{ks}(isnan(ssx{ks}))=0; %Does not exist
                        end
                        handles.ssxr=ssx;
                        axes(handles.trajectoryplot)
                        fst=handles.ssxr;
                        final_exist=sum(cat(3,fst{:}),3);
                        contourf(range1,range2,final_exist,'LineStyle','none','LevelStep',1)
                        set(handles.colormp,'Visible','on');
                        %LABELS
                        xlabel(param1)
                        ylabel(param2)
                        
                        set(handles.func_prog,'String','Running: Creating matrix of FP solutions','ForegroundColor','r')
                        
                        %Get Jacobian of system
                        sym_vars = symvar(eqs);
                        var_new=[]; var_new_char=[];
                        for k=1:length(sym_vars)
                            if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'S2') || strcmp(char(sym_vars(k)),'X2')
                                var_new=var_new;
                                var_new_char=var_new_char;
                            else
                                var_new=[var_new;sym_vars(k)];
                                var_new_char=strvcat(var_new_char,char(sym_vars(k)));
                            end
                        end
                        for kk=1:length(var_new_char)
                            eqs=subs(eqs,char(var_new(kk)),subs(var_new_char(kk,:)));
                        end
                        
                        var_order=strvcat('S1','X1','S2','X2');
                        jac_sys=jacobian(vpa(eqs),sym({'S1','X1','S2','X2'}));
                        
                        cpoly=charpoly(jac_sys);
                        [nrows,ncols]=size(S1new);
                        
                        %for the progess bar
                        total_steps=numel(S1new);
                        no_steps=zeros(nrows,ncols);
                        
                        set(handles.func_prog,'String','Running: Multiple operating points','ForegroundColor','r')
                        axes(handles.spmatrix)
                        
                        xln=repmat(1:s1a,1,s2a);
                        ylm=repmat(1:s1a,s2a,1);
                        yln=reshape(ylm,1,s1a*s2a);
                        
                        if exist('htxt')
                            delete(htxt)
                        end
                        
                        jacs=matlabFunction(jac_sys);
                        cpl=matlabFunction(cpoly);
                        
                        for ii=1:nrows
                            for kk=1:ncols
                                iter=kk+ncols*(ii-1);
                                total_tim=(iter)*100/(total_steps);
                                set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                                c = clock; % [year month day hour minute seconds]
                                TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                                set(handles.timestamp,'String',['Time: ',TSMP])
                                
                                %Update sparsity matrix
                                jac_P=jacs(S1new(ii,kk),X1new(ii,kk),S2new(ii,kk),X2new(ii,kk));
                                colormap(hot)
                                xtxt=[param1,' = ',num2str(range1(yln(kk)))];
                                ytxt=[param2,' = ',num2str(range2(xln(kk)))];
                                fpx=['Fixed Point ',num2str(ii)];
                                
                                if isnan(jac_P)==0
                                    minm=min(min(double(jac_P)));
                                    maxm=max(max(double(jac_P)));
                                    imagesc(double(jac_P),[minm maxm])
                                    text(-2.5,2.5,strvcat('Jacobian Matrix',xtxt,ytxt,fpx),'interpreter','latex','horiz','left','vert','middle','fontsize',14);
                                else
                                    imagesc(ones(4,4),[0 1])
                                    text(1.5,2.5,'No solution')
                                    text(-2.5,2.5,strvcat('Jacobian Matrix',xtxt,ytxt,fpx),'interpreter','latex','horiz','left','vert','middle','fontsize',14);
                                end
                                set(gca,'Xticklabel',cellstr(var_order),'Yticklabel',cellstr(var_order))
                                pause(0.000000000000000001)
                                %Characteristic polynomial
                                
                                if S1new(ii,kk)>=0 && X1new(ii,kk) >= 0 && S2new(ii,kk)>=0 && X2new(ii,kk) >= 0
                                    cpoly=cpl(S1new(ii,kk),X1new(ii,kk),S2new(ii,kk),X2new(ii,kk));
                                    
                                    eigen=roots(cpoly);
                                    
                                    if all(real(eigen<0))==1 && X1new(ii,kk)==0 && X2new(ii,kk)==0
                                        stz(ii,kk)=1; sta(ii,kk)=1;
                                    elseif all(real(eigen<0))==1 && X1new(ii,kk)>=0 && X2new(ii,kk)==0
                                        stz(ii,kk)=2; sta(ii,kk)=1;
                                    elseif all(real(eigen<0))==1 && X1new(ii,kk)>=0 && X2new(ii,kk)>=0
                                        stz(ii,kk)=3; sta(ii,kk)=1;
                                    elseif all(real(eigen<0))==0 && X1new(ii,kk)==0 && X2new(ii,kk)==0
                                        stz(ii,kk)=1; sta(ii,kk)=0.1;
                                    elseif all(real(eigen<0))==0 && X1new(ii,kk)>=0 && X2new(ii,kk)==0
                                        stz(ii,kk)=2; sta(ii,kk)=0.2;
                                    elseif all(real(eigen<0))==0 && X1new(ii,kk)>=0 && X2new(ii,kk)>=0
                                        stz(ii,kk)=3; sta(ii,kk)=0.3;
                                    end
                                else
                                    stz(ii,kk)=0; sta(ii,kk)=0;
                                end
                            end
                        end
                        
                        %Determine which points in SS are stable
                        [ln,wd]=size(fst{1});
                        [nr,nc]=size(stz);
                        
                        for km=1:nr
                            snr{km}=reshape(stz(km,:),ln,wd);
                            sna{km}=reshape(sta(km,:),ln,wd);
                        end
                        
                        z1=1:ln;
                        z2=1:wd;
                        
                        for km=1:nr
                            b=1;
                            combined_s{km}=snr{km};
                            combined_a{km}=sna{km};
                            for kk=z1
                                for kj=z2
                                    sz(km,b)=combined_s{km}(kk,kj);
                                    sza(km,b)=combined_a{km}(kk,kj);
                                    b=b+1;
                                end
                                
                            end
                            
                        end
                        
                        for kk=1:length(sz)
                            
                            comb_s=[sz(:,kk),sza(:,kk)];
                            
                            und=unique(comb_s,'rows');
                            
                            unind=find(und(:,1)>0);
                            
                            ssund=und(unind,:);
                            
                            if all(ssund(:,1)==1)
                                Js(kk)=1;
                            elseif all(ssund(:,1)==2)
                                Js(kk)=2;
                            elseif all(ssund(:,1)==3)
                                Js(kk)=3;
                            elseif any(ssund(:,1)==1) && any(ssund(:,1)==2) && ~ismember(3,ssund(:,1))
                                Js(kk)=4;
                            elseif any(ssund(:,1)==1) && any(ssund(:,1)==3) && ~ismember(2,ssund(:,1))
                                Js(kk)=5;
                            elseif any(ssund(:,1)==1) && any(ssund(:,1)==2) && any(ssund(:,1)==3)
                                Js(kk)=6;
                            elseif any(ssund(:,1)==2) && any(ssund(:,1)==3) && ~ismember(1,ssund(:,1))
                                Js(kk)=7;
                            end
                            
                            if all(ssund(:,2)==1)
                                Jsf(kk)=Js(kk);
                            elseif any(ssund(:,2)==0.1)
                                funstab=find(ssund(:,2)==0.1);
                                unss=ssund(funstab,1);
                                if length(unss)==1 || length(unss)==2
                                    Jsf(kk)=Js(kk)+0.51;
                                elseif length(unss)==3
                                    Jsf(kk)=0;
                                end
                            elseif any(ssund(:,2)==0.2)
                                funstab=find(ssund(:,2)==0.2);
                                unss=ssund(funstab,1);
                                if length(unss)==1 || length(unss)==2
                                    Jsf(kk)=Js(kk)+0.52;
                                elseif length(unss)==3
                                    Jsf(kk)=0;
                                end
                                
                            elseif any(ssund(:,2)==0.3)
                                funstab=find(ssund(:,2)==0.3);
                                unss=ssund(funstab,1);
                                if length(unss)==1 || length(unss)==2
                                    Jsf(kk)=Js(kk)+0.53;
                                elseif length(unss)==3
                                    Jsf(kk)=0;
                                end
                                
                            end
                            
                        end
                        
                        Jj=reshape(Jsf,ln,wd);
                        
                        final_stable=Jj';
                        
                        %Store values for plotting
                        handles.stzp=final_stable; handles.ssts=floor(Js');
                        handles.fe=final_exist;
                        handles.range1=range1;
                        handles.range2=range2; handles.siz1=s1a; handles.siz2=s2a; handles.prm1=param1; handles.prm2=param2;
                        bifstr=strvcat('Select fixed-point to view...','Steady-states','Stability');
                        set(handles.bifursp,'Visible','on','String',bifstr)
                        
                    case 'Mupad'
                        
                        %Get Jacobian of system
                        sym_vars = symvar(eqs);
                        var_new=[]; var_new_char=[];
                        for k=1:length(sym_vars)
                            if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'S2') || strcmp(char(sym_vars(k)),'X2') || strcmp(char(sym_vars(k)),strtrim(param1)) || strcmp(char(sym_vars(k)),strtrim(param2))
                                var_new_char=var_new_char;
                            else
                                var_new=[var_new;sym_vars(k)];
                                var_new_char=strvcat(var_new_char,char(sym_vars(k)));
                            end
                        end
                        for kk=1:length(var_new_char)
                            eqs=subs(eqs,char(var_new(kk)),subs(var_new_char(kk,:)));
                        end
                        
                        var_order=strvcat('S1','X1','S2','X2');
                        jac_sys=jacobian(vpa(eqs),sym({'S1','X1','S2','X2'}));
                        
                        xln=repmat(1:s1a,1,s2a);
                        ylm=repmat(1:s1a,s2a,1);
                        yln=reshape(ylm,1,s1a*s2a);
                        
                        jacs=matlabFunction(jac_sys);
                        set(handles.func_prog,'String','Running: Multiple operating points - Mupad','ForegroundColor','r')
                        sol=feval(symengine,'solve',[eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0],[S1,X1,S2,X2]);
                        
                        %for the progess bar
                        total_steps=s1a*s2a;
                        no_steps=zeros(s1a,s2a);
                        cla(handles.spmatrix)
                        axes(handles.spmatrix)
                        cc=0; dd=0;
                        
                        syms D S1in %Check this
                        [X,Y]=meshgrid(range2,range1);
                        tmp = subs(sol, {S1in D}, {X, Y});
                        % then further process the "tmp" array
                        %print tmp
                         
                        for ii=vpa(range2)
                            dd=dd+1; cc=0;
                            for kk=vpa(range1)
                                cc=cc+1;
                                iter=cc+s2a*(dd-1);
                                total_tim=(iter)*100/(total_steps);
                                set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                                c = clock; % [year month day hour minute seconds]
                                TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                                set(handles.timestamp,'String',['Time: ',TSMP])
                                pause(0.000000000000001)
                                
                                if ii>0 && kk>0
                                    temp=vpa(subs(sol,[sym(param1) sym(param2)],[kk ii]));
                                    members=children(temp);
                                    solutions=members(2);
                                else
                                    solutions=zeros(1,4);
                                end
                                
                                nS=zeros(4,numel(solutions));
                                for i=1:numel(solutions)
                                    if ii>0 && kk>0
                                        nS(:,i)=double(solutions(i));
                                    else
                                        nS(:,i)=solutions(i);
                                    end
                                    jac_P=jacs(range2(xln(dd)),nS(1,i),nS(3,i),nS(2,i),nS(4,i));
                                    
                                    colormap(hot) %Option?
                                    xtxt=[param1,' = ',num2str(range1(yln(cc)))];
                                    ytxt=[param2,' = ',num2str(range2(xln(dd)))];
                                    fpx=['Fixed Point ',num2str(i)];
                                    
                                    % FIX THIS                        if isnan(jac_P)==0
                                    %                                     minm=min(min(double(jac_P)));
                                    %                                     maxm=max(max(double(jac_P)));
                                    %                                     imagesc(double(jac_P),[minm maxm])
                                    %                                     text(-3.9,3.2,strvcat('Jacobian Matrix',xtxt,ytxt,fpx),'interpreter','latex','horiz','left','vert','middle','fontsize',14);
                                    %                                 else
                                    %                                     imagesc(ones(4,4),[0 1])
                                    %                                     text(1.5,2.5,'No solution')
                                    %                                     text(-3.9,3.2,strvcat('Jacobian Matrix',xtxt,ytxt,fpx),'interpreter','latex','horiz','left','vert','middle','fontsize',14);
                                    %                                 end
                                    set(gca,'Xticklabel',cellstr(var_order),'Yticklabel',cellstr(var_order),'Xtick',1:6,'Ytick',1:6)
                                    %Characteristic polynomial
                                    
                                    cpoly=charpoly(jac_P);
                                    
                                    eigen=roots(cpoly);
                                    
                                    if all(real(eigen<0))==1
                                        stab_k=1;
                                    else
                                        stab_k=0;
                                    end
                                    %Get existence and stability at point
                                    fpn=[nS(2,i),nS(4,i)];
                                    jk=0;
                                    if nS(1,i)>=0 && nS(2,i) >= 0 && nS(3,i)>=0 && nS(4,i) >= 0
                                        
                                        if nnz(fpn)==2
                                            jk=3;
                                        elseif nnz(fpn)==0
                                            jk=1;
                                        elseif nnz(fpn)==1
                                            nzu=find(fpn>0);
                                            if nzu==1
                                                jk=2;
                                            elseif nzu==2
                                                jk=4;
                                            end
                                            
                                        end
                                    end
                                    
                                    JmS = jk.*stab_k; jz=1-stab_k;
                                    JmU = jk.*jz;
                                    %Stable uni-stable
                                    ind = cc+(dd-1)*s2a;
                                    if stab_k==1 && jk>0
                                        Js(i,ind)=jk; Ss(i,ind) = {'S'}; Jn(i,ind)=jk;
                                    elseif stab_k==0 && jk>0
                                        Js(i,ind)=jk; Ss(i,ind) = {'U'}; Jn(i,ind)=jk+8.2;
                                    else
                                        Js(i,ind)=0; Ss(i,ind) = {'0'}; Jn(i,ind)=0;
                                    end
                                end
                            end
                            
                        end

                        sempS=cellfun(@isempty,Ss);
                        Ss(sempS)={'0'};
                        
                        %Unique Js
                        [unJs,unc,und]=unique(Js','rows');
                        
                        UniqueJs=Js(:,unc);
                        
                        %Unique Ss
                        [unSs,uns,unt]=unique(cell2mat(Ss'),'rows');
                        
                        UniqueSs=Ss(:,uns);
                        [m,n]=size(UniqueSs);
                        %Find unique positions in SS
                        
                        out = sprintf('%.0f', Js) - '0';
                        outU = unique(out);
                        outF = outU(unique(outU)>0); ssF=[];
                        for kF=1:length(outF)
                            ssF=strvcat(ssF,['SS',num2str(outF(kF))]);
                        end
                        jjF=[];
                        for k=1:n
                            cc=find(unt==k);
                            dd=Ss(:,cc);
                            ee=Js(:,cc);
                            [una,unb,unc]=unique(cell2mat(dd'),'rows');
                            qF=cell.empty(0,length(outF));
                            
                            %                  if unb==1
                            q=ee(:,1);
                            q1=unique(q(q>0));
                            q2=sort(q1);
                            q_stab=[];
                            for kk=1:length(q2)
                                N=find(outF==q2(kk));
                                q3=find(q==q1(kk));
                                q_stab=dd(q3,1);
                                if length(q3)>1
                                    jstb=strjoin(q_stab);
                                    q_stab=cellstr(strrep(jstb,' ','/'));
                                end
                                
                                qF{N}=cell2mat(q_stab);
                            end
                            
                            jjF=[jjF;strcat([['$\mathcal{J}_{',num2str(k),'}$'],qF])];
                            
                            % end
                            
                        end
                        
                        %%%%%MATCH JS AND SS%%%%%%%%%%
                        
                        sJ=sum(Jn,1);
                        rJ=reshape(sJ,length(range1),length(range2));
                        
                        axes(handles.trajectoryplot)
                        contourf(range1,range2,rJ','LineStyle','none','LevelStep',1)
                        set(handles.colormp,'Visible','on');
                        %LABELS
                        xlabel(param1)
                        ylabel(param2)
                        
                        %Populate stability matrix
                        cla(handles.stabmatrix)
                        axes(handles.stabmatrix)
                        
                        [a,b]=size(jjF);
                        
                        sm=cell(a+1,length(outF)+1);
                        ssc=cellstr(ssF);
                        ssct=ssc';
                        sm(1,2:length(outF)+1)=ssct;
                        sm(2:a+1,1:b)=jjF;
                        
                        axes(handles.stabmatrix)
                        [nm,nt]=size(sm);
                        yn=ones(1,nt);
                        xn=0;
                        for j=1:nt
                            xx(j)=(xn-0.16)+0.16;
                            xn=xn+0.16;
                        end
                        for k=1:nm
                            yy=yn*(0.8-(k-1)/10);
                            handles.stbtxt=text(xx+0.1,yy,sm(k,:),'interpreter','latex','horiz','left','vert','middle','fontsize',10);
                        end
                        
                        %Change search arrays to numeric
                        SsN=char(Ss); [rS,cS]=size(Ss);
                        UsN=char(UniqueSs); [rU,cU]=size(UniqueSs);
                        
                        SsN(SsN=='S')='1'; UsN(UsN=='S')='1';
                        SsN(SsN=='U')='2'; UsN(UsN=='U')='2';
                        
                        SsNu=str2num(SsN); UsNu=str2num(UsN);
                        %Reshape and transpose
                        SsNk=reshape(SsNu,rS,cS); UsNk=reshape(UsNu,rU,cU);
                        SsNr=SsNk'; UsNr=UsNk';
                        
                        %Print locations of Jn
                        xl=1:s1a;
                        yl=1:s2a;
                        [Xl Yl]=meshgrid(xl,yl); Xl=Xl'; Yl=Yl';
                        axes(handles.trajectoryplot)
                        for kj=1:cU
                            
                            Jna=ismember(SsNr,UsNr(kj,:),'rows');
                            Jnj=reshape(Jna,s1a,s2a);
                            
                            jnj0 = mean(Xl(Jnj==1)); jnj1 = mean(Yl(Jnj==1));
                            text(jnj0*range1(end)/length(range1),jnj1*range2(end)/length(range2),jjF(kj,1),'color','b','interpreter','latex','backgroundcolor','w','fontsize',14)
                        end
                        %Turn off options
                        set(handles.bifursp,'Visible','off')
                        
                end
                set(handles.func_prog,'String','Completed: Multiple operating points','ForegroundColor',[0 0.6 1])
                guidata(hObject,handles);
                
                
            % Five equations
            elseif noeq==5
                [S1_fix_num,X1_fix_num,S2_fix_num,X2_fix_num,S3_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0],'S1','X1','S2','X2','S3');
                set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
                set(handles.s2_check,'Enable','on','value',1); set(handles.x2_check,'Enable','on','value',1);
                set(handles.s3_check,'Enable','on','value',1); set(handles.x3_check,'Enable','off','value',0);
                %Substitute operating parameters
                r1=str2double(get(handles.min_p1,'String')); r2=str2double(get(handles.max_p1,'String')); s1a=str2double(get(handles.step_p1,'String'));
                r3=str2double(get(handles.min_p2,'String')); r4=str2double(get(handles.max_p2,'String')); s2a=str2double(get(handles.step_p2,'String'));
                range1=linspace(r1,r2,s1a);
                range2=linspace(r3,r4,s2a);
                S1n = matlabFunction(S1_fix_num); X1n = matlabFunction(X1_fix_num); S2n = matlabFunction(S2_fix_num); X2n = matlabFunction(X2_fix_num); S3n = matlabFunction(S3_fix_num);
                
                lfp=length(S1_fix_num);
                cn=1;
                for k=range1
                    for kk=range2
                        sss1(1:lfp,cn)=S1n(kk,k);
                        ssx1(1:lfp,cn)=X1n(kk,k);
                        sss2(1:lfp,cn)=S2n(kk,k);
                        ssx2(1:lfp,cn)=X2n(kk,k);
                        sss3(1:lfp,cn)=S3n(kk,k);
                        cn=cn+1;
                    end
                end
                S1new=sss1;
                X1new=ssx1;
                S2new=sss2;
                X2new=ssx2;
                S3new=sss3;
                
                %need to get rid of unrealistic fixed points (i.e. those with negative values)
                [nrows,ncols]=size(sss1);
                
                %Determine steady-states for each point
                
                %Convert positives to 1, negatives to 0
                ssx1(ssx1>0)=5.5; ssx1(ssx1<0)=nan;
                ssx2(ssx2>0)=6.5; ssx2(ssx2<0)=nan;
                
                %Check for existence
                sss1(sss1<0)=nan; sss1(sss1>0)=0;
                sss2(sss2<0)=nan; sss2(sss2>0)=0;
                sss3(sss3<0)=nan; sss3(sss3>0)=0;
                
                %Create handle of matrices
                for ks=1:size(ssx1,1);
                    ssxr1{ks}=reshape(ssx1(ks,:),s1a,s2a);
                    ssxr2{ks}=reshape(ssx2(ks,:),s1a,s2a);
                    sssr1{ks}=reshape(sss1(ks,:),s1a,s2a);
                    sssr2{ks}=reshape(sss2(ks,:),s1a,s2a);
                    sssr3{ks}=reshape(sss3(ks,:),s1a,s2a);
                    ssx{ks}=ssxr1{ks}+ssxr2{ks}+sssr1{ks}+sssr2{ks}+sssr3{ks};
                    ssx{ks}(ssx{ks}==0)=1; %SS1
                    ssx{ks}(ssx{ks}==5.5)=2; %SS2
                    ssx{ks}(ssx{ks}==6.5)=4; %SS4
                    ssx{ks}(ssx{ks}==12)=3; %SS3
                    ssx{ks}(isnan(ssx{ks}))=0; %Does not exist
                end
                handles.ssxr=ssx;
                axes(handles.trajectoryplot)
                fst=handles.ssxr;
                final_exist=sum(cat(3,fst{:}),3);
                contourf(range1,range2,final_exist,'LineStyle','none','LevelStep',1)
                set(handles.colormp,'Visible','on');
                %LABELS
                xlabel(param1)
                ylabel(param2)
                
                %Progress bar
                total_steps=ncols;
                
                set(handles.func_prog,'String','Running: Creating matrix of FP solutions','ForegroundColor','r')
                
                %Get Jacobian of system
                sym_vars = symvar(eqs);
                var_new=[]; var_new_char=[];
                for k=1:length(sym_vars)
                    if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'S2') || strcmp(char(sym_vars(k)),'X2') || strcmp(char(sym_vars(k)),'S3')
                        var_new=var_new;
                        var_new_char=var_new_char;
                    else
                        var_new=[var_new;sym_vars(k)];
                        var_new_char=strvcat(var_new_char,char(sym_vars(k)));
                    end
                end
                for kk=1:length(var_new_char)
                    eqs=subs(eqs,char(var_new(kk)),subs(var_new_char(kk,:)));
                end
                
                var_order=strvcat('S1','X1','S2','X2','S3');
                jac_sys=jacobian(vpa(eqs),sym({'S1','X1','S2','X2','S3'}));
                
                cpoly=charpoly(jac_sys);
                [nrows,ncols]=size(S1new);
                
                %for the progess bar
                total_steps=numel(S1new);
                no_steps=zeros(nrows,ncols);
                
                set(handles.func_prog,'String','Running: Multiple operating points','ForegroundColor','r')
                axes(handles.spmatrix)
                
                xln=repmat(1:s1a,1,s2a);
                ylm=repmat(1:s1a,s2a,1);
                yln=reshape(ylm,1,s1a*s2a);
                
                if exist('htxt')
                    delete(htxt)
                end
                
                jacs=matlabFunction(jac_sys);
                cpl=matlabFunction(cpoly);
                
                for ii=1:nrows
                    for kk=1:ncols
                        iter=kk+ncols*(ii-1);
                        total_tim=(iter)*100/(total_steps);
                        set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                        c = clock; % [year month day hour minute seconds]
                        TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                        set(handles.timestamp,'String',['Time: ',TSMP])
                        
                        %Update sparsity matrix
                        jac_P=jacs(S1new(ii,kk),X1new(ii,kk),S2new(ii,kk),X2new(ii,kk),S3new(ii,kk));
                        colormap(hot)
                        xtxt=[param1,' = ',num2str(range1(yln(kk)))];
                        ytxt=[param2,' = ',num2str(range2(xln(kk)))];
                        fpx=['Fixed Point ',num2str(ii)];
                        
                        if isnan(jac_P)==0
                            minm=min(min(double(jac_P)));
                            maxm=max(max(double(jac_P)));
                            imagesc(double(jac_P),[minm maxm])
                            text(-2.5,2.5,strvcat('Jacobian Matrix',xtxt,ytxt,fpx),'interpreter','latex','horiz','left','vert','middle','fontsize',14);
                        else
                            imagesc(ones(4,4),[0 1])
                            text(1.5,2.5,'No solution')
                            text(-2.5,2.5,strvcat('Jacobian Matrix',xtxt,ytxt,fpx),'interpreter','latex','horiz','left','vert','middle','fontsize',14);
                        end
                        set(gca,'Xticklabel',cellstr(var_order),'Yticklabel',cellstr(var_order))
                        pause(0.000000000000000001)
                        %Characteristic polynomial
                        
                        if S1new(ii,kk)>=0 && X1new(ii,kk) >= 0 && S2new(ii,kk)>=0 && X2new(ii,kk) >= 0 && S3new(ii,kk)>=0
                            cpoly=cpl(S1new(ii,kk),X1new(ii,kk),S2new(ii,kk),X2new(ii,kk),S3new(ii,kk));
                            
                            eigen=roots(cpoly);
                            
                            if all(real(eigen<0))==1 && X1new(ii,kk)==0 && X2new(ii,kk)==0
                                stz(ii,kk)=1; sta(ii,kk)=1;
                            elseif all(real(eigen<0))==1 && X1new(ii,kk)>=0 && X2new(ii,kk)==0
                                stz(ii,kk)=2; sta(ii,kk)=1;
                            elseif all(real(eigen<0))==1 && X1new(ii,kk)>=0 && X2new(ii,kk)>=0
                                stz(ii,kk)=3; sta(ii,kk)=1;
                            elseif all(real(eigen<0))==0 && X1new(ii,kk)==0 && X2new(ii,kk)==0
                                stz(ii,kk)=1; sta(ii,kk)=0.1;
                            elseif all(real(eigen<0))==0 && X1new(ii,kk)>=0 && X2new(ii,kk)==0
                                stz(ii,kk)=2; sta(ii,kk)=0.2;
                            elseif all(real(eigen<0))==0 && X1new(ii,kk)>=0 && X2new(ii,kk)>=0
                                stz(ii,kk)=3; sta(ii,kk)=0.3;
                            end
                        else
                            stz(ii,kk)=0; sta(ii,kk)=0;
                        end
                        
                        
                    end
                end
                
                %Determine which points in SS are stable
                [ln,wd]=size(fst{1});
                [nr,nc]=size(stz);
                
                for km=1:nr
                    snr{km}=reshape(stz(km,:),ln,wd);
                    sna{km}=reshape(sta(km,:),ln,wd);
                end
                
                z1=1:ln;
                z2=1:wd;
                
                for km=1:nr
                    b=1;
                    combined_s{km}=snr{km};
                    combined_a{km}=sna{km};
                    for kk=z1
                        for kj=z2
                            sz(km,b)=combined_s{km}(kk,kj);
                            sza(km,b)=combined_a{km}(kk,kj);
                            b=b+1;
                        end
                        
                    end
                    
                end
                
                for kk=1:length(sz)
                    
                    comb_s=[sz(:,kk),sza(:,kk)];
                    
                    und=unique(comb_s,'rows');
                    
                    unind=find(und(:,1)>0);
                    
                    ssund=und(unind,:);
                    
                    if all(ssund(:,1)==1)
                        Js(kk)=1;
                    elseif all(ssund(:,1)==2)
                        Js(kk)=2;
                    elseif all(ssund(:,1)==3)
                        Js(kk)=3;
                    elseif any(ssund(:,1)==1) && any(ssund(:,1)==2) && ~ismember(3,ssund(:,1))
                        Js(kk)=4;
                    elseif any(ssund(:,1)==1) && any(ssund(:,1)==3) && ~ismember(2,ssund(:,1))
                        Js(kk)=5;
                    elseif any(ssund(:,1)==1) && any(ssund(:,1)==2) && any(ssund(:,1)==3)
                        Js(kk)=6;
                    elseif any(ssund(:,1)==2) && any(ssund(:,1)==3) && ~ismember(1,ssund(:,1))
                        Js(kk)=7;
                    end
                    
                    if all(ssund(:,2)==1)
                        Jsf(kk)=Js(kk);
                    elseif any(ssund(:,2)==0.1)
                        funstab=find(ssund(:,2)==0.1);
                        unss=ssund(funstab,1);
                        if length(unss)==1 || length(unss)==2
                            Jsf(kk)=Js(kk)+0.51;
                        elseif length(unss)==3
                            Jsf(kk)=0;
                        end
                    elseif any(ssund(:,2)==0.2)
                        funstab=find(ssund(:,2)==0.2);
                        unss=ssund(funstab,1);
                        if length(unss)==1 || length(unss)==2
                            Jsf(kk)=Js(kk)+0.52;
                        elseif length(unss)==3
                            Jsf(kk)=0;
                        end
                        
                    elseif any(ssund(:,2)==0.3)
                        funstab=find(ssund(:,2)==0.3);
                        unss=ssund(funstab,1);
                        if length(unss)==1 || length(unss)==2
                            Jsf(kk)=Js(kk)+0.53;
                        elseif length(unss)==3
                            Jsf(kk)=0;
                        end
                        
                    end
                    
                end
                
                Jj=reshape(Jsf,ln,wd);
                
                final_stable=Jj';
                
                %Store values for plotting
                handles.stzp=final_stable; handles.ssts=floor(Js');
                handles.fe=final_exist;
                handles.range1=range1;
                handles.range2=range2; handles.siz1=s1a; handles.siz2=s2a; handles.prm1=param1; handles.prm2=param2;
                bifstr=strvcat('Select fixed-point to view...','Steady-states','Stability');
                set(handles.bifursp,'Visible','on','String',bifstr)
                set(handles.func_prog,'String','Completed: Multiple operating points','ForegroundColor',[0 0.6 1])
                guidata(hObject,handles);
            % Six Equations
            elseif noeq==6
                warning off
                set(handles.progress,'String',['Progress: ',num2str(0),'%'])
                syms S1 X1 S2 X2 S3 X3 D S1in
                
                %Make assumptions
                assume(S1>=0); assume(X1>=0); assume(S2>=0); assume(X2>=0); assume(S3>=0); assume(X3>=0);
                r1=str2double(get(handles.min_p1,'String')); r2=str2double(get(handles.max_p1,'String')); s1a=str2double(get(handles.step_p1,'String'));
                r3=str2double(get(handles.min_p2,'String')); r4=str2double(get(handles.max_p2,'String')); s2a=str2double(get(handles.step_p2,'String'));
                range1=linspace(r1,r2,s1a);
                range2=linspace(r3,r4,s2a);
                %Solve equations at equilibrium
                ssol=get(handles.solmeth_select,'String');
                vsol=get(handles.solmeth_select,'Value');
                calculationType=char(ssol(vsol,:));
                switch calculationType
                    case 'Matlab'
                        [S1_fix_num,X1_fix_num,S2_fix_num,X2_fix_num,S3_fix_num,X3_fix_num]=solve([eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0 eq6_numerical==0],[S1,X1,S2,X2,S3,X3]);
                        
                        set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
                        set(handles.s2_check,'Enable','on','value',1); set(handles.x2_check,'Enable','on','value',1);
                        set(handles.s3_check,'Enable','on','value',1); set(handles.x3_check,'Enable','on','value',1);
                        %Substitute operating parameters
                        
                        S1n = matlabFunction(S1_fix_num); X1n = matlabFunction(X1_fix_num); S2n = matlabFunction(S2_fix_num); X2n = matlabFunction(X2_fix_num); S3n = matlabFunction(S3_fix_num); X3n = matlabFunction(X3_fix_num);
                        %eq1n=matlabFunction(eq1_numerical); eq2n=matlabFunction(eq2_numerical); eq3n=matlabFunction(eq3_numerical); eq4n=matlabFunction(eq4_numerical); eq5n=matlabFunction(eq5_numerical); eq6n=matlabFunction(eq6_numerical);
                        
                        lfp=length(S1_fix_num);
                        cn=1;
                        for k=range1
                            for kk=range2
                                sss1(1:lfp,cn)=S1n(kk,k);
                                ssx1(1:lfp,cn)=X1n(kk,k);
                                sss2(1:lfp,cn)=S2n(kk,k);
                                ssx2(1:lfp,cn)=X2n(kk,k);
                                nai_s3=nargin(S3n);
                                switch nai_s3
                                    case 1
                                        sss3(1:lfp,cn)=S3n(kk);
                                    case 2
                                        sss3(1:lfp,cn)=S3n(kk,k);
                                end
                                ssx3(1:lfp,cn)=X3n(kk,k);
                                cn=cn+1;
                            end
                        end
                        
                        S1new=sss1;
                        X1new=ssx1;
                        S2new=sss2;
                        X2new=ssx2;
                        S3new=sss3;
                        X3new=ssx3;
                        
                        %need to get rid of unrealistic fixed points (i.e. those with negative values)
                        [nrows,ncols]=size(sss1);
                        
                        %Get Jacobian of system
                        sym_vars = symvar(eqs);
                        var_new=[]; var_new_char=[];
                        for k=1:length(sym_vars)
                            if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'S2') || strcmp(char(sym_vars(k)),'X2') || strcmp(char(sym_vars(k)),'S3') || strcmp(char(sym_vars(k)),'X3') || strcmp(char(sym_vars(k)),strtrim(param1)) || strcmp(char(sym_vars(k)),strtrim(param2))
                                var_new=var_new;
                                var_new_char=var_new_char;
                            else
                                var_new=[var_new;sym_vars(k)];
                                var_new_char=strvcat(var_new_char,char(sym_vars(k)));
                            end
                        end
                        for kk=1:length(var_new_char)
                            eqs=subs(eqs,char(var_new(kk)),subs(var_new_char(kk,:)));
                        end
                        
                        var_order=strvcat('S1','X1','S2','X2','S3','X3');
                        jac_sys=jacobian(vpa(eqs),sym({'S1','X1','S2','X2','S3','X3'}));
                        
                        %for the progess bar
                        total_steps=numel(S1new);
                        no_steps=zeros(nrows,ncols);
                        cla(handles.spmatrix)
                        set(handles.func_prog,'String','Running: Multiple operating points','ForegroundColor','r')
                        axes(handles.spmatrix)
                        
                        xln=repmat(1:s1a,1,s2a);
                        ylm=repmat(1:s1a,s2a,1);
                        yln=reshape(ylm,1,s1a*s2a);
                        
                        jacs=matlabFunction(jac_sys);
                        
                        for ii=1:nrows
                            for kk=1:ncols
                                iter=kk+ncols*(ii-1);
                                total_tim=(iter)*100/(total_steps);
                                set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                                c = clock; % [year month day hour minute seconds]
                                TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                                set(handles.timestamp,'String',['Time: ',TSMP])
                                
                                %Update sparsity matrix
                                jac_P=jacs(range2(xln(kk)),S1new(ii,kk),X1new(ii,kk),S2new(ii,kk),X2new(ii,kk),S3new(ii,kk),X3new(ii,kk));
                                colormap(hot)
                                xtxt=[param1,' = ',num2str(range1(yln(kk)))];
                                ytxt=[param2,' = ',num2str(range2(xln(kk)))];
                                fpx=['Fixed Point ',num2str(ii)];
                                
                                if isnan(jac_P)==0
                                    minm=min(min(double(jac_P)));
                                    maxm=max(max(double(jac_P)));
                                    imagesc(double(jac_P),[minm maxm])
                                    text(-3.9,3.2,strvcat('Jacobian Matrix',xtxt,ytxt,fpx),'interpreter','latex','horiz','left','vert','middle','fontsize',14);
                                else
                                    imagesc(ones(4,4),[0 1])
                                    text(1.5,2.5,'No solution')
                                    text(-3.9,3.2,strvcat('Jacobian Matrix',xtxt,ytxt,fpx),'interpreter','latex','horiz','left','vert','middle','fontsize',14);
                                end
                                set(gca,'Xticklabel',cellstr(var_order),'Yticklabel',cellstr(var_order),'Xtick',1:6,'Ytick',1:6)
                                pause(0.000000000000000001)
                                %Characteristic polynomial
                                
                                cpoly=charpoly(jac_P);
                                
                                eigen=roots(cpoly);
                                
                                if all(real(eigen<0))==1
                                    stab_k=1;
                                else
                                    stab_k=0;
                                end
                                %Get existence and stability at point
                                fpn=[X1new(ii,kk),X2new(ii,kk),X3new(ii,kk)];
                                jk=0;
                                if S1new(ii,kk)>=0 && X1new(ii,kk) >= 0 && S2new(ii,kk)>=0 && X2new(ii,kk) >= 0 && S3new(ii,kk)>=0 && X3new(ii,kk) >=0
                                    
                                    if nnz(fpn)==3
                                        jk=6;
                                    elseif nnz(fpn)==0
                                        jk=1;
                                    elseif nnz(fpn)==1
                                        nzu=find(fpn>0);
                                        if nzu==1
                                            jk=3;
                                        elseif nzu==2
                                            jk=7;
                                        elseif nzu==3
                                            jk=2;
                                        end
                                    elseif nnz(fpn)==2
                                        nzu2=find(fpn>0);
                                        if sum(nzu2)==3
                                            jk=4;
                                        elseif sum(nzu2)==4
                                            jk=5;
                                        elseif sum(nzu2)==5
                                            jk=8;
                                        end
                                    end
                                end
                                
                                
                                JmS = jk.*stab_k; jz=1-stab_k;
                                JmU = jk.*jz;
                                %Stable unistable
                                if stab_k==1 && jk>0
                                    Js(ii,kk)=jk; Ss(ii,kk) = {'S'}; Jn(ii,kk)=jk;
                                elseif stab_k==0 && jk>0
                                    Js(ii,kk)=jk; Ss(ii,kk) = {'U'}; Jn(ii,kk)=jk+8.2;
                                else
                                    Js(ii,kk)=0; Ss(ii,kk) = {'0'}; Jn(ii,kk)=0;
                                end
                            end
                        end
                        
                    case 'Mupad'
                        %Get Jacobian of system
                        sym_vars = symvar(eqs);
                        var_new=[]; var_new_char=[];
                        for k=1:length(sym_vars)
                            if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'S2') || strcmp(char(sym_vars(k)),'X2') || strcmp(char(sym_vars(k)),'S3') || strcmp(char(sym_vars(k)),'X3') || strcmp(char(sym_vars(k)),strtrim(param1)) || strcmp(char(sym_vars(k)),strtrim(param2))
                                var_new=var_new;
                                var_new_char=var_new_char;
                            else
                                var_new=[var_new;sym_vars(k)];
                                var_new_char=strvcat(var_new_char,char(sym_vars(k)));
                            end
                        end
                        for kk=1:length(var_new_char)
                            eqs=subs(eqs,char(var_new(kk)),subs(var_new_char(kk,:)));
                        end
                        
                        var_order=strvcat('S1','X1','S2','X2','S3','X3');
                        jac_sys=jacobian(vpa(eqs),sym({'S1','X1','S2','X2','S3','X3'}));
                        
                        xln=repmat(1:s1a,1,s2a);
                        ylm=repmat(1:s1a,s2a,1);
                        yln=reshape(ylm,1,s1a*s2a);
                        
                        jacs=matlabFunction(jac_sys);
                        set(handles.func_prog,'String','Running: Multiple op. points - Mupad solver','ForegroundColor','r')

                        sol=feval(symengine,'solve',[eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0 eq6_numerical==0],[S1,X1,S2,X2,S3,X3]);

                        syms D S1in

                        set(handles.func_prog,'String','Running: Multiple op. points - Calculating stability','ForegroundColor','r')

                        drange=range2;
                        srange=range1;
                        [X,Y]=meshgrid(srange,drange);
                        result_o = cell(s1a,s2a);tic
                        parfor k=1:s2a
                            temp = vpa( subs(sol,{S1in,D},{X(k,:),Y(k,:)}) );
                            members = children(temp); %Find children of result
                            solutions = cellfun(@(c)c(2),members,'UniformOutput',false); %Find numerical solutions
                            result_o(k,:) = cellfun(@symToDouble,solutions,'UniformOutput',false);
                        end
                        toc
                        rsl=cellfun(@cell2mat,result_o,'UniformOutput',false);
                        [rslw,rsll]=size(rsl);
                        total_steps=s1a*s2a;
                        no_steps=zeros(s1a,s2a);
                        
                        tic
                        for kk=1:rslw
                            for ii=1:rsll
                                ind=ii+s2a*(kk-1);
                                total_tim=(ind)*100/(total_steps);
                                set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                                c = clock; % [year month day hour minute seconds]
                                TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                                set(handles.timestamp,'String',['Time: ',TSMP])
                                pause(0.000000000000001)
                                Zr=rsl{kk,ii};
                                nS=zeros(6,numel(Zr)/6);
                                for i=1:numel(Zr)/6
                                    nS(:,i)=Zr(:,i);
                                    jac_P=jacs(range2(ylm(ii,kk)),nS(1,i),nS(3,i),nS(5,i),nS(2,i),nS(4,i),nS(6,i));
                                    cpoly=charpoly(jac_P);
                                    
                                    eigen=roots(cpoly);
                                    
                                    if all(real(eigen<0))==1
                                        stab_k=1;
                                    else
                                        stab_k=0;
                                    end
                                    %Get existence and stability at point
                                    fpn=[nS(2,i),nS(4,i),nS(6,i)];
                                    jk=0;
                                    if nS(1,i)>=0 && nS(2,i) >= 0 && nS(3,i)>=0 && nS(4,i) >= 0 && nS(5,i)>=0 && nS(6,i) >=0
                                        
                                        if nnz(fpn)==3
                                            jk=6;
                                        elseif nnz(fpn)==0
                                            jk=1;
                                        elseif nnz(fpn)==1
                                            nzu=find(fpn>0);
                                            if nzu==1
                                                jk=3;
                                            elseif nzu==2
                                                jk=7;
                                            elseif nzu==3
                                                jk=2;
                                            end
                                        elseif nnz(fpn)==2
                                            nzu2=find(fpn>0);
                                            if sum(nzu2)==3
                                                jk=4;
                                            elseif sum(nzu2)==4
                                                jk=5;
                                            elseif sum(nzu2)==5
                                                jk=8;
                                            end
                                        end
                                    end
                                    
                                    JmS = jk.*stab_k; jz=1-stab_k;
                                    JmU = jk.*jz;
                                    
                                    if stab_k==1 && jk>0
                                        Js(i,ind)=jk; Ss(i,ind) = {'S'}; Jn(i,ind)=jk;
                                    elseif stab_k==0 && jk>0
                                        Js(i,ind)=jk; Ss(i,ind) = {'U'}; Jn(i,ind)=jk+8.2;
                                    else
                                        Js(i,ind)=0; Ss(i,ind) = {'0'}; Jn(i,ind)=0;
                                    end
                                end
                            end
                        end
                        toc
                        
                        sempS=cellfun(@isempty,Ss);
                        Ss(sempS)={'0'};
                        
                        %Unique Js
                        [unJs,unc,und]=unique(Js','rows');
                        
                        UniqueJs=Js(:,unc);
                        
                        %Unique Ss
                        [unSs,uns,unt]=unique(cell2mat(Ss'),'rows');
                        
                        UniqueSs=Ss(:,uns);
                        [m,n]=size(UniqueSs);
                        %Find unique positions in SS
                        
                        out = sprintf('%.0f', Js) - '0';
                        outU = unique(out);
                        outF = outU(unique(outU)>0); ssF=[];
                        for kF=1:length(outF)
                            ssF=strvcat(ssF,['SS',num2str(outF(kF))]);
                        end
                        jjF=[];
                        for k=1:n
                            cc=find(unt==k);
                            dd=Ss(:,cc);
                            ee=Js(:,cc);
                            [una,unb,unc]=unique(cell2mat(dd'),'rows');
                            qF=cell.empty(0,length(outF));
                            
                            q=ee(:,1);
                            q1=unique(q(q>0));
                            q2=sort(q1);
                            q_stab=[];
                            for kk=1:length(q2)
                                N=find(outF==q2(kk));
                                q3=find(q==q1(kk));
                                q_stab=dd(q3,1);
                                if length(q3)>1
                                    jstb=strjoin(q_stab);
                                    q_stab=cellstr(strrep(jstb,' ','/'));
                                end
                                
                                qF{N}=cell2mat(q_stab);
                            end
                            
                            jjF=[jjF;strcat([['$\mathcal{J}_{',num2str(k),'}$'],qF])];
                            
                        end
                        
                        %%%%%MATCH JS AND SS%%%%%%%%%%
                        
                        sJ=sum(Jn,1);
                        rJ=reshape(sJ,length(range1),length(range2));
                        switch calculationType
                            case 'Matlab'
                                rJ=rJ';
                        end
                        
                        axes(handles.trajectoryplot)
                        contourf(range1,range2,rJ','LineStyle','none','LevelStep',1)
                        set(handles.colormp,'Visible','on');
                        %LABELS
                        xlabel(param1)
                        ylabel(param2)
                        
                        %Populate stability matrix
                        cla(handles.stabmatrix)
                        axes(handles.stabmatrix)
                        
                        [a,b]=size(jjF);
                        
                        sm=cell(a+1,length(outF)+1);
                        ssc=cellstr(ssF);
                        ssct=ssc';
                        sm(1,2:length(outF)+1)=ssct;
                        sm(2:a+1,1:b)=jjF;
                        [smr,smc]=size(sm);
                        axes(handles.stabmatrix)
                        [nm,nt]=size(sm);
                        yn=ones(1,nt);
                        xn=0;
                        for j=1:nt
                            xx(j)=(xn-0.16)+0.16;
                            xn=xn+0.16;
                        end
                        for k=1:nm
                            yy=yn*(0.8-(k-1)/10);
                            if smr<9
                                handles.stbtxt=text(xx+0.1,yy,sm(k,:),'interpreter','latex','horiz','left','vert','middle','fontsize',10);
                            else
                                handles.stbtxt=text(xx+0.1,yy,sm(k,:),'interpreter','latex','horiz','left','vert','middle','fontsize',6);
                            end
                        end
                        
                        %Change search arrays to numeric
                        SsN=char(Ss); [rS,cS]=size(Ss);
                        UsN=char(UniqueSs); [rU,cU]=size(UniqueSs);
                        
                        SsN(SsN=='S')='1'; UsN(UsN=='S')='1';
                        SsN(SsN=='U')='2'; UsN(UsN=='U')='2';
                        
                        SsNu=str2num(SsN); UsNu=str2num(UsN);
                        %Reshape and transpose
                        SsNk=reshape(SsNu,rS,cS); UsNk=reshape(UsNu,rU,cU);
                        SsNr=SsNk'; UsNr=UsNk';
                        
                        %Print locations of Jn
                        xl=1:s1a;
                        yl=1:s2a;
                        [Xl Yl]=meshgrid(xl,yl); Xl=Xl'; Yl=Yl';
                        axes(handles.trajectoryplot)
                        for kj=1:cU
                            
                            Jna=ismember(SsNr,UsNr(kj,:),'rows');
                            Jnj=reshape(Jna,s1a,s2a);
                            
                            jnj0 = mean(Xl(Jnj==1)); jnj1 = mean(Yl(Jnj==1));
                            text(jnj0*range1(end)/length(range1),jnj1*range2(end)/length(range2),jjF(kj,1),'color','b','interpreter','latex','backgroundcolor','w','fontsize',14)
                        end
                        %Turn off options
                        set(handles.bifursp,'Visible','off')
                        
                        set(handles.func_prog,'String','Completed: Multiple operating points','ForegroundColor',[0 0.6 1])
                        guidata(hObject,handles);
                        
                end
            end
            
            
        elseif handles.mpltopt==1
            axes(handles.trajectoryplot)
            %cla(handles.trajectoryplot)
            bfv=get(handles.bifursp,'Value');
            r1=handles.range1; r2=handles.range2;
            if bfv==2
                
                contourf(handles.range1,handles.range2,handles.fe,'LineStyle','none','LevelStep',1)
                set(handles.colormp,'Visible','on');
                sstn=handles.ssts;
                sm=reshape(sstn,length(r1),length(r2));
                sstn_matrix=sm';
                %Determine steady-states using previously created SS matrix
                xl=1:size(sstn_matrix,1);
                yl=1:length(sstn_matrix);
                [Xl Yl]=meshgrid(xl,yl);
                cx1 = mean(Xl(sstn_matrix==1));
                cx2 = mean(Xl(sstn_matrix==2));
                cx3 = mean(Xl(sstn_matrix==3));
                cx4 = mean(Xl(sstn_matrix==4));
                cx5 = mean(Xl(sstn_matrix==5));
                cx6 = mean(Xl(sstn_matrix==6));
                cx7 = mean(Xl(sstn_matrix==7));
                cy1 = mean(Yl(sstn_matrix==1));
                cy2 = mean(Yl(sstn_matrix==2));
                cy3 = mean(Yl(sstn_matrix==3));
                cy4 = mean(Yl(sstn_matrix==4));
                cy5 = mean(Yl(sstn_matrix==5));
                cy6 = mean(Yl(sstn_matrix==6));
                cy7 = mean(Yl(sstn_matrix==7));
                
                hold on
                text(cx1*r1(end)/length(r1),cy1*r2(end)/length(r2),'SS1','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx2*r1(end)/length(r1),cy2*r2(end)/length(r2),'SS2','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx3*r1(end)/length(r1),cy3*r2(end)/length(r2),'SS3','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx4*r1(end)/length(r1),cy4*r2(end)/length(r2),'SS4','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx5*r1(end)/length(r1),cy5*r2(end)/length(r2),'SS5','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx6*r1(end)/length(r1),cy6*r2(end)/length(r2),'SS6','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx7*r1(end)/length(r1),cy7*r2(end)/length(r2),'SS7','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                
                title('Steady-states Plot')
            elseif bfv==3
                cla(handles.stabmatrix)
                contourf(handles.range1,handles.range2,handles.stzp,'LineStyle','none','LevelStep',1)
                set(handles.colormp,'Visible','on');
                szp=handles.stzp;
                %Add steady-state fixed point labels
                
                xl=1:size(szp,1);
                yl=1:length(szp);
                [Xl Yl]=meshgrid(xl,yl);
                cx0 = mean(Xl(szp==0)); cy0 = mean(Yl(szp==0));
                cx1 = mean(Xl(szp==1)); cx151 = mean(Xl(szp==1.51)); cx152 = mean(Xl(szp==1.52)); cx153 = mean(Xl(szp==1.53));
                cx2 = mean(Xl(szp==2)); cx251 = mean(Xl(szp==2.51)); cx252 = mean(Xl(szp==2.52)); cx253 = mean(Xl(szp==2.53));
                cx3 = mean(Xl(szp==3)); cx351 = mean(Xl(szp==3.51)); cx352 = mean(Xl(szp==3.52)); cx353 = mean(Xl(szp==3.53));
                cx4 = mean(Xl(szp==4)); cx451 = mean(Xl(szp==4.51)); cx452 = mean(Xl(szp==4.52)); cx453 = mean(Xl(szp==4.53));
                cx5 = mean(Xl(szp==5)); cx551 = mean(Xl(szp==5.51)); cx552 = mean(Xl(szp==5.52)); cx553 = mean(Xl(szp==5.53));
                cx6 = mean(Xl(szp==6)); cx651 = mean(Xl(szp==6.51)); cx652 = mean(Xl(szp==6.52)); cx653 = mean(Xl(szp==6.53));
                cx7 = mean(Xl(szp==7)); cx751 = mean(Xl(szp==7.51)); cx752 = mean(Xl(szp==7.52)); cx753 = mean(Xl(szp==7.53));
                cy1 = mean(Yl(szp==1)); cy151 = mean(Yl(szp==1.51)); cy152 = mean(Yl(szp==1.52)); cy153 = mean(Yl(szp==1.53));
                cy2 = mean(Yl(szp==2)); cy251 = mean(Yl(szp==2.51)); cy252 = mean(Yl(szp==2.52)); cy253 = mean(Yl(szp==2.53));
                cy3 = mean(Yl(szp==3)); cy351 = mean(Yl(szp==3.51)); cy352 = mean(Yl(szp==3.52)); cy353 = mean(Yl(szp==3.53));
                cy4 = mean(Yl(szp==4)); cy451 = mean(Yl(szp==4.51)); cy452 = mean(Yl(szp==4.52)); cy453 = mean(Yl(szp==4.53));
                cy5 = mean(Yl(szp==5)); cy551 = mean(Yl(szp==5.51)); cy552 = mean(Yl(szp==5.52)); cy553 = mean(Yl(szp==5.53));
                cy6 = mean(Yl(szp==6)); cy651 = mean(Yl(szp==6.51)); cy652 = mean(Yl(szp==6.52)); cy653 = mean(Yl(szp==6.53));
                cy7 = mean(Yl(szp==7)); cy751 = mean(Yl(szp==7.51)); cy752 = mean(Yl(szp==7.52)); cy753 = mean(Yl(szp==7.53));
                
                hold on
                text(cx0*r1(end)/length(r1),cy0*r2(end)/length(r2),'$\not\exists$','color','w','interpreter','latex','backgroundcolor','k','fontsize',14)
                text(cx1*r1(end)/length(r1),cy1*r2(end)/length(r2),'$\mathcal{J}_1$','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx2*r1(end)/length(r1),cy2*r2(end)/length(r2),'$\mathcal{J}_2$','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx3*r1(end)/length(r1),cy3*r2(end)/length(r2),'$\mathcal{J}_3$','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx4*r1(end)/length(r1),cy4*r2(end)/length(r2),'$\mathcal{J}_4$','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx5*r1(end)/length(r1),cy5*r2(end)/length(r2),'$\mathcal{J}_5$','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx6*r1(end)/length(r1),cy6*r2(end)/length(r2),'$\mathcal{J}_6$','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx7*r1(end)/length(r1),cy7*r2(end)/length(r2),'$\mathcal{J}_7$','color','r','interpreter','latex','backgroundcolor','w','fontsize',14)
                text(cx151*r1(end)/length(r1),cy151*r2(end)/length(r2),'$\mathcal{J}_8$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx251*r1(end)/length(r1),cy251*r2(end)/length(r2),'$\mathcal{J}_9$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx351*r1(end)/length(r1),cy351*r2(end)/length(r2),'$\mathcal{J}_{10}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx451*r1(end)/length(r1),cy451*r2(end)/length(r2),'$\mathcal{J}_{11}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx551*r1(end)/length(r1),cy551*r2(end)/length(r2),'$\mathcal{J}_{12}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx651*r1(end)/length(r1),cy651*r2(end)/length(r2),'$\mathcal{J}_{13}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx751*r1(end)/length(r1),cy751*r2(end)/length(r2),'$\mathcal{J}_{14}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx152*r1(end)/length(r1),cy152*r2(end)/length(r2),'$\mathcal{J}_{15}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx252*r1(end)/length(r1),cy252*r2(end)/length(r2),'$\mathcal{J}_{16}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx352*r1(end)/length(r1),cy352*r2(end)/length(r2),'$\mathcal{J}_{17}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx452*r1(end)/length(r1),cy452*r2(end)/length(r2),'$\mathcal{J}_{18}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx552*r1(end)/length(r1),cy552*r2(end)/length(r2),'$\mathcal{J}_{19}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx652*r1(end)/length(r1),cy652*r2(end)/length(r2),'$\mathcal{J}_{20}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx752*r1(end)/length(r1),cy752*r2(end)/length(r2),'$\mathcal{J}_{21}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx153*r1(end)/length(r1),cy153*r2(end)/length(r2),'$\mathcal{J}_{22}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx253*r1(end)/length(r1),cy253*r2(end)/length(r2),'$\mathcal{J}_{23}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx353*r1(end)/length(r1),cy353*r2(end)/length(r2),'$\mathcal{J}_{24}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx453*r1(end)/length(r1),cy453*r2(end)/length(r2),'$\mathcal{J}_{25}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx553*r1(end)/length(r1),cy553*r2(end)/length(r2),'$\mathcal{J}_{26}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx653*r1(end)/length(r1),cy653*r2(end)/length(r2),'$\mathcal{J}_{27}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                text(cx753*r1(end)/length(r1),cy753*r2(end)/length(r2),'$\mathcal{J}_{28}$','color','y','interpreter','latex','backgroundcolor','b','fontsize',14)
                
                lstc=[cx1,cx2,cx3,cx4,cx5,cx6,cx7,cx151,cx251,cx351,cx451,cx551,cx651,cx751,cx152,cx252,cx352,cx452,cx552,cx652,cx752,cx153,cx253,cx353,cx453,cx553,cx653,cx753];
                fnan=find(~isnan(lstc));
                
                sstxt=[1,2,3,4,5,6,7];
                
                %Stability matrix
                rn=[]; ss=[]; stab=[];
                for k=1:length(fnan)
                    rn=strvcat(rn,['$\mathcal{J}_{',num2str(fnan(k)),'}$']);
                    if fnan(k)<8
                        ss=strvcat(ss,['SS',num2str(fnan(k))]);
                    elseif 8<fnan(k) && fnan(k)<14
                        ss=strvcat(ss,['SS',num2str(fnan(k)-7)]);
                    elseif 14<fnan(k) && fnan(k)<21
                        ss=strvcat(ss,['SS',num2str(fnan(k)-14)]);
                    elseif 21<fnan(k)
                        ss=strvcat(ss,['SS',num2str(fnan(k)-21)]);
                    end
                    
                end
                
                [a,b]=size(rn);
                %Get unique steady-states
                uss=unique(ss,'rows');
                
                sm=cell(a+1,length(uss)+1);
                ssc=cellstr(uss);
                ssct=ssc';
                sm(1,2:length(uss)+1)=ssct;
                sm(2:a+1,1)=cellstr(rn);
                
                for k=1:a
                    noJ=fnan(k);
                    switch noJ
                        case 1
                            sm(find(fnan==noJ)+1,2)={'S'};
                        case 2
                            sm(find(fnan==noJ)+1,3)={'S'};
                        case 3
                            sm(find(fnan==noJ)+1,4)={'S'};
                        case 4
                            sm(find(fnan==noJ)+1,(2:3))={'S'};
                        case 5
                            sm(find(fnan==noJ)+1,[2,4])={'S'};
                        case 6
                            sm(find(fnan==noJ)+1,(2:4))={'S'};
                        case 7
                            sm(find(fnan==noJ)+1,(3:4))={'S'};
                        case 8
                            sm(find(fnan==noJ)+1,2)={'U'};
                        case 9
                            sm(find(fnan==noJ)+1,3)={'U'}; %Cannot exist
                        case 10
                            sm(find(fnan==noJ)+1,4)={'U'}; %Cannot exist
                        case 11
                            sm(find(fnan==noJ)+1,[2,3])=[{'U'},{'S'}];
                        case 12
                            sm(find(fnan==noJ)+1,[2,4])=[{'U'},{'S'}];
                        case 13
                            sm(find(fnan==noJ)+1,(2:4))=[{'U'},{'S'},{'S'}];
                        case 14
                            sm(find(fnan==noJ)+1,(3:4))={'U'}; %Cannot exist
                            
                        case 15
                            sm(find(fnan==noJ)+1,2)={'U'}; %Cannot exist
                        case 16
                            sm(find(fnan==noJ)+1,3)={'U'};
                        case 17
                            sm(find(fnan==noJ)+1,4)={'U'}; %Cannot exist
                        case 18
                            sm(find(fnan==noJ)+1,[2,3])=[{'S'},{'U'}];
                        case 19
                            sm(find(fnan==noJ)+1,[2,4])={'U'}; %Cannot exist
                        case 20
                            sm(find(fnan==noJ)+1,(2:4))=[{'S'},{'U'},{'S'}]; %Cannot exist
                        case 21
                            sm(find(fnan==noJ)+1,(3:4))=[{'U'},{'S'}];
                            
                        case 22
                            sm(find(fnan==noJ)+1,2)={'U'}; %Cannot exist
                        case 23
                            sm(find(fnan==noJ)+1,3)={'U'}; %Cannot exist
                        case 24
                            sm(find(fnan==noJ)+1,4)={'U'};
                        case 25
                            sm(find(fnan==noJ)+1,[2,3])={'U'}; %Cannot exist
                        case 26
                            sm(find(fnan==noJ)+1,[2,4])=[{'S'},{'U'}];
                        case 27
                            sm(find(fnan==noJ)+1,(2:4))=[{'S'},{'S'},{'U'}];
                        case 28
                            sm(find(fnan==noJ)+1,(3:4))=[{'S'},{'U'}];
                    end
                    
                end
                
                cla(handles.stabmatrix)
                axes(handles.stabmatrix)
                [nm,nt]=size(sm);
                yn=ones(1,nt);
                xn=0;
                for j=1:nt
                    xx(j)=(xn-0.16)+0.16;
                    xn=xn+0.16;
                end
                for k=1:nm
                    yy=yn*(0.8-(k-1)/10);
                    
                    text(xx,yy,sm(k,:),'interpreter','latex','horiz','left','vert','middle','fontsize',10)
                end
                
                title('Stability Plot')
            end
            
            xlabel(handles.prm1)
            ylabel(handles.prm2)
        end
        handles.mpltopt=0;
        %Save Report
%         newfig_g=figure;
%         axesobject_g=copyobj(handles.trajectoryplot,newfig_g);
%         title('Operating Plot')
%         export_fig temp_fig/ExpandReport -pdf -r600 -append
%         close(newfig_g)
        
    case 'boa'
        basin_of_attraction;
        
    case 'pport'
        variables=set_variables(growth);
        out_x_end=[]; out_y_end=[]; out_z_end=[];
        xlns=[]; ylns=[]; zlns=[]; out_x=[]; out_y=[]; out_z=[];
        parameters=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time1 Ks1 Ks2 Ks3 Ks3c gamma0 gamma1 gamma2];
        
        %Check number of equations
        noeq=length(eqs);
        %For 3-dimensions
        is3D=get(handles.use_3v,'Value');
        
        %Options for ODE solver
        if jacobanaly==1
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[],'Jpattern',Jpat);
        else
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[]);
        end
        
        solver=char(strtrim(solver));
        
        %%Get initial conditions
        min_ic1=str2num(get(handles.min_p1,'String')); min_ic2=str2num(get(handles.min_p2,'String'));
        max_ic1=str2num(get(handles.max_p1,'String')); max_ic2=str2num(get(handles.max_p2,'String'));
        step_ic1=str2num(get(handles.step_p1,'String'));
        
        if is3D==1
            min_ic3=str2num(get(handles.min_p3,'String')); max_ic3=str2num(get(handles.max_p3,'String'));
        end
        
        noeq=length(eqs);
        
        str1=get(handles.simparam1,'String'); str2=get(handles.simparam2,'String');
        val1=get(handles.simparam1,'Value'); val2=get(handles.simparam2,'Value');
        param1=str1(val1,:); param2=str2(val2,:);
        
        if is3D==1
            str3=get(handles.simparam3,'String'); val3=get(handles.simparam3,'Value');
            param3=str2(val3,:); init3(val3)=1;
        end
        
        init2=zeros(1,noeq);
        init2(val1)=1; init2(val2)=1;
        
        fz=find(init2==0);
        if noeq==3
            s1i=str2num(get(handles.s1_init,'String')); x1i=str2num(get(handles.x1_init,'String'));
            x2i=str2num(get(handles.x2_init,'String')); s2i=[]; s3i=[]; x3i=[]
        elseif noeq==4
            s1i=str2num(get(handles.s1_init,'String')); s2i=str2num(get(handles.s2_init,'String'));
            x1i=str2num(get(handles.x1_init,'String')); x2i=str2num(get(handles.x2_init,'String')); s3i=[]; x3i=[];
        elseif noeq==5
            s1i=str2num(get(handles.s1_init,'String')); s2i=str2num(get(handles.s2_init,'String')); s3i=str2num(get(handles.s3_init,'String'));
            x1i=str2num(get(handles.x1_init,'String')); x2i=str2num(get(handles.x2_init,'String')); x3i=[];
        elseif noeq==6;
            s1i=str2num(get(handles.s1_init,'String')); s2i=str2num(get(handles.s2_init,'String')); s3i=str2num(get(handles.s3_init,'String'));
            x1i=str2num(get(handles.x1_init,'String')); x2i=str2num(get(handles.x2_init,'String')); x3i=str2num(get(handles.x3_init,'String'));
        end
        vrs={'s1i','x1i','s2i','x2i','s3i','x3i'};
        
        for kj=1:6
            if any(fz==kj)
                init2(kj)=eval(char(vrs(kj)));
            end
        end
        
        %% Dynamics
        %Run solver - Dynamics at different initial conditions
        zt = length(0:0.01:time1);cn=1;
        %Random or Fixed points
        grof=get(handles.step_p2,'Value');
        rng(0,'twister');
        if is3D==0
            tt=time1*step_ic1*step_ic1; flag=3;
            letr=param1(1);
            Sarray=[1,3,5];
            Xarray=[2,4,6];
            if strcmp(letr,'S')
                s1m=find(Sarray==val1);
                s2m=find(Sarray==val2);
                Sarray([s1m,s2m])=[];
                val3=Sarray;
            elseif strcmp(letr,'X')
                x1m=find(Xarray==val1);
                x2m=find(Xarray==val2);
                Xarray([x1m,x2m])=[];
                val3=Xarray;
            end
            
            p3=str1(val3,:);
            strp3=['handles.',lower(p3),'_init'];
            init_val3=str2num(get(eval(strp3),'String'));
            p3=str1(val3,:);
            strp3=['handles.',lower(p3),'_init'];
            init_val3=str2num(get(eval(strp3),'String'));
            
            
            if grof==1
                lsk1=linspace(min_ic1,max_ic1,step_ic1);
                lsk2=linspace(min_ic2,max_ic2,step_ic1);
                lsk3=linspace(init_val3,init_val3,step_ic1);
                xlns=repmat(lsk1,1,length(lsk1)); ynl=repmat(lsk2,length(lsk2),1);
                ylns=reshape(ynl,1,numel(ynl));
                znl=repmat(lsk3,length(lsk3),1);
                zlns=reshape(znl,1,numel(znl));
                
            elseif grof==2
                lsk1a = (max_ic1-min_ic1).*rand(1,step_ic1) + min_ic1; lsk1=sort(lsk1a);
                lsk2a = (max_ic2-min_ic2).*rand(1,step_ic1) + min_ic2; lsk2=sort(lsk2a);
                lsk3=linspace(init_val3,init_val3,step_ic1);
                xlns=repmat(lsk1,1,length(lsk1)); ynl=repmat(lsk2,length(lsk2),1);
                ylns=reshape(ynl,1,numel(ynl));
                znl=repmat(lsk3,length(lsk3),1);
                zlns=reshape(znl,1,numel(znl));
            end
            
        elseif is3D==1
            tt=time1*step_ic1*step_ic1*step_ic1; flag=3;
            if grof==1
                lsk1=linspace(min_ic1,max_ic1,step_ic1);
                lsk2=linspace(min_ic2,max_ic2,step_ic1);
                lsk3=linspace(min_ic3,max_ic3,step_ic1);
                xlns=repmat(lsk1,1,length(lsk1)); ynl=repmat(lsk2,length(lsk2),1);
                ylns=reshape(ynl,1,numel(ynl)); zl=repmat(lsk3,1,length(lsk3)); zl2a=repmat(lsk3,length(lsk3),1);
                zl2b=reshape(zl2a,1,numel(zl2a)); zlns=meshgrid(zl,zl2b);
                
            elseif grof==2
                lsk1a = (max_ic1-min_ic1).*rand(1,step_ic1) + min_ic1; lsk1=sort(lsk1a);
                lsk2a = (max_ic2-min_ic2).*rand(1,step_ic1) + min_ic2; lsk2=sort(lsk2a);
                lsk3a = (max_ic3-min_ic3).*rand(1,step_ic1) + min_ic3; lsk3=sort(lsk3a);
                xlns=repmat(lsk1,1,length(lsk1)); ynl=repmat(lsk2,length(lsk2),1);
                ylns=reshape(ynl,1,numel(ynl));zl=repmat(lsk3,1,length(lsk3)); zl2a=repmat(lsk3,length(lsk3),1);
                zl2b=reshape(zl2a,1,numel(zl2a)); zlns=meshgrid(zl,zl2b);
            end
        end
        
        
        set(handles.func_prog,'String','Running: Dynamics','ForegroundColor','r')
        if is3D==0
            for k=lsk1
                for kk=lsk2
                    
                    init2(val1)=k;
                    init2(val2)=kk;
                    
                    h=handles.timestamp;
                    h1=handles.progress;
                    cno=(cn-1)*time1;
                    eval(['[tout,yout]=',solver,'(@model_gen, [0:0.01:time1], init2, options, S1in, D, Y1, kdec1, Y2, kdec2, Y3, kdec3, km1, Ks1, km2, Ks2, km3, Ks3, Ks3c, KI2,gamma0, gamma1, gamma2, S2in, S3in,h,h1,tt,motif,flag,cno);']);
                    
                    out_x(1:zt,cn)=yout(:,val1);
                    out_y(1:zt,cn)=yout(:,val2);
                    out_z(1:zt,cn)=yout(:,val3);
                    
                    out_x_end(cn)=yout(end,val1);
                    out_y_end(cn)=yout(end,val2);
                    out_z_end(cn)=yout(end,val3);
                    
                    cn=cn+1;
                end
            end
            
        elseif is3D==1
            for k=lsk1
                for kk=lsk2
                    for kkk=lsk3
                        
                        init2(val1)=k;
                        init2(val2)=kk;
                        init2(val3)=kkk;
                        
                        h=handles.timestamp;
                        h1=handles.progress;
                        cno=(cn-1)*time1;
                        eval(['[tout,yout]=',solver,'(@model_gen, [0:0.01:time1], init2, options, S1in, D, Y1, kdec1, Y2, kdec2, Y3, kdec3, km1, Ks1, km2, Ks2, km3, Ks3, Ks3c, KI2,gamma0, gamma1, gamma2, S2in, S3in,h,h1,tt,motif,flag,cno);']);
                        
                        out_x(1:zt,cn)=yout(:,val1);
                        out_y(1:zt,cn)=yout(:,val2);
                        out_z(1:zt,cn)=yout(:,val3);
                        
                        out_x_end(cn)=yout(end,val1);
                        out_y_end(cn)=yout(end,val2);
                        out_z_end(cn)=yout(end,val3);
                        cn=cn+1;
                    end
                end
            end
        end
        
        set(handles.func_prog,'String','Completed: Dynamics','ForegroundColor',[0 0.6 1])
        axes(handles.trajectoryplot)
        
        %Separate branches
        rnd_out=round(out_x(end,:),5);
        un_rnd=unique(round(out_x(end,:),5));
        n=length(un_rnd);
        
        %Overlay plots?
        ovopt=get(handles.overlay,'Value');
        
        if n<7
            if ovopt==1
                hold on
                colr=jet(6);
                colre='y';
                mrk='-.';
            else
                hold off
                colr=lines(6);
                colre='g';
                mrk='-';
            end
        else
            if ovopt==1
                hold on
                colr=jet(n);
                colre='y';
                mrk='-.';
            else
                hold off
                colr=lines(n);
                colre='g';
                mrk='-';
            end
        end
        
        
        for k=1:n
            if is3D==0
                h=plot(out_x(:,rnd_out==un_rnd(k)),out_y(:,rnd_out==un_rnd(k)),'color',colr(k,:),'Linestyle',mrk);hold on
            elseif is3D==1
                h=plot3(out_x(:,rnd_out==un_rnd(k)),out_y(:,rnd_out==un_rnd(k)),out_z(:,rnd_out==un_rnd(k)),'color',colr(k,:),'Linestyle',mrk);hold on
            end
            set(h,'linewidth',2)
        end
        %Plot initial conditions
        if is3D==0
            plot(xlns,ylns,'ks','markersize',7,'markerfacecolor','r')
            %Plot steady-states (if available), otherwise plot final values
            plot(out_x_end,out_y_end,'ko','markersize',10,'markerfacecolor',colre)
            axis tight
            xlabel(param1); ylabel(param2)
            hold off
            box off
            
            %Save Report
            newfig_h=figure;
            axesobject_h=copyobj(handles.trajectoryplot,newfig_h);
            title('2D Phase Portrait')
            export_fig temp_fig/ExpandReport -pdf -r600 -append
            close(newfig_h)
            
        elseif is3D==1
            plot3(xlns,ylns,zlns,'ks','markersize',7,'markerfacecolor','r')
            %Plot steady-states (if available), otherwise plot final values
            plot3(out_x_end,out_y_end,out_z_end,'ko','markersize',10,'markerfacecolor',colre)
            axis tight
            xlabel(param1); ylabel(param2); zlabel(param3);
            hold off
            grid on
            box on
            set(handles.Dimensions,'Visible','off')
            %Save Report
            newfig_i=figure;
            axesobject_i=copyobj(handles.trajectoryplot,newfig_i);
            title('3D Phase Portrait')
            export_fig temp_fig/ExpandReport -pdf -r600 -append
            close(newfig_i)
            
        end
        
        set(handles.overlay,'Enable','on')
        handles.yout_phase=out_x; handles.yout_phase2=out_y; handles.yout_phase3=out_z;
        handles.yout_phase_end=out_x_end; handles.yout_phase_end2=out_y_end; handles.yout_phase_end3=out_z_end;
        handles.colr_phase=colr; handles.val3=val3; handles.mrk_phase=mrk;
        handles.colre_phase=colre;
        handles.xlns=xlns; handles.ylns=ylns; handles.zlns=zlns;
        set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
        set(gca,'uicontextmenu',handles.Dimensions)
        
end
handles.plothandle=get(gca,'Children');
guidata(hObject,handles);




