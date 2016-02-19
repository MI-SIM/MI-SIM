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

%Clear figure axes and text information
if ovra==0
    try
        cla(handles.solutionplot)
        delete(handles.legend1)
    catch
        try
            for k=1:length(handles.hsub)
                cla(handles.hsub(k))
                delete(handles.hlegend(k))
            end
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
        
        %Substitute numerical values into equations
        set(handles.func_prog,'String','Running: Initialisation','ForegroundColor','r','Visible','on')
        set(handles.progress,'String',['Progress: ',num2str(0),'%'])
        drawnow
        for k=1:noeq
            eval(['eq',num2str(k),'_numerical=subs(eqs(',num2str(k),'),variables,parameters);']);
        end
        
        %% Fixed points solutions
        
        if noeq==3
            syms S1 X1 X2 real
            assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(X2>=0); %Make assumptions (X,S real non-negative)
            %Use mupad to get all solutions
            sol=feval(symengine,'solve',[eq1_numerical==0 eq2_numerical==0 eq3_numerical==0],[S1,X1,X2]);
            temp=vpa(sol);
            members=children(children(temp));
            solutions = cellfun(@(c)c(:),members,'UniformOutput',false);
            num_sol=children(solutions{2});
            
            %Check for valid imaginary parts
            s_imag=sum(any(imag(double(vpa(num_sol)))>1e-12));
            if s_imag==0
                doub_sol=real(double(vpa(num_sol)));
            else
                doub_sol=double(vpa(num_sol));
            end
            set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
            set(handles.s2_check,'Enable','off','value',0); set(handles.x2_check,'Enable','on','value',1);
            set(handles.s3_check,'Enable','off','value',0); set(handles.x3_check,'Enable','off','value',0);
        elseif noeq==4
            syms S1 X1 S2 X2 real
            assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(S2>=0); assumeAlso(X2>=0); %Make assumptions (X,S real non-negative)
            %Use mupad to get all solutions
            sol=feval(symengine,'solve',[eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0],[S1,X1,S2,X2]);
            temp=vpa(sol);
            members=children(children(temp));
            solutions = cellfun(@(c)c(:),members,'UniformOutput',false);
            num_sol=children(solutions{2});
            
            %Check for valid imaginary parts
            s_imag=sum(any(imag(double(vpa(num_sol)))>1e-12));
            if s_imag==0
                doub_sol=real(double(vpa(num_sol)));
            else
                doub_sol=double(vpa(num_sol));
            end
            set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
            set(handles.s2_check,'Enable','on','value',1); set(handles.x2_check,'Enable','on','value',1);
            set(handles.s3_check,'Enable','off','value',0); set(handles.x3_check,'Enable','off','value',0);
        elseif noeq==5
            syms S1 X1 S2 X2 S3 real
            assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(S2>=0); assumeAlso(X2>=0); assumeAlso(S3>=0); %Make assumptions (X,S real non-negative)
            %Use mupad to get all solutions
            sol=feval(symengine,'solve',[eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0],[S1,X1,S2,X2,S3]);
            temp=vpa(sol);
            members=children(children(temp));
            solutions = cellfun(@(c)c(:),members,'UniformOutput',false);
            num_sol=children(solutions{2});
            
            %Check for valid imaginary parts
            s_imag=sum(any(imag(double(vpa(num_sol)))>1e-12));
            if s_imag==0
                doub_sol=real(double(vpa(num_sol)));
            else
                doub_sol=double(vpa(num_sol));
            end
            set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
            set(handles.s2_check,'Enable','on','value',1); set(handles.x2_check,'Enable','on','value',1);
            set(handles.s3_check,'Enable','on','value',1); set(handles.x3_check,'Enable','off','value',0);
        elseif noeq==6
            syms S1 X1 S2 X2 S3 X3 real
            assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(S2>=0); assumeAlso(X2>=0); assumeAlso(S3>=0); assumeAlso(X3>=0); %Make assumptions (X,S real non-negative)
            %Use mupad to get all solutions
            sol=feval(symengine,'solve',[eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0 eq6_numerical==0],[S1,X1,S2,X2,S3,X3]);
            temp=vpa(sol);
            members=children(children(temp));
            solutions = cellfun(@(c)c(:),members,'UniformOutput',false);
            num_sol=children(solutions{2});
            
            %Check for valid imaginary parts
            s_imag=sum(any(imag(double(vpa(num_sol)))>1e-12));
            if s_imag==0
                doub_sol=real(double(vpa(num_sol)));
            else
                doub_sol=double(vpa(num_sol));
            end
            set(handles.s1_check,'Enable','on','value',1); set(handles.x1_check,'Enable','on','value',1);
            set(handles.s2_check,'Enable','on','value',1); set(handles.x2_check,'Enable','on','value',1);
            set(handles.s3_check,'Enable','on','value',1); set(handles.x3_check,'Enable','on','value',1);
        end
        
        %Remove invalid fixed points (fp<0)
        [ii,jj]=find(doub_sol<0);
        iu = unique(ii);
        doub_sol(iu,:)=[];
        
        format shorteng
        format compact
        fixed_numerical1=doub_sol;
        
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
        lfn=length(fixed_numerical1(:,1));
        if lfn<=5
            for k=1:5
                if k<=lfn
                    eval(['set(handles.fp',num2str(k),',','''Visible'',''on'')'])
                else
                    eval(['set(handles.fp',num2str(k),',','''Visible'',''off'')'])
                end
            end
        else
            for k=1:5
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
            
            warning off
            set(handles.progress,'String',['Progress: ',num2str(0),'%'])
            syms S1 X1 X2 real
            eval(['syms ',param1]); eval(['syms ',param2]);
            %Make assumptions
            assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(X2>=0); 
            r1=str2double(get(handles.min_p1,'String')); r2=str2double(get(handles.max_p1,'String')); s1a=str2double(get(handles.step_p1,'String'));
            r3=str2double(get(handles.min_p2,'String')); r4=str2double(get(handles.max_p2,'String')); s2a=str2double(get(handles.step_p2,'String'));
            range1=linspace(r1,r2,s1a);
            range2=linspace(r3,r4,s2a);
            %Solve equations at equilibrium
            
            %Get Jacobian of system
            sym_vars = symvar(eqs);
            var_new=[]; var_new_char=[];
            for k=1:length(sym_vars)
                if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'X2') || strcmp(char(sym_vars(k)),strtrim(param1)) || strcmp(char(sym_vars(k)),strtrim(param2))
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
            jac_sys=jacobian(vpa(eqs),sym({'S1','X1','X2'}));
            
            xln=repmat(1:s1a,1,s2a);
            ylm=repmat(1:s1a,s2a,1);
            yln=reshape(ylm,1,s1a*s2a);
            
            jacs=matlabFunction(jac_sys);
            set(handles.func_prog,'String','Running: Multiple op. points - Mupad solver','ForegroundColor','r')
            drawnow
            sol=feval(symengine,'solve',[eq1_numerical==0 eq2_numerical==0 eq3_numerical==0],[S1,X1,X2]);
            eval(['syms ',param1]);
            eval(['syms ',param2]);

            drange=range2;
            srange=range1;
            [X,Y]=meshgrid(srange,drange);
            result_o = cell(s1a,s2a);
            parfor k=1:s2a
                temp = vpa( subs(sol,{param1,param2},{X(k,:),Y(k,:)}) );
                members = children(temp); %Find children of result
                solutions = cellfun(@(c)c(2),members,'UniformOutput',false); %Find numerical solutions
                r={};
                try
                    result_o(k,:) = cellfun(@symToDouble,solutions,'UniformOutput',false);
                catch
                    for kn=1:s2a
                        r{kn}={[X(k,kn);zeros(2,1)]};
                    end
                    result_o(k,:) = r;
                end
 
            end
          
            rsl=cellfun(@cell2mat,result_o,'UniformOutput',false);
            [rslw,rsll]=size(rsl);
            total_steps=s1a*s2a;
            no_steps=zeros(s1a,s2a);
            
            set(handles.func_prog,'String','Running: Multiple op. points - Calculating stability','ForegroundColor','r')
            drawnow
            for kk=1:rslw
                for ii=1:rsll
                    ind=ii+s2a*(kk-1);
                    total_tim=(ind)*100/(total_steps);
                    set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                    c = clock; % [year month day hour minute seconds]
                    TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                    set(handles.timestamp,'String',['Time: ',TSMP])
                    drawnow
                    Zr=rsl{kk,ii};
                    nS=zeros(3,numel(Zr)/3);
                    for i=1:numel(Zr)/3
                        nS(:,i)=Zr(:,i);
                        jac_P=jacs(range2(ylm(ii,kk)),nS(1,i),nS(2,i),nS(3,i));
                        cpoly=charpoly(jac_P);
                        
                        eigen=roots(cpoly);
                        
                        if all(real(eigen<=0))==1
                            stab_k=1;
                        else
                            stab_k=0;
                        end
                        %Get existence and stability at point
                        fpn=[nS(2,i),nS(3,i)];
                        jk=0;
                        if nS(1,i)>=0 && nS(2,i) >= 0 && nS(3,i)>=0
                            
                            if nnz(fpn)==2
                                jk=4;
                            elseif nnz(fpn)==0
                                jk=1;
                            elseif nnz(fpn)==1
                                nzu=find(fpn>0);
                                if nzu==1
                                    jk=2;
                                elseif nzu==2
                                    jk=3;
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
            
            axes(handles.trajectoryplot)
            contourf(range1,range2,rJ','LineStyle','none','LevelStep',1)
            axis tight
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
            
            set(handles.func_prog,'String','Completed: Multiple operating points','ForegroundColor',[0 0.6 1])
            guidata(hObject,handles);

            set(handles.func_prog,'String','Completed: Multiple operating points','ForegroundColor',[0 0.6 1])
            guidata(hObject,handles);
            
            %Four equations
        elseif noeq==4
            warning off
            set(handles.progress,'String',['Progress: ',num2str(0),'%'])
            syms S1 X1 S2 X2 real
            eval(['syms ',param1]); eval(['syms ',param2]);
            %Make assumptions
            assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(S2>=0); assumeAlso(X2>=0);
            r1=str2double(get(handles.min_p1,'String')); r2=str2double(get(handles.max_p1,'String')); s1a=str2double(get(handles.step_p1,'String'));
            r3=str2double(get(handles.min_p2,'String')); r4=str2double(get(handles.max_p2,'String')); s2a=str2double(get(handles.step_p2,'String'));
            range1=linspace(r1,r2,s1a);
            range2=linspace(r3,r4,s2a);
            %Solve equations at equilibrium
            
            %Get Jacobian of system
            sym_vars = symvar(eqs);
            var_new=[]; var_new_char=[];
            for k=1:length(sym_vars)
                if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'S2') || strcmp(char(sym_vars(k)),'X2') || strcmp(char(sym_vars(k)),strtrim(param1)) || strcmp(char(sym_vars(k)),strtrim(param2))
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
            
            xln=repmat(1:s1a,1,s2a);
            ylm=repmat(1:s1a,s2a,1);
            yln=reshape(ylm,1,s1a*s2a);
            
            jacs=matlabFunction(jac_sys);
            set(handles.func_prog,'String','Running: Multiple op. points - Mupad solver','ForegroundColor','r')
            drawnow
            sol=feval(symengine,'solve',[eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0],[S1,X1,S2,X2]);
            eval(['syms ',param1]);
            eval(['syms ',param2]);

            drange=range2;
            srange=range1;
            [X,Y]=meshgrid(srange,drange);
            result_o = cell(s1a,s2a);
            parfor k=1:s2a
  
                temp = vpa( subs(sol,{param1,param2},{X(k,:),Y(k,:)}) );
                members = children(temp); %Find children of result
                solutions = cellfun(@(c)c(2),members,'UniformOutput',false); %Find numerical solutions
                r={};
                try
                    result_o(k,:) = cellfun(@symToDouble,solutions,'UniformOutput',false);
                catch
                    for kn=1:s2a
                        r{kn}={[X(k,kn);zeros(3,1)]};
                    end
                    result_o(k,:) = r;
                end
 
            end
            
            rsl=cellfun(@cell2mat,result_o,'UniformOutput',false);
            [rslw,rsll]=size(rsl);
            total_steps=s1a*s2a;
            no_steps=zeros(s1a,s2a);
            set(handles.func_prog,'String','Running: Multiple op. points - Calculating stability','ForegroundColor','r')
            drawnow
            for kk=1:rslw
                for ii=1:rsll
                    ind=ii+s2a*(kk-1);
                    total_tim=(ind)*100/(total_steps);
                    set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                    c = clock; % [year month day hour minute seconds]
                    TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                    set(handles.timestamp,'String',['Time: ',TSMP])
                    drawnow
                    Zr=rsl{kk,ii};
                    nS=zeros(4,numel(Zr)/4);
                    for i=1:numel(Zr)/4
                        nS(:,i)=Zr(:,i);
                        jac_P=jacs(range2(ylm(ii,kk)),nS(1,i),nS(3,i),nS(2,i),nS(4,i));
                        cpoly=charpoly(jac_P);
                        
                        eigen=roots(cpoly);
                        
                        if all(real(eigen<=0))==1
                            stab_k=1;
                        else
                            stab_k=0;
                        end
                        %Get existence and stability at point
                        fpn=[nS(2,i),nS(4,i)];
                        jk=0;
                        if nS(1,i)>=0 && nS(2,i) >= 0 && nS(3,i)>=0 && nS(4,i) >= 0
                            
                            if nnz(fpn)==2
                                jk=4;
                            elseif nnz(fpn)==0
                                jk=1;
                            elseif nnz(fpn)==1
                                nzu=find(fpn>0);
                                if nzu==1
                                    jk=2;
                                elseif nzu==2
                                    jk=3;
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
            
            axes(handles.trajectoryplot)
            contourf(range1,range2,rJ','LineStyle','none','LevelStep',1)
            axis tight
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
            
            set(handles.func_prog,'String','Completed: Multiple operating points','ForegroundColor',[0 0.6 1])
            guidata(hObject,handles);
            
            % Five equations
        elseif noeq==5
            warning off
            set(handles.progress,'String',['Progress: ',num2str(0),'%'])
            syms S1 X1 S2 X2 S3 real
            eval(['syms ',param1]); eval(['syms ',param2]);
            %Make assumptions
            assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(S2>=0); assumeAlso(X2>=0); assumeAlso(S3>=0); 
            r1=str2double(get(handles.min_p1,'String')); r2=str2double(get(handles.max_p1,'String')); s1a=str2double(get(handles.step_p1,'String'));
            r3=str2double(get(handles.min_p2,'String')); r4=str2double(get(handles.max_p2,'String')); s2a=str2double(get(handles.step_p2,'String'));
            range1=linspace(r1,r2,s1a);
            range2=linspace(r3,r4,s2a);
            %Solve equations at equilibrium
            
            %Get Jacobian of system
            sym_vars = symvar(eqs);
            var_new=[]; var_new_char=[];
            for k=1:length(sym_vars)
                if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'S2') || strcmp(char(sym_vars(k)),'X2') || strcmp(char(sym_vars(k)),'S3') || strcmp(char(sym_vars(k)),strtrim(param1)) || strcmp(char(sym_vars(k)),strtrim(param2))
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
            
            xln=repmat(1:s1a,1,s2a);
            ylm=repmat(1:s1a,s2a,1);
            yln=reshape(ylm,1,s1a*s2a);
            
            jacs=matlabFunction(jac_sys);
            set(handles.func_prog,'String','Running: Multiple op. points - Mupad solver','ForegroundColor','r')
            drawnow
            sol=feval(symengine,'solve',[eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0],[S1,X1,S2,X2,S3]);
            eval(['syms ',param1]);
            eval(['syms ',param2]);

            drange=range2;
            srange=range1;
            [X,Y]=meshgrid(srange,drange);
            result_o = cell(s1a,s2a);
            parfor k=1:s2a
  
                temp = vpa( subs(sol,{param1,param2},{X(k,:),Y(k,:)}) );
                members = children(temp); %Find children of result
                solutions = cellfun(@(c)c(2),members,'UniformOutput',false); %Find numerical solutions
                r={};
                try
                    result_o(k,:) = cellfun(@symToDouble,solutions,'UniformOutput',false);
                catch
                    for kn=1:s2a
                        r{kn}={[X(k,kn);zeros(4,1)]};
                    end
                    result_o(k,:) = r;
                end
 
            end
            
            rsl=cellfun(@cell2mat,result_o,'UniformOutput',false);
            [rslw,rsll]=size(rsl);
            total_steps=s1a*s2a;
            no_steps=zeros(s1a,s2a);
            set(handles.func_prog,'String','Running: Multiple op. points - Calculating stability','ForegroundColor','r')
            drawnow
            for kk=1:rslw
                for ii=1:rsll
                    ind=ii+s2a*(kk-1);
                    total_tim=(ind)*100/(total_steps);
                    set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                    c = clock; % [year month day hour minute seconds]
                    TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                    set(handles.timestamp,'String',['Time: ',TSMP])
                    drawnow
                    Zr=rsl{kk,ii};
                    nS=zeros(5,numel(Zr)/5);
                    for i=1:numel(Zr)/5
                        nS(:,i)=Zr(:,i);
                        jac_P=jacs(range2(ylm(ii,kk)),nS(1,i),nS(3,i),nS(5,i),nS(2,i),nS(4,i));
                        cpoly=charpoly(jac_P);
                        
                        eigen=roots(cpoly);
                        
                        if all(real(eigen<=0))==1
                            stab_k=1;
                        else
                            stab_k=0;
                        end
                        %Get existence and stability at point
                        fpn=[nS(2,i),nS(4,i)];
                        jk=0;
                        if nS(1,i)>=0 && nS(2,i) >= 0 && nS(3,i)>=0 && nS(4,i) >= 0 && nS(5,i)>=0
                            
                            if nnz(fpn)==2
                                jk=4;
                            elseif nnz(fpn)==0
                                jk=1;
                            elseif nnz(fpn)==1
                                nzu=find(fpn>0);
                                if nzu==1
                                    jk=2;
                                elseif nzu==2
                                    jk=3;
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
            
            axes(handles.trajectoryplot)
            contourf(range1,range2,rJ','LineStyle','none','LevelStep',1)
            axis tight
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
            
            set(handles.func_prog,'String','Completed: Multiple operating points','ForegroundColor',[0 0.6 1])
            guidata(hObject,handles);
            

            % Six Equations
        elseif noeq==6
            warning off
            set(handles.progress,'String',['Progress: ',num2str(0),'%'])
            syms S1 X1 S2 X2 S3 X3 real
            eval(['syms ',param1]); eval(['syms ',param2]);
            %Make assumptions
            assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(S2>=0); assumeAlso(X2>=0); assumeAlso(S3>=0); assumeAlso(X3>=0);
            r1=str2double(get(handles.min_p1,'String')); r2=str2double(get(handles.max_p1,'String')); s1a=str2double(get(handles.step_p1,'String'));
            r3=str2double(get(handles.min_p2,'String')); r4=str2double(get(handles.max_p2,'String')); s2a=str2double(get(handles.step_p2,'String'));
            range1=linspace(r1,r2,s1a);
            range2=linspace(r3,r4,s2a);
            %Solve equations at equilibrium
            
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
            drawnow
            sol=feval(symengine,'solve',[eq1_numerical==0 eq2_numerical==0 eq3_numerical==0 eq4_numerical==0 eq5_numerical==0 eq6_numerical==0],[S1,X1,S2,X2,S3,X3]);
            eval(['syms ',param1]);
            eval(['syms ',param2]);
          
            drange=range2;
            srange=range1;
            [X,Y]=meshgrid(srange,drange);
            result_o = cell(s1a,s2a);
            
            parfor k=1:s2a
  
                temp = vpa( subs(sol,{param1,param2},{X(k,:),Y(k,:)}) );
                members = children(temp); %Find children of result
                solutions = cellfun(@(c)c(2),members,'UniformOutput',false); %Find numerical solutions
                r={};
                try
                    result_o(k,:) = cellfun(@symToDouble,solutions,'UniformOutput',false);
                catch
                    for kn=1:s2a
                        r{kn}={[X(k,kn);zeros(5,1)]};
                    end
                    result_o(k,:) = r;
                end
 
            end
            
            rsl=cellfun(@cell2mat,result_o,'UniformOutput',false);
            [rslw,rsll]=size(rsl);
            total_steps=s1a*s2a;
            no_steps=zeros(s1a,s2a);
            set(handles.func_prog,'String','Running: Multiple op. points - Calculating stability','ForegroundColor','r')
            drawnow
            for kk=1:rslw
                for ii=1:rsll
                    ind=ii+s2a*(kk-1);
                    total_tim=(ind)*100/(total_steps);
                    set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                    c = clock; % [year month day hour minute seconds]
                    TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                    set(handles.timestamp,'String',['Time: ',TSMP])
                    drawnow
                    Zr=rsl{kk,ii};
                    nS=zeros(6,numel(Zr)/6);
                    for i=1:numel(Zr)/6
                        nS(:,i)=Zr(:,i);
                        jac_P=jacs(range2(ylm(ii,kk)),nS(1,i),nS(3,i),nS(5,i),nS(2,i),nS(4,i),nS(6,i));
                        cpoly=charpoly(jac_P);
                        
                        eigen=roots(cpoly);
                        
                        if all(real(eigen<=0))==1
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
            
            axes(handles.trajectoryplot)
            contourf(range1,range2,rJ','LineStyle','none','LevelStep',1)
            axis tight
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
            
            set(handles.func_prog,'String','Completed: Multiple operating points','ForegroundColor',[0 0.6 1])
            guidata(hObject,handles);
            
        end
        
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
%Get axes labels
handles.xlabelhandle=get(gca,'Xlabel');
handles.ylabelhandle=get(gca,'Ylabel');
handles.zlabelhandle=get(gca,'Zlabel');
guidata(hObject,handles)




