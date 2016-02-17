%NUMERICAL METHOD
%this code introduces a pertubation in the parameter space and sees if it remains
%at the fixed points.

%NOTE: the commented code introduces a perturbation in the opposite
%direction however this uses more memory and slows down the calculation
    if exist('handles.htxt')
        delete(handles.htxt)
    end
stabanaly=get(handles.lsanaly,'Value');
stabanaly2=get(handles.routhcrit,'Value');
jacobanaly=get(handles.jacobian_but,'Value');
set(handles.eigplot,'enable','off')
s={'-'}; s2={'-'};
%Get Jacobian of system
sym_vars = symvar(eqs);
var_new=[]; var_new_char=[];
for k=1:length(sym_vars)
    if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'S2') || strcmp(char(sym_vars(k)),'X2') || strcmp(char(sym_vars(k)),'S3') || strcmp(char(sym_vars(k)),'X3')
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
kj=1; ind=[];
for k=1:length(sym_vars)
    if strmatch(char(sym_vars(k)),var_order)>0
        ind(kj)=strmatch(char(sym_vars(k)),var_order); kj=kj+1;
    end
end
%Create string elements
ind1=sort(ind);
var_order1=var_order(ind1,:);

%Create symbolic elements
vv=[];
for kh=1:length(ind1);
    eval(['v',num2str(kh),'= sym(var_order1(',num2str(kh),',:));']);
    eval(['vv=[vv,v',num2str(kh),'];']);
end

jac_sys=jacobian(vpa(eqs),vv);

%Create JPattern
[rws,cols]=size(jac_sys);
nuE=numel(jac_sys);

zons=ones(nuE,1);

zzer=find(jac_sys==0);

zons(zzer)=0;

Jpat=reshape(zons,rws,cols);

set(handles.spmatrix,'Visible','on');
axes(handles.spmatrix)
spy(Jpat,'r.',30)
axis tight
nz = nnz(Jpat);
pct = 100 / nuE;
xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*pct));
set(gca,'Xticklabel',cellstr(var_order1),'Yticklabel',cellstr(var_order1),'Xtick',1:length(Jpat),'Ytick',1:length(Jpat))

axes(handles.txax)

htxt=text(1.4,0.5,'Jacobian Sparsity Matrix','interpreter','latex','horiz','left','vert','middle','fontsize',14);
handles.htxt=htxt;
if stabanaly==1 %Linear stability analysis
    for i=1:length(fixed_numerical1(:,1))
        %Options for ODE solver (Jacobian not for ODE15 or ODE45)
        if jacobanaly==1
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[],'JPattern',Jpat);
        else
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[]);
        end
        
        %Get perturbation
        prt_s = get(handles.pert_p,'String');
        
        if isempty(str2num(prt_s))
            msgbox('The value entered for perturbation was not valid so the default value (1e-5) was used','Error','error')
            set(handles.pert_p,'String',0.00001);
            prt = 0.00001;
        else
            prt = str2num(prt_s);
        end
        
        %initial is the fixed point coordinates
        if noeq ==3
            initial=double([fixed_numerical1(i,1) fixed_numerical1(i,2) fixed_numerical1(i,3)]);
            init = initial+[prt 0 prt];
        elseif noeq==4
            initial=double([fixed_numerical1(i,1) fixed_numerical1(i,2) fixed_numerical1(i,3) fixed_numerical1(i,4)]);
            init = initial+[prt prt prt prt];
        elseif noeq==5
            initial=double([fixed_numerical1(i,1) fixed_numerical1(i,2) fixed_numerical1(i,3) fixed_numerical1(i,4) fixed_numerical1(i,5)]);
            init = initial+[prt 0 prt 0 prt];
        elseif noeq==6
            initial=double([fixed_numerical1(i,1) fixed_numerical1(i,2) fixed_numerical1(i,3) fixed_numerical1(i,4) fixed_numerical1(i,5) fixed_numerical1(i,6)]);
            init = initial+[prt prt prt prt prt prt];
        end
        
        solver=char(strtrim(solver));
        set(handles.func_prog,'String',['Running: Stability Analysis ',num2str(i)],'ForegroundColor','r')
        if growth==1
            h=handles.timestamp;
            h1=handles.progress;
            tt=time1;
            flag=[];
            eval(['[~,yout2]=',solver,'(@model_gen, 0:0.01:time1, init, options, S1in, D, Y1, kdec1, Y2, kdec2, Y3, kdec3, km1, Ks1, km2, Ks2, km3, Ks3, Ks3c, KI2, gamma0,gamma1,gamma2, S2in, S3in,h,h1,tt,motif,flag);']);
              
            %[~,yout3]=ode23s(@model_gen, 0:0.01:100, abs(initial-[0.00001 0 0.00001 0]), options, Spin1, D1, Yp1, kdecp1, YH1, kdecH1, kmp1, Ksp1, kmH1, KsH1, KIH1);
            
        elseif growth==2
            [~,yout2]=ode23s(@four_mod2, [0:0.01:time1], init, options, Spin1, D, Yp1, kdecp1, YH1, kdecH1, kmp1, kmH1, Ksxp1, KsxH1, KIH1); %CHANGE FOR ADDITIONAL substrate 3
            %[~,yout3]=ode23s(@four_mod2, [0:0.01:time1], abs(initial-[0.00001 0 0.00001 0]), options, Spin1, D1, Yp1, kdecp1, YH1, kdecH1, kmp1, kmH1, Ksxp1, KsxH1, KIH1);
        end
        
        %calculate the difference between the trajectories
        dif1(:,i)=yout2(end,:)-initial;
        %dif2(:,i)=yout3(end,:)-initial;
        
        error=str2double(get(handles.error_val,'String'));
        if error<0 || isnumeric(error)==0
            set(handles.error_val,'String',1e-4);
            error=str2double(get(handles.error_val,'String'));
            msgbox('The value entered for Error was not valid so the default value was reset','Error','error')
        else
        end
        
        if norm(dif1(:,i))<error
            s(:,i)={'Stable'};
        else
            s(:,i)={'Unstable'};
        end
        set(handles.func_prog,'String',['Completed: Stability Analysis ',num2str(i)],'ForegroundColor',[0 0.6 1])
    end
end

if stabanaly2==1 %Routh-Hurwitz stability criterion
    
    for i=1:length(fixed_numerical1(:,1))
        jac_sysa=jac_sys;
        for kk=1:length(fixed_numerical1(i,:))
            var_n = fixed_numerical1(i,kk);
            jac_sysa=subs(jac_sysa,cellstr(var_order1(kk,:)),subs(var_n));
            
        end
        cpoly=charpoly(jac_sysa);
        eigen=roots(cpoly);
        if all(real(eigen)>0) && all(isreal(eigen))
            s2(:,i)={'Unstable node'};
        elseif all(real(eigen)<0) && all(isreal(eigen))
            s2(:,i)={'Stable node'};
        elseif any(real(eigen)>0) && any(real(eigen)<0) && all(isreal(eigen))
            s2(:,i)={'Unstable saddle point'};
        elseif all(~isreal(eigen)) && all(real(eigen)>0)
            s2(:,i)={'Unstable spiral'};
        elseif all(~isreal(eigen)) && all(real(eigen)<0)
            s2(:,i)={'Stable spiral'};
        elseif all(~isreal(eigen)) && any(real(eigen)==0)
            s2(:,i)={'Circle'};
        else
            s2(:,i)={'Unstable'};
        end
        
        qeig(i,:)=double(eigen);

        set(handles.eigplot,'enable','on')
        clear jac_sysa
        
    end
            handles.eigen=qeig;
end
guidata(hObject,handles)
