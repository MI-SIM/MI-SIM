%plot_results - Plot results for algorithm output
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: SOFTWARE HOSTED SITE
% September 2015; Last revision: 05-Jan-2016

yout=handles.yout;
tout=handles.tout;
fixed_numerical=handles.fn;
%Check that data is available
if exist('yout')==0
    return
end

%For single-point (dynamic) analysis, check which variables to plot
s1_check=get(handles.s1_check,'Value');
x1_check=get(handles.x1_check,'Value');
s2_check=get(handles.s2_check,'Value');
x2_check=get(handles.x2_check,'Value');
s3_check=get(handles.s3_check,'Value');
x3_check=get(handles.x3_check,'Value');
labels={'s1','x1','s2','x2','s3','x3'};
checkind=[s1_check,x1_check,s2_check,x2_check,s3_check,x3_check];
fullcheck=sum(checkind);
indices_check = find(checkind);
checkenb=[];
if strcmp(get(handles.s1_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.x1_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.s2_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.x2_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.s3_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.x3_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end

ly=length(yout);
[lfn,sfn]=size(fixed_numerical);
kk=1;
for k=1:length(checkind)
    if checkenb(k)==1
        yout_n(:,k)=yout(:,kk);
        fixed_numerical_n(:,k)=fixed_numerical(:,kk);
        kk=kk+1;
    else
        yout_n(:,k)=zeros(ly,1);
        fixed_numerical_n(:,k)=zeros(lfn,1);
    end
end

%% Solutions plot
%Plots according to routine
switch which_plot
    case 'sol' %Single-point solutions
        switch handles.plotsolution
            
            case 'time' %Time series plot
                
                axes(handles.solutionplot);
                %reset the plot
                cla reset
                
                %Variables to plot
                legd=[];
                if s1_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,1),'linewidth',2,'color','b')
                    hold off
                    legd=[legd,{'S1'}];
                else
                end
                if x1_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,2),'linewidth',2,'color','r')
                    hold off
                    legd=[legd,{'X1'}];
                else
                end
                if s2_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,3),'linewidth',2,'color',[0 0.5 0])
                    hold off
                    legd=[legd,{'S2'}];
                else
                end
                if x2_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,4),'linewidth',2,'color',[1 0.6 0])
                    hold off
                    legd=[legd,{'X2'}];
                else
                end
                if s3_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,5),'linewidth',2,'color',[0.5 0 0.5])
                    hold off
                    legd=[legd,{'S3'}];
                else
                end
                if x3_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,6),'linewidth',2,'color',[0.302 0.745 0.933])
                    hold off
                    legd=[legd,{'X3'}];
                else
                end
                hold on
                %adds axis labels
                xlabel(handles.solutionplot,'Time (days)')
                ylabel(handles.solutionplot,'Concentration (kgCOD m^{-3})')
                ylimit=get(handles.solutionplot,'Ylim');
                set(handles.solutionplot,'Ylim',[0 ylimit(2)]);
                leg=legend(legd);
                hold off
                
                %Generate image of parameters
%                 newfig_par=figure;
%                 axesobject_par1=copyobj(handles.values.Children,newfig_par);
% 
%                 set(gca,'Visible','off')
%                 text(0,1,'Parameters','interpreter','latex','fontsize',14,'fontweight','bold')
%                 text(0.4,0.6,handles.eqtx,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
%                 text(0.4,0.2,handles.fnctx,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
%                 set(gcf,'Position',[1000 911 747 427])
%                 axesobject_par2=copyobj(handles.motif_image,newfig_par);
%                 set(axesobject_par2,'Position',[30 25 100 4])
%                 export_fig temp_fig/ExpandReport -pdf -r600 -append
%                 close(newfig_par)
%                 
%                 %Save Report
%                 newfig_a=figure;
%                 axesobject_a=copyobj([leg,handles.solutionplot],newfig_a);
%                 title('Plot of system solutions')
%                 export_fig temp_fig/ExpandReport -pdf -r600 -append
%                 close(newfig_a)
                
            case 'phase' %Phase plot
                
                axes(handles.solutionplot);
                %reset the plot
                cla reset
                
                if fullcheck == 0 || fullcheck == 1
                    msgbox('Too few variables selected, please select 2 or 3 variables','Error','error')
                    return
                elseif fullcheck == 4
                    msgbox('Too many variables selected, please select 2 or 3 variables','Error','error')
                    return
                elseif fullcheck == 2
                    plot(handles.solutionplot,yout_n(:,indices_check(1)),yout_n(:,indices_check(2)),'linewidth',2,'color','b')
                    ylimit=get(handles.solutionplot,'Ylim');
                    set(handles.solutionplot,'Ylim',[0 ylimit(2)]);
                    xlabel(handles.solutionplot,labels(indices_check(1)))
                    ylabel(handles.solutionplot,labels(indices_check(2)))
                elseif fullcheck == 3
                    plot3(handles.solutionplot,yout_n(:,indices_check(1)),yout_n(:,indices_check(2)),yout_n(:,indices_check(3)),'linewidth',2,'color','r')
                    grid on
                    xlabel(handles.solutionplot,labels(indices_check(1)))
                    ylabel(handles.solutionplot,labels(indices_check(2)))
                    zlabel(handles.solutionplot,labels(indices_check(3)))
                end
                
                %Save Report
                newfig_b=figure;
                axesobject_b=copyobj(handles.solutionsplot,newfig_b);
                title('Phase Plot')
                export_fig temp_fig/ExpandReport -pdf -r600 -append
                close(newfig_b)
                
            case 'subp' %Subplotting each variable separately
                
                if isfield(handles,'solutionplot') && ishandle(handles.solutionplot)
                    axes(handles.solutionplot)
                elseif isfield(handles,'hsub') 
                   hold on; axes(handles.hsub(1))
                end
                nstat=get(handles.normeig,'Value');

                chk=[s1_check,x1_check,s2_check,x2_check,s3_check,x3_check];
                vars={'S1','X1','S2','X2','S3','X3'};
                vc=vars(find(chk));
                
                N = sum(chk);
                
                if N == 3 || N == 4
                    sa = 2; sb=2;
                    
                elseif N == 5 || N ==6
                    sa = 3; sb=2;
                end
                list=lines(6);
                if handles.clcyc>6
                    list=lines(handles.clcyc);
                end
                
                if nstat==0
                    handles.clcyc=1;
                    for k=1:N
                        handles.hsub(k)=subplot(sa,sb,k);
                        plot(handles.hsub(k),tout,yout_n(:,k),'linewidth',2,'color',list(handles.clcyc,:))
                        ylim=get(gca,'Ylim');
                        handles.hlegend(k)=legend(vc(k)); hold off; set(gca,'Xlim',[0 tout(end)],'Ylim',[0 ylim(2)]);

                        handles.plotaxes_sub{k}=get(gca,'Children');
                            
                    end
                elseif nstat==1
                    handles.clcyc=handles.clcyc+1; hold all
                    for k=1:N

                        axes(handles.hsub(k))
                        hold on
                        
                        plot(handles.hsub(k),tout,yout_n(:,k),'linewidth',2,'color',list(handles.clcyc,:))
                        %legend(vc(k)); 
                        ylim=get(gca,'Ylim');
                        set(gca,'Xlim',[0 tout(end)],'Ylim',[0 ylim(2)]);

                        handles.plotaxes_sub{k}=get(gca,'Children');
                    end
                end
                
                %adds axis labels
                %suplabel('Time (days)');
                %suplabel('Concentration (kgCOD m^{-3})','y');
                
                %Save Report
                newfig_c=figure;
                axesobject_c=copyobj(handles.solutionplot,newfig_c);
                title('Individual subplots of solutions')
                export_fig temp_fig/ExpandReport -pdf -r600 -append
                close(newfig_c)
                
            case 'eig' %Eigenvalue plot
                axes(handles.solutionplot);
                %reset the plot
                cla reset
                nstat=get(handles.normeig,'Value');
                eig=handles.eigen';
                %Get real and imaginary parts
                if nstat==0
                    re = real(eig);
                    ime = imag(eig);
                    xlabel('Real')
                    ylabel('Img')
                else
                    re=normr(real(eig));
                    ime=imag(eig);
                end
                
                [m,n]=size(eig);
                col=num2cell(jet(n),2);
                markers = {'o','s','d','v','x','+','*','^','.','>','<','p','h'};
                mrk=markers(1:n)';
                lh=plot(re,ime,'.','markersize',10);
                set(lh,{'marker'},mrk,{'markerfacecolor'},col,'markeredgecolor',[0,0,0]);

                if nstat==0
                    xlabel('Real')
                    ylabel('Img')
                else
                    xlabel('Norm(Real)')
                    ylabel('Img')
                end
                
               %Legend
               for k=1:n
                   lb(k)=cellstr(['FP',num2str(k)]);
               end
               legend(lb,'location','Best')
               xl=get(gca,'Xlim');
               yl=get(gca,'Ylim');
               hold on
               hx=line(xl,[0 0]);
               hy=line([0 0],yl);
               set(hx,'linestyle','--','color','k')
               set(hy,'linestyle','--','color','k')
               
               %Save Report
               newfig_d=figure;
               axesobject_d=copyobj(handles.solutionplot,newfig_d);
               title('Eigenvalues Plot')
               export_fig temp_fig/ExpandReport -pdf -r600 -append
               close(newfig_d)
        end

    case 'traj'  %Trajectory plot
        pben=get(handles.plot_button,'enable');
        if strcmp(pben,'on')==0
            return
        end
        axes(handles.trajectoryplot);
        
        switch handles.plotdim
            case 'two' %2D plot
                
                %check if plot is to be overlayed
                multi=get(handles.overlay,'Value');
                if multi==0
                    cla reset
                else
                end
                
                %Check correct number of variables are selected
                if fullcheck == 0 || fullcheck == 1 || fullcheck == 3 || fullcheck == 4 || fullcheck==5
                    msgbox('Too few variables selected, please select 2 variables','Error','error')
                    return
                else
                    X_value=yout_n(:,indices_check(1));
                    Y_value=yout_n(:,indices_check(2));
                    hold on
                    %plot the fixed points
                    z=length(fixed_numerical(:,1));
                    list=hsv(z);
                    for i=1:z
                        plot(handles.trajectoryplot,fixed_numerical_n(i,indices_check(1)),fixed_numerical_n(i,indices_check(2)),'s','color',list(i,:),'MarkerFaceColor',list(i,:))
                        eval(['set(handles.fp',num2str(i),',''ForegroundColor'',list(',num2str(i),',:));'])
                    end
                    %plot the trajectory
                    plot(handles.trajectoryplot,X_value,Y_value,'linewidth',2)
                    xlabel(handles.trajectoryplot,labels(indices_check(1)))
                    ylabel(handles.trajectoryplot,labels(indices_check(2)))
                    %plot the start point (green dot)
                    plot(handles.trajectoryplot,yout_n(1,indices_check(1)),yout_n(1,indices_check(2)),'o','color','g','MarkerFaceColor','g')
                    %plot the end point (red dot)
                    plot(handles.trajectoryplot,yout_n(end,indices_check(1)),yout_n(end,indices_check(2)),'o','color','r','MarkerSize',10,'linewidth',2)
                    hold off
                end
                set(handles.overlay,'enable','on')
                
                %Save Report
                newfig_e=figure;
                axesobject_e=copyobj(handles.trajectoryplot,newfig_e);
                title('2D Trajectory Plot')
                export_fig temp_fig/ExpandReport -pdf -r600 -append
                close(newfig_e)
                
            case 'three' %3D plot
                %if user wants to overlay multiple trajectories onto the plot of
                %trajectories
                multi=get(handles.overlay,'Value');
                if multi==0
                    cla reset
                else
                    hold on
                end
                
                %Check correct number of variables are selected
                if fullcheck == 0 || fullcheck == 1 || fullcheck == 2 || fullcheck == 4 || fullcheck == 5
                    msgbox('Too few variables selected, please select 3 variables','Error','error')
                    return
                else
                    X_value=yout_n(:,indices_check(1));
                    Y_value=yout_n(:,indices_check(2));
                    Z_value=yout_n(:,indices_check(3));
                    
                    %plot the fixed points
                     z=length(fixed_numerical(:,1));
                    list=hsv(z);
                    for i=1:z
                        plot3(handles.trajectoryplot,squeeze(fixed_numerical_n(i,indices_check(1))),squeeze(fixed_numerical_n(i,indices_check(2))),squeeze(fixed_numerical_n(i,indices_check(3))),'s','color',list(i,:),'MarkerFaceColor',list(i,:))
                        eval(['set(handles.fp',num2str(i),',''ForegroundColor'',list(',num2str(i),',:));'])
                        hold on
                    end
                    %plot the trajectory
                    plot3(handles.trajectoryplot,X_value,Y_value,Z_value,'linewidth',2)
                    xlabel(handles.trajectoryplot,labels(indices_check(1)))
                    ylabel(handles.trajectoryplot,labels(indices_check(2)))
                    zlabel(handles.trajectoryplot,labels(indices_check(3)))
                    %plot the start point (green dot)
                    plot3(handles.trajectoryplot,yout_n(1,indices_check(1)),yout_n(1,indices_check(2)),yout_n(1,indices_check(3)),'o','color','g','MarkerFaceColor','g')
                    %plot the end point (red dot)
                    plot3(handles.trajectoryplot,yout_n(end,indices_check(1)),yout_n(end,indices_check(2)),yout_n(end,indices_check(3)),'o','color','r','MarkerSize',10,'linewidth',2)
                    grid on
                    hold off
                end
                set(handles.overlay,'enable','on')

                %Save Report
                newfig_f=figure;
                axesobject_f=copyobj(handles.trajectoryplot,newfig_f);
                title('3D Trajectory Plot')
                export_fig temp_fig/ExpandReport -pdf -r600 -append
                close(newfig_f)
                
        end
        
end
handles.plothandle=get(gca,'Children');
guidata(hObject,handles);


