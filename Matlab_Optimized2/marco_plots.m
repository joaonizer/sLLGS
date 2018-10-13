%%
text_font_size=30;
plot_linewidth=6;
axis_linewidth=6;
i=1;
angles=m(:,2,i)*90;
figure('Position',[50,50,1400,700]);
for j=1:50:length(angles)
    pause(0.01);
    subplot(1,3,1);
    fill(25+px(i,:)+d_or(i,1),50+py(i,:)+d_or(i,2),cor(i,:),'EdgeColor','none','LineStyle','none')
    hold on
    text(25+d_or(i,1),25+d_or(i,2)+py(i,1),['$P_{' num2str(i) '}$'],'Interpreter','latex','fontsize',text_font_size,'FontWeight','bold','Color','black','HorizontalAlignment','Center');
    title(['Tempo: ' num2str(j*time_step/1e-9) 'ns']);
    if angles(j)>0
        text(25+d_or(i,1),25+d_or(i,2),['$\rightarrow$'],...
            'Interpreter','latex','fontsize',text_font_size*2,...
            'Rotation',angles(j),'FontWeight','bold',...
            'Color','black','HorizontalAlignment','Center');
    else
        text(25+d_or(i,1),25+d_or(i,2),['$\rightarrow$'],...
            'Interpreter','latex','fontsize',text_font_size*2,...
            'Rotation',angles(j),'FontWeight','bold',...
            'Color','white','HorizontalAlignment','Center');
    end
    xlabel('nm','Interpreter','latex')
    ylabel('nm','Interpreter','latex')
    daspect([1 1 1]);
    hold off
    
    subplot(1,3,2)
    plot(h_app(1:j,1,1)*1000,m(1:j,1,1)*Ms,'LineWidth',2);
    hold on
    plot(h_app(1:j,2,1)*1000,m(1:j,2,1)*Ms,'LineWidth',2);
    xlabel('H (mT)');
    ylabel('M (A/m)');
    xlim([min(min(h_app(:,1:2,1)))*1000-10 max(max(h_app(:,1:2,1)))*1000+10]);
    ylim([min(min(m(:,1:2,1)))*Ms-1000 max(max(m(:,1:2,1)))*Ms+1000]);
    legend('M_x \times H_x','M_y \times H_y');
    hold off
    
    subplot(1,3,3)
    plot([0:(j-1)]*time_step/1e-9,h_app(1:j,1,1)*1000,'LineWidth',2);
    hold on;
    plot([0:(j-1)]*time_step/1e-9,h_app(1:j,2,1)*1000,'LineWidth',2);
    xlabel('Tempo (ns)');
    ylabel('H_x (mT)');
    xlim([0 N*time_step/1e-9]);
    ylim([min(min(h_app(:,1:2,1)))*1000-10 max(max(h_app(:,1,1)))*1000+10]);
    %sdf('P1')
    hold off
    
end

%%
clc
Ex=-trapz(h_app(:,1,1),m(:,1,1)*Ms)*V(1)/q;
Ey=-trapz(h_app(:,2,1),m(:,2,1)*Ms)*V(1)/q;
fprintf('Energia dissipada: %2.2f meV\n',(Ex+Ey)*1000);
%% area artigo
base=325000;
altura=500*0.7977*1000;
altura=500e-3;
area=base*altura/2;
volume=2*pi*(10e-9)^2;
eex=area*volume/q