function []=plot_M_and_H(m,h_app,t,part_n,a,jj,cols,rows,cor,grid)

%% Plot Place
nn=size(grid,1);
mm=size(grid,2);
count=1;
for i=1:mm
    for j=1:nn
        if grid(j,i)
            plot_place(count)=i+(j-1)*mm;
            count=count+1;
        end
    end
end

%%
h=figure(1);
W = 4*cols; H = 2*rows;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
    text_font_size=18;
    plot_linewidth=1;
    axis_linewidth=2;
    for j=1:1%part_n
        %subplot(rows,cols,plot_place(j));
        plot(t,squeeze(m(:,1:3,j)),'linewidth',plot_linewidth); % Plota a MagnetizaÃ§Ã£o
        hold on
        plot(t,squeeze(h_app(:,1:2,j))/a,'--','linewidth',plot_linewidth); % Plota o Campo aplicado em X (1) normalizado por a
        ylim([-1.2 1.2]);
        if j==part_n
            title(['$Output P_{' num2str(j) '}$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size); % cosiderando como ultima particula
        elseif j==1
            title(['$P_{' num2str(j) '} - Y$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size);
        elseif j==2
            title(['$P_{' num2str(j) '} - \overline{X}$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size);
        elseif j==3
            title(['$P_{' num2str(j) '} - \overline{Y}$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size);
        elseif j==4
            title(['$P_{' num2str(j) '} - X$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size);
        else
            title(['$P_{' num2str(j) '}$'], 'color', cor(j,:),'Interpreter','latex');
        end
        set(gca,'fontsize',text_font_size,'linewidth',axis_linewidth);
        %xlim([0 0.5])
        %hl=legend('m_x','m_y','h_{app_x}','h_{app_y}');
        %set(hl,'Orientation','Horizontal','Location','Best')
    end
    %sdf('P1');
    
    %print( '-dpdfwrite', ['XOR_' num2str(jj) '.pdf'])
    print( '-dpng', '-r300' ,['XOR_' num2str(jj) '.png'])
    close all
end