function []=plot_M_and_H(m,h_app,t,part_n,a,jj,cols,rows,cor,grid,name,eps)

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
h=figure('Visible','off');
W = 4*cols; H = 2*rows;
set(h,'Units','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
set(h,'Position',[0,0,W,H])
    text_font_size=18;
    plot_linewidth=1;
    axis_linewidth=2;
    for j=1:part_n
        subplot(rows,cols,plot_place(j));
        plot(t,squeeze(m(:,1:3,j)),'linewidth',plot_linewidth); % Plota a MagnetizaÃ§Ã£o
        hold on
        plot(t,squeeze(h_app(:,1:2,j))/a,'--','linewidth',plot_linewidth); % Plota o Campo aplicado em X (1) normalizado por a
        if j==part_n
            title(['$P_{' num2str(j) '}$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size); % cosiderando como ultima particula
        elseif j==1
            title(['$P_{' num2str(j) '}$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size);
        elseif j==2
            title(['$P_{' num2str(j) '}$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size);
        elseif j==3
            title(['$P_{' num2str(j) '}$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size);
        elseif j==4
            title(['$P_{' num2str(j) '}$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size);
        else
            title(['$P_{' num2str(j) '}$'], 'color', cor(j,:),'Interpreter','latex','fontsize',text_font_size);
        end
        set(gca,'fontsize',text_font_size,'linewidth',axis_linewidth);
        %set(gca, 'LooseInset', get(gca,'TightInset'))
        ylim([-1.2 1.2]);
        xlim([0 max(t)]);
        %xlim([0 0.5])
        %hl=legend('m_x','m_y','h_{app_x}','h_{app_y}');
        %set(hl,'Orientation','Horizontal','Location','Best')
    end
    sdf('P1');
    set(gca, 'LooseInset', get(gca,'TightInset'))
    %print( '-dpdfwrite', ['XOR_' num2str(jj) '.pdf'])
    if ~eps
    print( '-dpng', '-r300' ,[name '_' num2str(jj) '.png'])
    else
        print( '-depsc', [name '_' num2str(jj) '.eps'])
    end
    %close all
end