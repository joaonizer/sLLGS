function []=plot_Particles(px,py,d_or,dx,dy,cor,jj,rows,cols,angles,name,eps)
%figure(%'Position',[0 0 100*cols 100*rows*2], ...
%    'Name','Alocação das Partículas');
h=figure('Visible','on');
%W = cols/2; H = rows*2;
W = 10; H = 10;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
part_n=size(px,1);
text_font_size=30;
    plot_linewidth=6;
    axis_linewidth=6;
for i=1:part_n
    fill(25+px(i,:)+d_or(i,1),50+py(i,:)+d_or(i,2),cor(i,:),'EdgeColor','none','LineStyle','none')
    hold on
    text(25+d_or(i,1),25+d_or(i,2)+py(i,1),['$P_{' num2str(i) '}$'],'Interpreter','latex','fontsize',text_font_size,'FontWeight','bold','Color','black','HorizontalAlignment','Center');
    if angles(i)>0
        text(25+d_or(i,1),25+d_or(i,2),['$\rightarrow$'],...
        'Interpreter','latex','fontsize',text_font_size*2,...
        'Rotation',angles(i),'FontWeight','bold',...
        'Color','black','HorizontalAlignment','Center');
    else
        text(25+d_or(i,1),25+d_or(i,2),['$\rightarrow$'],...
        'Interpreter','latex','fontsize',text_font_size*2,...
        'Rotation',angles(i),'FontWeight','bold',...
        'Color','white','HorizontalAlignment','Center');
    end
        
end
set(gca,'fontsize',text_font_size,'linewidth',axis_linewidth);


ylim([-30+min(d_or(:,2)) 30+max(d_or(:,2))+2*max(max(py(:,1:2)))])
xlim([-30+min(d_or(:,1)) 30+max(d_or(:,1))+2*max(max(px(:,2:3)))])
xlabel('nm','Interpreter','latex')
ylabel('nm','Interpreter','latex')
daspect([1 1 1])

%title(['XOR Architecture'],'Interpreter','tex')
%sdf('P1');
%print('-dpdfwrite',['XOR_Architecture.pdf'])
if ~eps
%print('-dpng','-r300',[name '_Architecture_' num2str(jj) '.png'])
else
 %       print( '-depsc', [name '_Architecture_' num2str(jj) '.eps'])
    end
%close all
end