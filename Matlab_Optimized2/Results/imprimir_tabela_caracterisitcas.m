clc
for i=1:part_n
    fprintf('P_{%d} & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d & %d \\\\\n',...
        i,...
        th(i),...
        px(i,1),px(i,2),px(i,3),px(i,4),...
        py(i,1),py(i,2),py(i,3),py(i,4),...
        cortes_y(i,1),cortes_y(i,2),cortes_y(i,3),cortes_y(i,4),...
        d_or(i,1),d_or(i,2));
end