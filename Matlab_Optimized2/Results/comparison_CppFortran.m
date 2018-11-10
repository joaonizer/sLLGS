data_p=[
    num2str(px(1)) ' ' num2str(py(1)) ' '...
    num2str(px(2)) ' ' num2str(py(2)) ' '...
    num2str(px(3)) ' ' num2str(py(3)) ' '...
    num2str(px(4)) ' ' num2str(py(4)) ' '...
    num2str(th)...
    ];
tic;
for i=1:100
[flag,Nstr]=system(['echo '' ' data_p ' '' | ./Cpp/demag_cut']); % roda o executavel para gerar o OUT_demag3D2
end
tcpp=toc;

tic;
for i=1:100
[flag,Nstr]=system(['echo '' ' data_p ' '' | ./Fortran/demag3D3_noDAT']); % roda o executavel para gerar o OUT_demag3D2
end
tfortran=toc;

fprintf('Cpp:\t\t%f s\n',tcpp/100);
fprintf('Fortran:\t%f s\n',tfortran/100);

fprintf('Improvement:\t%2.2f %%\n',(1-tcpp/tfortran)*100);