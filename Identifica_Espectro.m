clear all
clc
raiz = pwd;



%identifica os espectros das folhas
cd([raiz '/Numeros_folhas_W']);
files=dir;
b1=[];
b2=[];
for k=1:(length(files)-2)
  cd(files(k+2).name);
  arquivos=dir('*.txt'); 
 for j=1:length(arquivos)
    m=importdata(arquivos(j).name,'\t');
    b1=[b1;m];
end
n_folhas{k}(:,1)=b1;
b1=[];
cd ..
end

clear k j

for k=1:(length(files)-2)
   for i=1:length(n_folhas{k}())
      bi = [i*ones(n_folhas{k}(i,1),1), [1:n_folhas{k}(i,1)]'];
      b2=[b2;bi];
   end
   folhas{k}(:,1:2)=b2;
   b2=[];
end

clear k i 
%carrega numero dos espectros para o calculo do beta
cd([raiz '/Folhas_beta_Eq1/compoe_o_beta']);
beta_espectros=[];
n1=[];
n2=[];
b3=[];
arquivos1=dir('*.txt');

for j=1:length(arquivos1)
  m1=importdata(arquivos1(j).name,'\t');
  beta_folhas{j}(:,1:2)=m1;
end

clear k j
for k=1:3
  for j=1:size(beta_folhas{k},1)
    for i=1:size(folhas{k},1)
      if folhas{k}(i,1) == beta_folhas{k}(j,1) && folhas{k}(i,2) == beta_folhas{k}(j,2)
      n1 = i;
      beta_espectros = [beta_espectros;n1];
      endif
    end         
  end
 beta{k}(:,1)=beta_espectros;
 beta_espectros=[]; 
end
cd([raiz '/E_beta_FIL/compoe_o_beta']);
dlmwrite("Assintomatica.txt",beta{1})
dlmwrite("Sadia.txt",beta{2})
dlmwrite("Sintomatica.txt",beta{3})
cd ..