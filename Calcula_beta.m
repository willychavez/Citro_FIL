clc
clear all
%salva a raiz
raiz = pwd ;

ym= dlmread("YMedios_FIL.txt");

cd('DiasAgrupados_FIL');

%% CARREGA DADOS

files=dir;
data1=[];
for k=1:(length(files)-2)
cd(files(k+2).name);
arquivos=dir('*.txt'); % identifica os arquivos .txt na pasta

    for j=1:length(arquivos)
       data=dlmread(arquivos(j).name,'',2,0);     
       data1=[data1,data(:,2:(end-2))];
    end
    YYx{k}(:,:)=data1;
    data1=[];
    
cd ..
end

xx=data(:,1);

%Comprimentos de onda inicial x1 e final x2 do recorte
x1=444.2869;
x2=819.9622;
ind = find(xx>x1&xx<x2);

clear k j

%Corte
for i=1:length(YYx)
  for j=1:size(YYx{i},2)
    YYr{i}(:,j) = YYx{i}(ind,j);
  end
end 
x=xx(ind,:);

clear j i

%Retira o offset
for i=1:length(YYr)
  for j=1:size(YYr{i},2)
    m=min(YYr{i}(:,j));
    YYoff{i}(:,j)=YYr{i}(:,j)-m;
  end
end

clear j i

%Produto Escalar
A1 = [];
A2 = [];
A3 = [];

for j=1:size(ym,2)
norma_ym = norm(ym(:,j));
  for i=1:size(YYoff{j},2)
    norma_y = norm(YYoff{j}(:,i));
    costheta = dot(YYoff{j}(:,i),ym(:,j))/(norma_y * norma_ym);
    if j==1 &&  costheta < 0.996
      A1 = [A1;i];
    endif
    if j==2 &&  costheta < 0.9973
      A2 = [A2;i];
    endif
    if j==3 &&  costheta < 0.989
      A3 = [A3;i];
    endif
  end
end

clear k j

%identifica os espectros das folhas
cd([raiz '/Numeros_folhas_W']);
files1=dir;
b1=[];
b2=[];
for k=1:(length(files)-2)
cd(files1(k+2).name);
arquivos1=dir('*.txt'); 
 for j=1:length(arquivos1)
    m=importdata(arquivos1(j).name,'\t');
    b1=[b1;m];
end
n_folhas{k}(:,1)=b1;
b1=[];
cd ..
end

clear k j

for k=1:(length(files1)-2)
   for i=1:length(n_folhas{k}())
      bi = [i*ones(n_folhas{k}(i,1),1), [1:n_folhas{k}(i,1)]'];
      b2=[b2;bi];
   end
   folhas{k}(:,1:2)=b2;
   b2=[];
end

%Remove outliers
clear s j a

for s=1:length(A1)
   j = s - 1;
   a = A1(s) - j;
   YYoff{1}(:,a)=[];
   folhas{1}(a,:)=[];
end  

clear s j a

for s=1:length(A2)
   j = s - 1;
   a = A2(s) - j;
   YYoff{2}(:,a)=[];
   folhas{2}(a,:)=[];
end  


clear s j a 

for s=1:length(A3)
   j = s - 1;
   a = A3(s) - j;
   YYoff{3}(:,a)=[];
   folhas{3}(a,:)=[];
end

clear j i

%% NORMALIZA PELA AREA
for i=1:length(YYoff)
    for j=1:size(YYoff{i},2)
        %Calcula area embaixo da reta
        int{i}=trapz(x,YYoff{i});
        % Divido cada valor de intensidade do pico pela area abaixo da RETA
        YYcorrigido{i}(:,j)=YYoff{i}(:,j)/int{i}(1,j);
    end
end


%figure;plot(x,YYcorrigido{1})
%figure;plot(x,YYcorrigido{2})
%figure;plot(x,YYcorrigido{3})

clear k i 
%carrega numero dos espectros para o calculo do beta
cd([raiz '/Folhas_beta_FIL/Classificacao_FIL']);
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
 RR{k}(:,1)=beta_espectros;
 beta_espectros=[]; 
end

%carrega espectros para o calculo do beta
clear j k
for k = 1:3
  for j = 1: length(RR{k})
    Beta_Esp{k}(:,j) = YYcorrigido{k}(:,RR{k}(j));
  end
end  

clear j i
%Agrupa espectros em uma matriz
ym=[];
for i=1:3
  ym = [ym;Beta_Esp{i}'];
end

clear k i 
%carrega beta
cd([raiz '/Beta_FIL']);
arquivos1=dir('*.txt');

for j=1:length(arquivos1)
  n1=importdata(arquivos1(j).name,'\t');
  beta{j}(:,1)=n1;
end

%dlmwrite("EspectroFIL.txt",ym,'delimiter', '\t','precision',4)

ymcorrigido=[ones(size(ym,1),1) ym];

clear j i
R=[];
R2=[];
for j=1:3
  for i=1:size(ymcorrigido,1)
    r1=ymcorrigido(i,:)*beta{j};
    R= [R;r1];
  end
  R2(:,j)=R;
  R=[]; 
end

folhasid = [beta_folhas{1};beta_folhas{2};beta_folhas{3}];
resultados =[folhasid,R2];

cd([raiz '/Resultados']);
dlmwrite("Resultado_FIL.txt",resultados,'delimiter', '\t','precision',4)
cd ..