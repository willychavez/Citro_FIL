clc
clear all
%salva a raiz
raiz = pwd ;

ym= dlmread("YMedios_FIL.txt");

cd('Agrupados_FIL/FIL_Citrosuco_09_22_2017');

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
    if j==1 &&  costheta < 0.991
      A1 = [A1;i];
    endif
    if j==2 &&  costheta < 0.995
      A2 = [A2;i];
    endif
    if j==3 &&  costheta < 0.990
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

cd([raiz '/Medidas']);

Espectros_1=[x,YYcorrigido{1}];
Espectros_2=[x,YYcorrigido{2}];
Espectros_3=[x,YYcorrigido{3}];
dlmwrite ('EspectrosMMO_exclusao1_ASS_agrupados_FIL.dat', Espectros_1, 'delimiter','\t', 'precision',4)
dlmwrite ('EspectrosMMO_exclusao1_SADIA_agrupados_FIL.dat', Espectros_2, 'delimiter','\t', 'precision',4)
dlmwrite ('EspectrosMMO_exclusao1_SINT_agrupados_FIL.dat', Espectros_3, 'delimiter','\t', 'precision',4)


figure;plot(x,YYcorrigido{1})
figure;plot(x,YYcorrigido{2})
figure;plot(x,YYcorrigido{3})

cd([raiz '/weka_FIL']);

%% ARFF
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%nome com que o arquivo arff deve ser salvo!
%cd ..

fid = (fopen('FIL_Citrosuco_09_22_2017.arff','w'));
fprintf(fid,'@RELATION Diagnostico \n');
fprintf(fid,'\n');
%o primeiro atributo sera o nome da folha
fprintf(fid,'@ATTRIBUTE PLANTA string \n');
fprintf(fid,'@ATTRIBUTE FOLHA string \n');
fprintf(fid,'\n');
for i=1:size(x,1)
    fprintf(fid,'@ATTRIBUTE %s NUMERIC \n',num2str(x(i)));
end

fprintf(fid,'@ATTRIBUTE class {Assintomatica,Sadia,Sintomatica} \n');

fprintf(fid,'\n');
fprintf(fid,'@DATA \n');

    
        %fprintf(fid,'%d,',k);    
        for l=1:length(YYcorrigido{1}(1,:))
        fprintf(fid,'%i,',folhas{1}(l,1:end));
            for j=1:size(x,1)
                fprintf(fid,'%0.5f,',YYcorrigido{1}(j,l));
            end
        fprintf(fid,'Assintomatica');
        fprintf(fid,'\n');
        end                            

        %fprintf(fid,'%d,',k);    
        for l=1:length(YYcorrigido{2}(1,:))
        fprintf(fid,'%i,',folhas{2}(l,1:end));
            for j=1:size(x,1)
                fprintf(fid,'%0.5f,',YYcorrigido{2}(j,l));
            end
            fprintf(fid,'Sadia');
            fprintf(fid,'\n');
        end
    

        %fprintf(fid,'%d,',k);    
        for l=1:length(YYcorrigido{3}(1,:))
        fprintf(fid,'%i,',folhas{3}(l,1:end));
            for j=1:size(x,1)
                fprintf(fid,'%0.5f,',YYcorrigido{3}(j,l));
            end
            fprintf(fid,'Sintomatica');
            fprintf(fid,'\n');
        end 
fclose(fid);    

cd(raiz)


