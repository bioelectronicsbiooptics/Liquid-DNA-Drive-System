%% 검사파일 오픈 및 바이너리로 변경.
addpath('Original Sequence') % addpath 함수, 특정 폴더에 대한 경로 추가
filename = 1
fid = fopen(filename+".txt"); % encoded TXT sequnces file
[scan,count] = fscanf(fid,'%s');
orgseq = reshape(scan,[length(scan)/count count])';
trim_seq =orgseq(:,21:90); % 21+78
binSeq = zeros(size(trim_seq,1),size(trim_seq,2)*2);
d = {'00','01','10','11'};
for r = 1 : size(trim_seq,1)
    Nnum = strfind(trim_seq(r,:),'N');
    if  numel(Nnum) >  0
        for i = 1 : numel(Nnum)
            trim_seq(r,Nnum(i)) = 'A';
        end
    end
    [~,x] = ismember(trim_seq(r,:),'ATGC');
    chr_out = cell2mat(d(x));
    binSeq(r,:) = chr_out-'0';
end
% A = zeros(length(dN),1);
% for a = 1 : length(dN)
%     if sum(dN{a}(21:end) == orgseq(floor((numel(Temp_mat)*2/3)/s),21:length(dNd))) == length(dN{1}(21:end))
%         A(a) = 1;
%     end
% end
% B = zeros(length(xN),1);
% for a = 1 : length(xN)
%     if sum(xN{a}(21:end) == orgseq(end-1,21:length(xNd))) == length(xN{1}(21:end))
%         B(a) = 1;
%     end
% end
  
%% Consensus 한뒤.
C = zeros(length(SeqtoBin),1);
for a = 1 : length(SeqtoBin)
        if sum(SeqtoBin(a,1:156) == binSeq(a,:)) >=156
            C(a) = C(a)+1;
        end
end
C= find(C==0)
%% GRS_한뒤.
D = zeros(length(rsECC),1);
for a = 1 : length(rsECC)
        if sum(rsECC(a,:) == binSeq(a,:)) >=156
            D(a) = D(a)+1;
        end
end
D = find(D==0)
%% XOR_한뒤. 변수 XOR로 맞추기 (index 빼고) 146 140
XOR = GRS_XOR_matrix;
E = zeros(length(XOR),1);
for a = 1 : length(XOR)
        if sum(XOR(a,1:matsz) == binSeq(a,1:matsz)) >=matsz
            E(a) = E(a)+1;
        end
end
E = find(E==0)
 
%% Depthcut_한뒤.  (index 빼고)
% DC = DC;
F = zeros(length(DC),1);
for a = 1 : length(DC)
        if sum(DC(a,1:140) == binSeq(a,1:140)) >=140
            F(a) = F(a)+1;
        end
end
F = find(F==0)
%% Depthcut_XOR 한뒤.  (index 빼고)
DX = Depth_XOR_matrix;
G = zeros(length(DX),1);
for a = 1 : length(DX)
        if sum(DX(a,1:140) == binSeq(a,1:140)) >=140
            G(a) = G(a)+1;
        end
end
G = find(G==0)
%% LRS_Depthcut_한뒤. Step1 (index 빼고)
% RS_DC
H = zeros(length(RS_DC),1);
for a = 1 : length(RS_DC)
        if sum(RS_DC(a,1:140) == binSeq(a,1:140)) >=140
            H(a) = H(a)+1;
        end
end
H = find(H==0)
%% LRS_Depthcut_XOR 한뒤. Step2 (index 빼고)
LDX=RS_DC_XOR_matrix;
I = zeros(length(LDX),1);
for a = 1 : length(LDX)
        if sum(LDX(a,1:140) == binSeq(a,1:140)) >=140
            I(a) = I(a)+1;
        end
end
I = find(I==0)
%% LRS  Step 3(index 빼고)
LRS=LRS_matrix;
J = zeros(length(LRS),1);
for a = 1 : length(LRS)
        if sum(LRS(a,1:140) == binSeq(a,1:140)) >=140
            J(a) = J(a)+1;
        end
end
J = find(J==0)
%% LRS XOR Step 4(index 빼고)
% LDX=LRS_matrix
 
% C(floor((numel(Temp_mat)*2/3)/s)) = sum(A);
% C(end-1) = sum(B);
%%
 
 
  % csvwrite(filename+"_Depth.xcsv",Depth,C,D)
% csvwrite(filename+"_PerfectMatching.csv",C)