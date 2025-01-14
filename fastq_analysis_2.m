%% FASTQ Read
tic
[~,S1,Q1] = fastqread('LDY1_S1_L001_R1_001.fastq');
[~,S2,Q2] = fastqread('LDY1_S1_L001_R2_001.fastq');


S1S = [S1 S2];
S1Q = [Q1 Q2];
toc
clearvars -except S1S S1Q


for a = 1 : length(S1S)
    if mean(double(S1Q{a}-33)) < 30
        S1S{a} = {};
    end
end



Qz = S1S(~cellfun('isempty',S1S));
A = zeros(length(Qz),1);
for a= 1 : length(Qz)
    A(a) = length(Qz{a});
end
histogram(A)
B = Qz(find(A==100));
mode(cell2mat(B'),1)
%%
% A
F_primers = 'GCACTTACCTCAACATACCC';
R_primers = 'GCGTCTAACCTGGGATGCAC';

% B or C
% F_primers = 'GTACTGCCATCGACACTGGC';
% R_primers = 'GTACACCTCACGCCTTCACC';

% D
%  F_primers = 'ATGCACTCAGACGTACAACG';
%  R_primers = 'CATTACACCTCACGCCTTCG';

% B'
% F_primers = 'CGACACTGGCTACCAGTACG';
% R_primers = 'GTACACCTCACGCCTTCACC';

Reverse = R_primers;
A_F = strfind(Reverse,'A');
T_F = strfind(Reverse,'T');
G_F = strfind(Reverse,'G');
C_F = strfind(Reverse,'C');
S_comp = Reverse;
S_comp(A_F) = 'T';
S_comp(T_F) = 'A';
S_comp(G_F) = 'C';
S_comp(C_F) = 'G';
rev = flip(S_comp);
i = 0;
S = cell(1,length(Qz));
for a = 1 : length(Qz)
    if sum(Qz{a}(1:20) == rev) >= length(rev)-2
        Reverse = Qz{a};
        A_F = strfind(Reverse,'A');
        T_F = strfind(Reverse,'T');
        G_F = strfind(Reverse,'G');
        C_F = strfind(Reverse,'C');
        S_comp = Reverse;
        S_comp(A_F) = 'T';
        S_comp(T_F) = 'A';
        S_comp(G_F) = 'C';
        S_comp(C_F) = 'G';
        S{a} = flip(S_comp);
         i = i + 1;
    else
        S{a} = Qz{a};
         end
end



DNAmatrix = S(find(A==100));
mode(cell2mat(DNAmatrix'),1)
%% Codon Bar
W = cell2mat(DNAmatrix');
for b = 1 : 100
    U(b,1) = length(strfind(W(:,b)','A')); % A가 포함되어 있는 개수
    U(b,2) = length(strfind(W(:,b)','T')); % T가 포함되어 있는 개수
    U(b,3) = length(strfind(W(:,b)','C')); % C가 포함되어 있는 개수
    U(b,4) = length(strfind(W(:,b)','G')); % G가 포함되어 있는 개수
end



for c = 1 : length(U)
    K(c,1) = U(c,1)/sum(U(c,:))*100; % codon 비율 백분율로 표현
    K(c,2) = U(c,2)/sum(U(c,:))*100;
    K(c,3) = U(c,3)/sum(U(c,:))*100;
    K(c,4) = U(c,4)/sum(U(c,:))*100;
end
figure
bar(K,'stacked')
title("title")
legend('A','T','C','G')

