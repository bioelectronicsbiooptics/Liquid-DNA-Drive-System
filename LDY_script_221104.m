%% FASTQ Read
tic
[~,S1,~] = fastqread('WSD_1.fastq');
[~,S2,~] = fastqread('WSD_2.fastq');
S1S = [S1 S2];
"Fastq read time: "+toc
%% 원본의 정보
OligoLen = input("원본 sequence의 길이: "); %스크립트 실행하면 명령창에 Oligo의 길이 입력할 수 있음
% Seq = 'GCACTTACCTCAACATACCCAACCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCGCGTCTAACCTGGGATGCAC'; %원본 sequence 입력란
% WSD
Seq = 'GCACTTACCTCAACATACCCTCCTAGTCTTGGGTGCGTCTAACCTGGGATGCACGTACTGCCATCGACACTGGCTGTTGTTATACAATGAATAGTTGTTGCTAATTAGTCAGCGGTGTACACCTCACGCCTTCACCATGCACTCAGACGTACAACGGCGCTCACGGCTCCGCGTACCCATTACACCTCACGCCTTCG';
% WLD
% Seq = 'GCACTTACCTCAACATACCCTCCTAGTCTTGGGTGCGTCTAACCTGGGATGCACGTACTGCCATCGACACTGGCGTATTCTACTTAGTTCGAGAAGTACACCTCACGCCTTCACCATGCACTCAGACGTACAACGGCGCTCACGGCTCCGCGTACCCATTACACCTCACGCCTTCG';
%% fastq 길이 분석
Qz = S1S(~cellfun('isempty',S1S));
A = zeros(length(Qz),1);
for a= 1 : length(Qz)
    A(a) = length(Qz{a});
end
histogram(A)
title('길이 Histogram')
B = Qz(find(A==OligoLen));
"뒤집기 전 최빈값 sequence: "+mode(cell2mat(B'),1)
%% F primer 잡고 뒤집기
% F_primers = 'GCACTTACCTCAACATACCC';
% R_primers = 'GCGTCTAACCTGGGATGCAC';

% F_primers = 'GTACTGCCATCGACACTGGC';
% R_primers = 'GTACACCTCACGCCTTCACC';

% WSD, WLD
F_primers = 'GCACTTACCTCAACATACCC';
R_primers = 'CATTACACCTCACGCCTTCG';


Reverse = F_primers;
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
    if length(Qz{a}) >= OligoLen
        if sum(Qz{a}(OligoLen-19:OligoLen) == rev) >= length(rev) %primer 완전 일치
            Reverse = Qz{a}(1:OligoLen);
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
            S{a} = Qz{a}(1:OligoLen);
        end
    end
end
"F primer로 뒤집힌 sequence 개수: " + i
%% R primer 잡고 뒤집기
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
j = 0;
T = cell(1,length(S));
for a = 1 : length(S)
    if length(S{a}) >= OligoLen
        if sum(S{a}(1:20) == rev) >= length(rev) %primer 완전 일치
            Reverse = S{a}(1:OligoLen);
            A_F = strfind(Reverse,'A');
            T_F = strfind(Reverse,'T');
            G_F = strfind(Reverse,'G');
            C_F = strfind(Reverse,'C');
            S_comp = Reverse;
            S_comp(A_F) = 'T';
            S_comp(T_F) = 'A';
            S_comp(G_F) = 'C';
            S_comp(C_F) = 'G';
            T{a} = flip(S_comp);
            j = j + 1;
        else
            T{a} = S{a}(1:OligoLen);
        end
    end
end
"R primer로 뒤집힌 sequence 개수: " + j
T = T(~cellfun('isempty',T));

"뒤집은 후 최빈값 sequence: "+mode(cell2mat(T'),1)
%% Primer Sorting
PrimerSortingSeq = cell(length(T),1);
for a = 1 : length(T)
    W = T{a};
    if length(W) >= OligoLen
        if sum(W(1:20) == F_primers) >=20
            if sum(W(OligoLen-19:OligoLen) == R_primers) >= 20
                PrimerSortingSeq{a} = W(1:OligoLen);
            end
        end
    end
end
PrimerSortingSeq = PrimerSortingSeq(~cellfun('isempty',PrimerSortingSeq));

if sum(Seq==mode(cell2mat(PrimerSortingSeq),1)) == OligoLen
    "최빈값 원본 매칭 완료"
end
%% Codon Bar
W = cell2mat(PrimerSortingSeq);
for b = 1 : 100
    U(b,1) = length(strfind(W(:,b)','A')); % A가 포함되어 있는 개수
    U(b,2) = length(strfind(W(:,b)','T')); % T가 포함되어 있는 개수
    U(b,3) = length(strfind(W(:,b)','C')); % C가 포함되어 있는 개수
    U(b,4) = length(strfind(W(:,b)','G')); % G가 포함되어 있는 개수
            %     N{b,1} = [double(H(strfind(W(:,b)','A')))-33 double(H(strfind(W(:,b)','T')))-33 ...

        %         double(H(strfind(W(:,b)','C')))-33 double(H(strfind(W(:,b)','G')))-33];
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