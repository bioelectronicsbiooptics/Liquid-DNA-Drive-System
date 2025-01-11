%% FASTQ query
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Encoded data lists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. DNAdata Store Apple
% 2. JamesWatson Likes Orange
% 3. FrancisCrick Likes Grape
% 4. JohnVonNeumann Has Apple
% 5. DNAdata Store Apple Orange
% 6. DNAdata Store Parallelly [Apple, Orange, Grape]
% 7. We synthesize DNA.
% 8. We love DNA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Primer Number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. DNAdata
% 2. Storage
% 3. RandomAccess
% 4. Apple
% 5. Orange
% 6. Grape
% 7. JamesWatson
% 8. FrancisCrick
% 9. JohnVonNeumann
% 10. Store
% 11. Likes
% 12. Has
% 13. Parallelly
% 14. DNA.
% 15. synthesize, love
% 16. We
%%%%%%%%%%%%%%%%%%%%%%%%%%%% DNA Encoding lists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNAdata = 'TTCGTTCGTCGTTGATTGGTGCGCTCACGGCAGTCGACAATCGATTAGTAACAAACGGAGCCATGAGTTTGT';
% Store = 'TCCTCAGCCGATGAAATTCCAATTTGCACTAGTGGTCCCCGGACCTTGTACCATCCGTTTGACTGG';
% Apple = 'CTGTCCATAGCCTTGTTCGTCGGCGATGGTCACTCCATTATTCGCGCGGAAACGTAGTGAAGGTA';
% J_W = 'GTCCAGGCAAAGATCCAGTTAAGTGGTCTATAGTTGTTTAACGACAATTGTTTCTATACCCCTGACACCACCGTTAGGCTAAAGTG';
% Likes = 'TACCGCATCCTTATTCGAGCAAACCAGCAGCGAATTCCCCGGATACTCTGGTGCAAGCCAATGAAA';
% Orange = 'TGTATTTCCTTCGGTGCTCCTTGCCATCGACTCTTCATAGTTACTCACTTTCGACAACGGTCTGGTTT';
% F_C = 'ATCCTGCAAACGCATTTCCTACGCCATCGACTCTTGATGCTTGTTCGAAGTTGCTTGATCCTCTCGGGAATGCCTTTCCGAAGTTTCCA';
% Grape = 'AGCCTTGTGTCCATCAATCCTCGCCATCGACGATCCATTAGATTGTGCGCTATGGTTTGGCTAAT';
% J_V_N = 'TAGCCTCCAGAATGAAACGGAAGTTCTGAACTCTACCTTCTATACTCACCATAATCTATGGTATACTCTAGTGTTGTTCAAGCCAAACCGTGTGTA';
% Has = 'GAAGAGTTTAGCCACCTGGTAACAGAGATTCCTGGTAAGGCCAATTCGCGGTTATT';
% Parallelly = 'CAAGATTGTGGACGATTGGCAAAAATTAGTAGAAAAACCGAACTCGGAATTCCCCCCCCCGTTAGTTTGCAATGTTTCCGTCGGTTT';
% We = 'GCACTTACCTCAACATACCCTCCTAGTCTTGGGTGCGTCTAACCTGGGATGCAC';
% love = 'GTACTGCCATCGACACTGGCGTATTCTACTTAGTTCGAGAAGTACACCTCACGCCTTCACC';
% synthesize = 'GTACTGCCATCGACACTGGCTGTTGTTATACAATGAATAGTTGTTGCTAATTAGTCAGCGGTGTACACCTCACGCCTTCACC;
% DNA. = ATGCACTCAGACGTACAACGGCGCTCACGGCTCCGCGTACCCATTACACCTCACGCCTTCG;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear,clc
filename = input("Name of Fastq file?: ",'s');
DNA1 = input("First word: ",'s');
DNA2 = input("Second word: ",'s');
DNA3 = input("Third word: ",'s');
DNA4 = input("Fourth word(enter if not): ",'s');
PrimerEquality = input("Number of primer matching?: ");
DNAinf = [{DNA1} DNA2 DNA3 DNA4];

[~,Fastq1,~] = fastqread(filename+"_1.fastq");
[~,Fastq2,~] = fastqread(filename+"_2.fastq");
%%
Primer1 = cell(1,2);
Primer2 = cell(1,2);
Primer3 = cell(1,2);
Primer4 = cell(1,2);

if sum(ismember(DNAinf,'DNAdata')) > 0
    eval("Primer"+find(ismember(DNAinf,'DNAdata'))+"{1} = 'TTCGTTCGTCGTTGATTGGT';")
    eval("Primer"+find(ismember(DNAinf,'DNAdata'))+"{2} = 'AAACGGAGCCATGAGTTTGT';")
end
if sum(ismember(DNAinf,'Store')) > 0
    eval("Primer"+find(ismember(DNAinf,'Store'))+"{1} = 'TCCTCAGCCGATGAAATTCC';")
    eval("Primer"+find(ismember(DNAinf,'Store'))+"{2} = 'TGTACCATCCGTTTGACTGG';")
end
if sum(ismember(DNAinf,'Apple')) > 0
    eval("Primer"+find(ismember(DNAinf,'Apple'))+"{1} = 'CTGTCCATAGCCTTGTTCGT';")
    eval("Primer"+find(ismember(DNAinf,'Apple'))+"{2} = 'GCGGAAACGTAGTGAAGGTA';")
end
if sum(ismember(DNAinf,'J_W')) > 0
    eval("Primer"+find(ismember(DNAinf,'JamesWatson'))+"{1} = 'GTCCAGGCAAAGATCCAGTT';")
    eval("Primer"+find(ismember(DNAinf,'JamesWatson'))+"{2} = 'ACCACCGTTAGGCTAAAGTG';")
end
if sum(ismember(DNAinf,'Likes')) > 0
    eval("Primer"+find(ismember(DNAinf,'Likes'))+"{1} = 'TACCGCATCCTTATTCGAGC';")
    eval("Primer"+find(ismember(DNAinf,'Likes'))+"{2} = 'TCTGGTGCAAGCCAATGAAA';")
end
if sum(ismember(DNAinf,'Orange')) > 0
    eval("Primer"+find(ismember(DNAinf,'Orange'))+"{1} = 'TGTATTTCCTTCGGTGCTCC';")
    eval("Primer"+find(ismember(DNAinf,'Orange'))+"{2} = 'TTTCGACAACGGTCTGGTTT';")
end
if sum(ismember(DNAinf,'F_C')) > 0
    eval("Primer"+find(ismember(DNAinf,'FrancisCrick'))+"{1} = 'ATCCTGCAAACGCATTTCCT';")
    eval("Primer"+find(ismember(DNAinf,'FrancisCrick'))+"{2} = 'ATGCCTTTCCGAAGTTTCCA';")
end
if sum(ismember(DNAinf,'Grape')) > 0
    eval("Primer"+find(ismember(DNAinf,'Grape'))+"{1} = 'AGCCTTGTGTCCATCAATCC';")
    eval("Primer"+find(ismember(DNAinf,'Grape'))+"{2} = 'TGCGCTATGGTTTGGCTAAT';")
end
if sum(ismember(DNAinf,'J_V_N')) > 0
    eval("Primer"+find(ismember(DNAinf,'JohnVonNeumann'))+"{1} = 'TAGCCTCCAGAATGAAACGG';")
    eval("Primer"+find(ismember(DNAinf,'JohnVonNeumann'))+"{2} = 'TTCAAGCCAAACCGTGTGTA';")
end
if sum(ismember(DNAinf,'Has')) > 0
    eval("Primer"+find(ismember(DNAinf,'Has'))+"{1} = 'GAAGAGTTTAGCCACCTGGT';")
    eval("Primer"+find(ismember(DNAinf,'Has'))+"{2} = 'AAGGCCAATTCGCGGTTATT';")
end
if sum(ismember(DNAinf,'Parallelly')) > 0
    eval("Primer"+find(ismember(DNAinf,'Parallelly'))+"{1} = 'CAAGATTGTGGACGATTGGC';")
    eval("Primer"+find(ismember(DNAinf,'Parallelly'))+"{2} = 'TGCAATGTTTCCGTCGGTTT';")
end
if sum(ismember(DNAinf,'DNA.')) > 0
    eval("Primer"+find(ismember(DNAinf,'DNA.'))+"{1} = 'ATGCACTCAGACGTACAACG';")
    eval("Primer"+find(ismember(DNAinf,'DNA.'))+"{2} = 'CATTACACCTCACGCCTTCG';")
end
if sum(ismember(DNAinf,'synthesize')) > 0
    eval("Primer"+find(ismember(DNAinf,'synthesize'))+"{1} = 'GTACTGCCATCGACACTGGC';")
    eval("Primer"+find(ismember(DNAinf,'synthesize'))+"{2} = 'GTACACCTCACGCCTTCACC';")
end
if sum(ismember(DNAinf,'love')) > 0
    eval("Primer"+find(ismember(DNAinf,'love'))+"{1} = 'GTACTGCCATCGACACTGGC';")
    eval("Primer"+find(ismember(DNAinf,'love'))+"{2} = 'GTACACCTCACGCCTTCACC';")
end
if sum(ismember(DNAinf,'We')) > 0
    eval("Primer"+find(ismember(DNAinf,'We'))+"{1} = 'GCACTTACCTCAACATACCC';")
    eval("Primer"+find(ismember(DNAinf,'We'))+"{2} = 'GCGTCTAACCTGGGATGCAC';")
end

Fastq = [Fastq1 Fastq2];
if isempty(DNA4)
    Rv_seq = Primer3{2};
else
    Rv_seq = Primer4{2};
end
Complement_rv = seqrcomplement(Rv_seq);
%% Find Complementary Sequence
SortingFastq1 = cell(1,length(Fastq));
for a = 1 : length(Fastq)
    if sum(Fastq{a}(1:length(Complement_rv)) == Complement_rv) >= length(Complement_rv)-2
        SortingFastq1{a} = seqrcomplement(Fastq{a});
    else
        SortingFastq1{a} = Fastq{a};
    end
end
%% Primer Sorting
SortingFastq2 = cell(length(SortingFastq1),1);
for a = 1 : length(SortingFastq1)
    if sum(SortingFastq1{a}(1:20) == Primer1{1}) >= PrimerEquality
        SortingFastq2{a} = SortingFastq1{a};
    end
end
empty = ~cellfun('isempty',SortingFastq2);
SortingFastq2 = SortingFastq2(empty);

SortingFastq3 = cell(length(SortingFastq2),4);
PMidx = zeros(length(SortingFastq2),4);
for b = 1 : length(SortingFastq2)
    for c = 40 : length(SortingFastq2{1})
        if sum(SortingFastq2{b}(c-19:c) == Primer1{2}) >= PrimerEquality
            SortingFastq3{b,1} = SortingFastq2{b}(1:c);
            PMidx(b,1) = 1;
            break
        end
    end
    for d = 60 : length(SortingFastq2{1})
        if sum(SortingFastq2{b}(d-19:d) == Primer2{1}) >= PrimerEquality
            for e = 80 : length(SortingFastq2{1})
                if sum(SortingFastq2{b}(e-19:e) == Primer2{2}) >= PrimerEquality
                    SortingFastq3{b,2} = SortingFastq2{b}(d-19:e);
                    PMidx(b,2) = 1;
                break
                end
            end
            break
        end
    end
    for f = 100 : length(SortingFastq2{1})
        if sum(SortingFastq2{b}(f-19:f) == Primer3{1}) >= PrimerEquality
            for g = 120 : length(SortingFastq2{1})
                if sum(SortingFastq2{b}(g-19:g) == Primer3{2}) >= PrimerEquality
                    SortingFastq3{b,3} = SortingFastq2{b}(f-19:g);
                    PMidx(b,3) = 1;
                break
                end
            end
            break
        end
    end
    if ~isempty(DNA4)
        for h = 120 : length(SortingFastq2{1})
            if sum(SortingFastq2{b}(h-19:h) == Primer4{1}) >= PrimerEquality
                for i = 100 : length(SortingFastq2{1})
                    if sum(SortingFastq2{b}(i-19:i) == Primer4{2}) >= PrimerEquality
                        SortingFastq3{b,4} = SortingFastq2{b}(h-19:i);
                        PMidx(b,4) = 1;
                        break
                    end
                end
            end
            break
        end
    end
end

str1 = SortingFastq3(:,1);
str2 = SortingFastq3(:,2);
str3 = SortingFastq3(:,3);
str4 = SortingFastq3(:,4);
%% Primer perfect matching
PMsort = zeros(length(SortingFastq3),1);
for a = 1 : length(SortingFastq3)
    if ~isempty(str1(a))&&length(str1{a})==length(DNA1)
        if ~isempty(str2(a))&&length(str2{a})==length(DNA2)
            if ~isempty(str3(a))&&length(str3{a})==length(DNA3)
                PMsort(a) = 1;
            end
        end
    end
end

PMAll = SortingFastq3(find(PMsort),:);

PMstr1 = PMAll(:,1);
PMstr2 = PMAll(:,2);
PMstr3 = PMAll(:,3);
PMstr4 = PMAll(:,4);

[PCodondata1,Pdata1] = strmode(PMstr1);
[PCodondata2,Pdata2] = strmode(PMstr2);
[PCodondata3,Pdata3] = strmode(PMstr3);
Pdatacell = [{Pdata1} {Pdata2} {Pdata3}];
if ~isempty(DNA4)
    [PCodondata4,Pdata4] = strmode(PMstr4);
    Pdatacell = [{Pdata1} {Pdata2} {Pdata3} {Pdata4}];
end


writematrix(PCodondata1,[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"A2:D"+length(DNA1)+1)
writematrix(PCodondata2,[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"F2:I"+length(DNA2)+1)
writematrix(PCodondata3,[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"K2:N"+length(DNA3)+1)

writematrix(["A" "T" "G" "C"],[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"A1:D1")
writematrix(["A" "T" "G" "C"],[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"F1:I1")
writematrix(["A" "T" "G" "C"],[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"K1:N1")
if ~isempty(DNA4)
    writematrix(PCodondata4,[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"P2:S"+length(eval(DNA4))+1)
    writematrix(["A" "T" "G" "C"],[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"P1:S1")
end
%% file Mode & Matching comfirm
[Codondata1,data1] = strmode(str1);
[Codondata2,data2] = strmode(str2);
[Codondata3,data3] = strmode(str3);
datacell = [{data1} {data2} {data3}];

if ~isempty(DNA4)
    [Codondata4,data4] = strmode(str4);
    datacell = [{data1} {data2} {data3} {data4}];
end

csvwrite("Codon number data_"+DNA1+".csv",Codondata1)
csvwrite("Codon number data_"+DNA2+".csv",Codondata2)
csvwrite("Codon number data_"+DNA3+".csv",Codondata3)
if ~isempty(DNA4)
    csvwrite("Codon number data_"+DNA4+".csv",Codondata4)
end
%% Decode to Letter
pn1 = primercheck(Primer1{1});
pn2 = primercheck(Primer2{1});
pn3 = primercheck(Primer3{1});

decode1 = strDecode(data1,pn1);
decode2 = strDecode(data2,pn2);
decode3 = strDecode(data3,pn3);

if ~isempty(DNA4)
    pn4 = primercheck(Primer4{1});
    decode4 = strDecode(data4,pn4);
    "Decode Complete: "+decode1+" "+decode2+" "+decode3+" "+decode4
else
    "Decode Complete: "+decode1+" "+decode2+" "+decode3
end
%% functions
function [Cdata,data] = strmode(str)
    matrix = zeros(length(str),1);
    for a = 1 : length(matrix)
        matrix(a) = length(str{a});
    end
    sortmatrix = matrix(matrix ~= 0);
    [C,~,ic] = unique(sortmatrix);
    icounts = accumarray(ic,1);
    Smat = cell2mat(str(matrix==C(icounts==max(icounts))));
    data = mode(Smat,1);
    
    Cdata = zeros(size(Smat,2),4);
    for i = 1 : size(Smat,2)
        Cdata(i,1) = numel(strfind(transpose(Smat(:,i)),'A'));
        Cdata(i,2) = numel(strfind(transpose(Smat(:,i)),'T'));
        Cdata(i,3) = numel(strfind(transpose(Smat(:,i)),'G'));
        Cdata(i,4) = numel(strfind(transpose(Smat(:,i)),'C'));
    end
end

function pn = primercheck(P)
    FwPrimers = ["TTCGTTCGTCGTTGATTGGT";"AAATCCTTTGTGCCTGCCAT";"AATCATGGCCTTCAAACCGT";...
        "CTGTCCATAGCCTTGTTCGT";"TGTATTTCCTTCGGTGCTCC";"AGCCTTGTGTCCATCAATCC";...
        "GTCCAGGCAAAGATCCAGTT";"ATCCTGCAAACGCATTTCCT";"TAGCCTCCAGAATGAAACGG";...
        "TCCTCAGCCGATGAAATTCC";"TACCGCATCCTTATTCGAGC";"GAAGAGTTTAGCCACCTGGT";...
        "CAAGATTGTGGACGATTGGC";'ATGCACTCAGACGTACAACG';'GTACTGCCATCGACACTGGC';...
        'GCACTTACCTCAACATACCC'];
    pn = find(ismember(FwPrimers,P));
end

function decode = strDecode(seq,primernumber)
    payload = seq(21:end-20);
    m = 7;
    switch primernumber
        case 10
            m = 6;
        case 11
            m = 6;
        case 12
            m = 5;
    end
    binmapping = {'00','01','10','11'};
    if primernumber < 10
        [~,y] = ismember(payload,'GACT');
        chr_out = cell2mat(binmapping(y));
        binarydata = chr_out-'0';
        if rem(numel(binarydata),m) == 1 && binarydata(end) == 0
            binarydata(end) = [];
        end
        bimatrix = reshape(binarydata,[m numel(binarydata)/m])';
        dematrix = bi2de(bimatrix)';
        n = length(dematrix);
        k = n-2;
        code = rsdec(gf(dematrix,m),n,k);
        decode = char(code.x);
    else
        [~,y] = ismember(payload,'ATGC');
        chr_out = cell2mat(binmapping(y));
        binarydata = chr_out-'0';
        payldpart = binarydata(1:end-m*2);
        redunpart = binarydata(end-m*2+1:end);
        RSredundancy = bi2de(reshape(redunpart,[2 length(redunpart)/2]));
        dematrix = [payldpart RSredundancy'];
        k = length(payldpart);
        n = k+2;
        code = rsdec(gf(dematrix,m),n,k);
        code = code.x;
        cmatrix = reshape(code,[numel(code)/8 8]);
        cmatrix(:,1) = [];
        decode = char(bi2de(cmatrix))';
    end
end