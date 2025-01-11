%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Letter Encoding %%%%%%%%%%%%%%%%%%%%%%%%%%
% File Selection
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
% The encoding technique for words that are nouns and other cases of speech is different.
% Two redundancy are added to the RS code, noun words can recover one ascii code, and verbs and adverb words can recover one binary data.
%%
% 다른 단어를 쓰고 싶으면 "Text for DNA Encoding"에 원하는 단어 입력
% 프라이머를 바꾸고 싶으면 primer select.m 켜서 primer 바꾸거나 추가
%%
clear % 작업공간 모든 변수 제거
clc % 명령창 깨끗하게
scan = input('Text for DNA Encoding: ','s'); % input: 명령창 질문? 정의? 물어보는거? 내가 직접 쓰는것/ output이 내가 명령창에 써놓은 것/ 's'는 문자형태로 변수에 저장
j = input('File select(number) : ');
strform = input('is it a noun?(Y/N): ','s');
int_scan = double(scan); % 문자형태를 ascii코드표로 숫자 변환하는것

bindata = [zeros(length(int_scan),1) de2bi(int_scan)]; % de2bi: 십진수를 이진수로 변환하는 것

if strcmp(strform,'Y') % 만약에 strcmp(문자형태 두개를 넣었을때 서로 일치하는지 확인)/ true(1) false(0)
    k = length(int_scan); % k는 데이터의 크기
    n = k+2;
    m = length(de2bi(max(int_scan))); % length: 행렬 중 더 큰 것의 길이/ max: 최대값/ 괄호 가장 안쪽부터 풀기
    payload = int_scan;
else
    k = numel(bindata); % numel: 데이터의 총개수
    n = k+2;
    m = size(de2bi(n),2); % size: 첫번째 변수(행렬)의 행 or 열을 출력 
    payload = reshape(bindata,[1 numel(bindata)]); % 행렬 변환
end
gf_data = gf(payload,m); % galois field 형태로 만드는 것(RS code를 위한 것)
code = rsenc(gf_data,n,k); % RS 인코딩/ 에러 고치는 코드, 데이터가 뒤에 중복되어서 붙음
% code라는 변수까지 실행하면 데이터 끝에 RS 중복이 붙는다.

if strcmp(strform,'Y')
    bin_code = de2bi(code.x); % RS 코드 붙은 데이터를 이진수로 변환
    payload = reshape(bin_code',[1 numel(bin_code)]); % 행렬 변환
    if rem(numel(payload),2) ~=0 % rem: 나머지 출력/ 데이터 개수를 짝수로 만들려고
        payload(end+1) = 0;
    end
    bases = 'GACT';
else
    bin_code = code.x;
    parity_file = de2bi(double(bin_code(k+1:end))); % 십진수를 이진수로 변환할때 큰 수 기준으로 행렬이 만들어진다.(작은 수는 0이 뒤에 붙는다)
    parity_bin = reshape(parity_file,[1 numel(parity_file)]);
    payload = [bin_code(1:k) parity_bin];
    bases = 'ATGC';
end

log_inf = logical(payload); % logical 형태
result = repmat(' ',1,length(log_inf)/2); % repmat: 반복하여 행렬을 생성한다.
for e = 1 : 2 : length(log_inf) % for: e를 순서대로 넣으면서 루프한다.
    module = 2 * log_inf(e) + log_inf(e+1) + 1;
    result((e+1)/2) = bases(module);
end

[Fw,Rv] = primerselect(j); % 파일번호에 따라 프라이머 정하는 함수
EncSequence = [Fw result Rv]; % []: 행렬 합치기
writematrix(EncSequence,"DNA_encoded_"+scan+".txt")
EncSequence