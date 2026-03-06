# SecretSharing_Powerof2
Secret sharing implementation over GF(2^8)

# 說明
本程式為將Galois Field GF(2^8^)應用於secret sharing之測試實作，可用於測試並比較GF(2^k^)與GF(p)查表與否、以及不同解密法之效率差異，功能較為簡陋，僅為測試之用。

## define
**TEST_MODE**: 若為true則不顯示解密結果，用於計算花費時間。  
**TYPE**: 本次測試選用的四則運算方式，即GF(2^k^)或GF(p)。  
**TYPE_STR**: 同上，用於將測試結果輸出的字串上。  
**USE_TABLE**: 是否於乘、除法中使用查表法。  
**TEST_TIME**: 本次測試運算次數，運算越多次時間越準確。  
**PRIME_MAX**: GF(p)中p的數值，需為質數。  
**GF_LENGTH**: GF(2^k^)中k的數值，  
**GF_MAX_EXP**: GF(2^k^)中2^k^計算後的數值減一。  
**GF_EXP**: GF(2^k^)中2^k^的十六進位表示法。  
**GF_IRPOLY**: GF(2^k^)之irreducible polynomial。  

## 其他參數
如同範例程式，需設定**人數**、**share數**、**secret**，拉格朗日插值法的secret為一int陣列。
