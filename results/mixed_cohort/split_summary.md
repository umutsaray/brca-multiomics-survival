# Leakage-Safe Imputation ve Train/Test Split Özeti

## Örnekler ve etiketler

- Başlangıç örnek sayısı: 853
- Eksik survival kontrol bayrağıyla çıkarılan örnek: 11
- Kalan örnek sayısı: 842
- Train örnek sayısı: 673
- Test örnek sayısı: 169
- Genel event oranı: 15.676960%
- Train event oranı: 15.750371%
- Test event oranı: 15.384615%

## Özellikler ve leakage kontrolleri

- RNA özelliği: 60660; doğrulanan NaN sayısı: 0
- miRNA özelliği: 1881; doğrulanan NaN sayısı: 0
- Başlangıç methylation probe sayısı: 486427
- Yalnız train varyansıyla seçilen methylation probe sayısı: 10000
- Varyans ön filtresinde elenen probe sayısı: 476427
- Varyans tarama süresi: 65.17 saniye
- KNNImputer komşu sayısı: 5
- KNN fit(train) + transform(train/test) süresi: 24.74 saniye
- Methylation NaN train: 726559 -> 0
- Methylation NaN test: 191480 -> 0
- Clinical kovaryatlar: ['age_at_index.demographic', 'gender.demographic', 'race.demographic', 'ethnicity.demographic', 'ajcc_pathologic_stage.diagnoses']
- Clinical median/mode imputation ve one-hot encoder yalnız train setinde fit edildi.
- Train/test ayrımı, varyans filtresi ve tüm imputer fit işlemlerinden önce yapıldı.

## Nihai matrisler

- X_train: 673 x 72564
- X_test: 169 x 72564
- y_train: 673 x 2
- y_test: 169 x 2
- test_size=0.2, random_state=42, stratify=event
