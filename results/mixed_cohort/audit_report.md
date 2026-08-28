# Faz 1 Pipeline Denetimi

Bu rapor mevcut scriptler, kaynak TSV başlıkları ve üretilmiş Parquet dosyaları doğrudan incelenerek hazırlanmıştır. Aşağıdaki sayılar çalıştırılan koddan alınmıştır.

## 1. Sağkalım sonlanım noktası

- `TCGA-BRCA.survival.tsv.gz` gerçek kolonları: `sample`, `OS.time`, `OS`, `_PATIENT`.
- `prepare_dataset.py`, event için `OS`, zaman için `OS.time` kullanmıştır. Bu sonlanım noktası Overall Survival'dır.
- Zaman değerleri kaynak dosyada ve TCGA/Xena survival tanımında gün cinsindedir; pipeline herhangi bir birim dönüşümü yapmamıştır.
- Eksik etiketli 11 örnek çıkarıldıktan sonra analiz matrisinde 842 örnek vardır: event=1 olan 132, censored (event=0) olan 710 örnek.
- Censored örneklerde medyan takip süresi: 809.0 gün.
- Event olan örneklerde medyan gözlenen sağkalım zamanı: 1223.0 gün.
- Bu iki medyan, Kaplan–Meier medyan sağkalımı değildir; istenen alt gruplardaki ham `OS.time` medyanlarıdır.

## 2. Hasta/örnek eşleştirme

- `barcode()` ilk 15 karakteri kullanır. Bu anahtar `TCGA-XX-XXXX-SS` biçimindedir: ilk 12 karakter participant, son iki karakter sample type kodudur. Vial harfi (örn. `A`) 16. karakterdir ve anahtara dahil değildir.
- Aynı 15-karakter anahtara ait birden fazla vial/aliquot varsa `prepare_dataset.py` omik değerlerini sütun bazında ortalama ile birleştirmiştir. Metadata tekrarlarında ilk dolu değer korunmuş ve çatışmalar sayılmıştır.
- Participant düzeyinde tek kayda indirgeme yapılmamıştır. Nihai 842 satır 769 benzersiz participant içerir; 71 participant birden fazla sample key ile temsil edilir ve toplam 73 ek satır vardır.
- Nihai sample-type dağılımı: `01` (primary solid tumor) 762, `11` (solid tissue normal) 75, `06` (metastatic) 5.
- Normal doku dışlanmamıştır. Mevcut `train_test_split` sample key düzeyindedir; bu nedenle “hasta bazlı split” iddiası doğru değildir ve aynı participant'ın farklı sample type'ları farklı setlere/fold'lara gidebilir.

## 3. Ön işleme ve ölçeklendirme

- RNA STAR-TPM değerlerine `log2(TPM+1)` veya başka bir dönüşüm uygulanmamıştır; kaynak sayısal değerler doğrudan `float32` olarak taşınmıştır.
- Methylation değerleri M-value'ya çevrilmemiş, kaynak beta değerleri olarak bırakılmıştır.
- miRNA değerlerine ek log dönüşümü veya normalizasyon uygulanmamıştır.
- `impute_and_split.py`: methylation top-10.000 varyans seçimi yalnız ilk train split üzerinde öğrenilmiş; `KNNImputer(k=5)` yalnız train'de fit edilip test'e transform uygulanmıştır.
- Clinical yaş için median imputasyon, kategorikler için mode imputasyon ve one-hot encoder yalnız ilk train split'te fit edilmiştir. Bu aşamada yaş standardize edilmemiştir.
- CoxPH'de `StandardScaler` train üzerinde fit edilip test'e uygulanmıştır. DeepSurv'de ortalama/standart sapma yalnız iç fit alt kümesinden öğrenilmiştir. RSF ölçeklendirilmemiştir.
- Genel sıfır-varyans filtresi yoktur. CoxPH öncesinde event ve censored alt gruplarından herhangi birinde varyansı `<=1e-12` olan kolonlar complete-separation ön kontrolü olarak çıkarılmıştır. RSF/DeepSurv için aynı filtre uygulanmamıştır.

## 4. Methylation varyans filtresi ve CV leakage

- `cross_validate.py`, `X_train.parquet` ve `X_test.parquet` dosyalarını birleştirir. Bu iki dosyadaki 10.000 methylation probu, daha önce `impute_and_split.py` tarafından yalnız orijinal 673 örneklik train split varyansına göre seçilmiştir.
- Outer 5-fold CV içinde univariate Cox seçimi yeniden yapılmaktadır; ancak 486.427→10.000 methylation varyans ön filtresi her fold'un training kısmında yeniden hesaplanmamaktadır.
- Bu nedenle mevcut `cv_results.md` tamamen fold-local bir pipeline sonucu değildir. Özellikle orijinal test örnekleri outer training'e girdiğinde, sabit prob listesi başka bir örnek kümesinden öğrenilmiş kalır. Clinical encoder ve methylation KNN imputasyonu da orijinal split üzerinde fit edilmiş temsillerden taşınmıştır.
- Sonuç: leakage/metodolojik uyumsuzluk bulundu. Yeni Faz-1 CV scriptleri ham methylation ve clinical girdilerinden outer-fold training-only fit yapmalıdır; eski CV sayıları düzeltilmiş sonuç olarak kullanılmamalıdır.

## 5. Yazılım sürümleri

Mevcut `.venv` eksik/bozuk olduğu için yeniden üretim amacıyla kısa yollu `C:\moenv` ortamı oluşturulmuş ve gereksinimler kurulmuştur. Import edilebilir analiz ortamındaki sürümler:

- pandas: 2.3.3
- numpy: 2.5.2
- scikit-learn: 1.9.0
- scikit-survival: 0.28.0
- lifelines: 0.30.3
- torch: 2.13.0+cpu
- shap: 0.52.0
- statsmodels: 0.15.0

Windows uzun-yol sınırı nedeniyle proje içindeki eski `.venv` PyTorch kurulumu tamamlanamamıştı; yeni sonuçlar üretilirken kullanılan Python yürütücüsü açıkça loglanmalıdır.
