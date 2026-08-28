# Clinical-only vs Omics: Paired Participant Bootstrap

ΔC-index yönü: belirtilen konfigürasyon eksi Clinical-only. Full karşılaştırmaları önceden belirlenmiş primer testlerdir; diğerleri keşifseldir.

## mixed

Kayıt: 842; katılımcı: 769; seed: 42.

| Rol | Model | Konfigürasyon | Δ C-index | %95 GA | Ham p | BH-FDR p | Yorum |
|---|---|---|---:|---|---:|---:|---|
| exploratory | coxph | RNA-seq only | -0.1563960439082709 | [-0.24334127983056134, -0.07125598026223061] | 0.001998001998001998 | 0.011988011988011988 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | Methylation only | -0.13268123030105428 | [-0.21868634435320505, -0.051663950181638954] | 0.001998001998001998 | 0.011988011988011988 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | miRNA only | -0.1791653081186828 | [-0.2669924749427176, -0.09121734076389625] | 0.001998001998001998 | 0.011988011988011988 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | Clinical + RNA-seq | -0.07294859254428865 | [-0.1442296238418519, -0.0039158654202873025] | 0.03996003996003996 | 0.11988011988011989 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | Clinical + RNA-seq + Methylation | -0.07101402021519398 | [-0.14371967723024062, 0.0004722211625071038] | 0.055944055944055944 | 0.14918414918414918 | anlamlı fark bulunamadı |
| exploratory | coxph | Clinical + RNA-seq + miRNA | -0.07842625801543313 | [-0.14767740850548122, -0.011570942878943506] | 0.02197802197802198 | 0.07535321821036108 | bootstrap %95 GA sıfırı içermiyor |
| pre-specified primary | coxph | Full | -0.0681665036409086 | [-0.13907661504282098, -0.0018592512863339966] | 0.04595404595404595 | uygulanmadı | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | Clinical + Methylation | -0.0354309314204978 | [-0.10265668102557592, 0.02706129312453069] | 0.25574425574425574 | 0.3836163836163836 | anlamlı fark bulunamadı |
| exploratory | coxph | Clinical + miRNA | -0.030996630800999903 | [-0.08123965425663078, 0.018061120677652982] | 0.1958041958041958 | 0.3132867132867133 | anlamlı fark bulunamadı |
| exploratory | rsf | RNA-seq only | -0.11941093359417454 | [-0.2041718368279563, -0.0374265590233537] | 0.007992007992007992 | 0.03836163836163836 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | rsf | Methylation only | -0.0970003260515161 | [-0.1686053339689026, -0.023875190093585793] | 0.011988011988011988 | 0.04795204795204795 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | rsf | miRNA only | -0.16497119878274102 | [-0.2543798002942152, -0.07414220380391623] | 0.001998001998001998 | 0.011988011988011988 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | rsf | Clinical + RNA-seq | -0.05522225845016848 | [-0.12484033211725189, 0.01182205247366438] | 0.1038961038961039 | 0.21178821178821178 | anlamlı fark bulunamadı |
| exploratory | rsf | Clinical + RNA-seq + Methylation | -0.05733072492120417 | [-0.12482792693390762, 0.013171308738934398] | 0.10589410589410589 | 0.21178821178821178 | anlamlı fark bulunamadı |
| exploratory | rsf | Clinical + RNA-seq + miRNA | -0.05237474187588298 | [-0.12185432066221526, 0.016243291691114426] | 0.13786213786213786 | 0.2545147160531776 | anlamlı fark bulunamadı |
| pre-specified primary | rsf | Full | -0.06089555483099662 | [-0.1296039553060451, 0.011509323077989799] | 0.08591408591408592 | uygulanmadı | anlamlı fark bulunamadı |
| exploratory | rsf | Clinical + Methylation | -0.060786870992283415 | [-0.12440968409542077, 0.006384329134261078] | 0.08191808191808192 | 0.1966033966033966 | anlamlı fark bulunamadı |
| exploratory | rsf | Clinical + miRNA | -0.03989783719160955 | [-0.10186042753141687, 0.020582969527228313] | 0.18981018981018982 | 0.3132867132867133 | anlamlı fark bulunamadı |
| exploratory | deepsurv | RNA-seq only | -0.03959352244321268 | [-0.14804579076591434, 0.06538480255839832] | 0.46153846153846156 | 0.4816053511705686 | anlamlı fark bulunamadı |
| exploratory | deepsurv | Methylation only | -0.04850559721769376 | [-0.15461950689202011, 0.05694813958767828] | 0.38961038961038963 | 0.468622286804105 | anlamlı fark bulunamadı |
| exploratory | deepsurv | miRNA only | -0.04492989892402999 | [-0.1516501800856419, 0.050370521255080764] | 0.3776223776223776 | 0.468622286804105 | anlamlı fark bulunamadı |
| exploratory | deepsurv | Clinical + RNA-seq | 0.036267796978589284 | [-0.05600517522436696, 0.1235811210760073] | 0.4115884115884116 | 0.468622286804105 | anlamlı fark bulunamadı |
| exploratory | deepsurv | Clinical + RNA-seq + Methylation | 0.03505053798500157 | [-0.06022702805032554, 0.1363993984994544] | 0.5314685314685315 | 0.5314685314685315 | anlamlı fark bulunamadı |
| exploratory | deepsurv | Clinical + RNA-seq + miRNA | 0.04707097054667975 | [-0.03690669724716238, 0.13284064094577663] | 0.30569430569430567 | 0.4315684315684315 | anlamlı fark bulunamadı |
| pre-specified primary | deepsurv | Full | -0.01613955004890777 | [-0.11995386938222488, 0.08500012221209195] | 0.7612387612387612 | uygulanmadı | anlamlı fark bulunamadı |
| exploratory | deepsurv | Clinical + Methylation | -0.0412889903271384 | [-0.13651254898435858, 0.053188290325760425] | 0.42957042957042957 | 0.468622286804105 | anlamlı fark bulunamadı |
| exploratory | deepsurv | Clinical + miRNA | 0.044375611346592736 | [-0.04684884504964585, 0.14109647655094035] | 0.4155844155844156 | 0.468622286804105 | anlamlı fark bulunamadı |

## primary

Kayıt: 762; katılımcı: 762; seed: 42.

| Rol | Model | Konfigürasyon | Δ C-index | %95 GA | Ham p | BH-FDR p | Yorum |
|---|---|---|---:|---|---:|---:|---|
| exploratory | coxph | RNA-seq only | -0.2088985266525175 | [-0.292934245759764, -0.11868188563929245] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | miRNA only | -0.2680784965088739 | [-0.341862296894609, -0.19526757941259235] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | Clinical + RNA-seq | -0.08610931684660572 | [-0.14818281605213007, -0.02389129626601684] | 0.017982017982017984 | 0.02157842157842158 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | Clinical + RNA-seq + Methylation | -0.08847774599366365 | [-0.15134199760596284, -0.02625480308160069] | 0.00999000999000999 | 0.014103543515308221 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | Clinical + RNA-seq + miRNA | -0.11446894897111737 | [-0.1787471415539899, -0.051937713511235106] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| pre-specified primary | coxph | Full | -0.12286610685614097 | [-0.19313868919883884, -0.055208255040384624] | 0.001998001998001998 | uygulanmadı | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | Clinical + Methylation | -0.025637476546399585 | [-0.05839227099157473, 0.004883377634756285] | 0.12787212787212787 | 0.13949686676959405 | anlamlı fark bulunamadı |
| exploratory | coxph | Clinical + miRNA | -0.030358955430469625 | [-0.053953817876300644, -0.007144132756775835] | 0.00999000999000999 | 0.014103543515308221 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | coxph | Methylation only | -0.14936530324400565 | [-0.26472693895904625, -0.03065237007521544] | 0.017982017982017984 | 0.02157842157842158 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | rsf | RNA-seq only | -0.2223862692627111 | [-0.3106437409576079, -0.12336145699449837] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | rsf | miRNA only | -0.23033742425640558 | [-0.299117302225224, -0.14887158247291718] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | rsf | Clinical + RNA-seq | -0.14921103626464893 | [-0.21910692720088049, -0.07104794971800378] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | rsf | Clinical + RNA-seq + Methylation | -0.1420442311832919 | [-0.21195927848296578, -0.06420880026178846] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | rsf | Clinical + RNA-seq + miRNA | -0.1440743133093415 | [-0.21322295844873773, -0.0725808242348785] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| pre-specified primary | rsf | Full | -0.14727323059887432 | [-0.21714908383911835, -0.07092354432363333] | 0.001998001998001998 | uygulanmadı | bootstrap %95 GA sıfırı içermiyor |
| exploratory | rsf | Clinical + Methylation | -0.040109501399526404 | [-0.09298041791795895, 0.012848089481378474] | 0.14185814185814186 | 0.14802588715632192 | anlamlı fark bulunamadı |
| exploratory | rsf | Clinical + miRNA | -0.01828611854449269 | [-0.0485820262612538, 0.014257227184770491] | 0.2857142857142857 | 0.2857142857142857 | anlamlı fark bulunamadı |
| exploratory | rsf | Methylation only | -0.17066290550070518 | [-0.29934988660080086, -0.020501348710712648] | 0.01998001998001998 | 0.02283430854859426 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | deepsurv | RNA-seq only | -0.20756051797853026 | [-0.31360238030379456, -0.10208453910143263] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | deepsurv | miRNA only | -0.2629110147334748 | [-0.3295932855679559, -0.18960722727931453] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | deepsurv | Clinical + RNA-seq | -0.17980068284580597 | [-0.25076520496372945, -0.1029943755658281] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | deepsurv | Clinical + RNA-seq + Methylation | -0.21020577650641314 | [-0.28278340348089626, -0.13883859697320272] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | deepsurv | Clinical + RNA-seq + miRNA | -0.24548614315154865 | [-0.3135287311159154, -0.17232574639165463] | 0.001998001998001998 | 0.003688619073234458 | bootstrap %95 GA sıfırı içermiyor |
| pre-specified primary | deepsurv | Full | -0.214850358340254 | [-0.2785581623379153, -0.15463663160779445] | 0.001998001998001998 | uygulanmadı | bootstrap %95 GA sıfırı içermiyor |
| exploratory | deepsurv | Clinical + Methylation | -0.07077604503091262 | [-0.12387887139659694, -0.01892626205400655] | 0.007992007992007992 | 0.012787212787212787 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | deepsurv | Clinical + miRNA | -0.09009258404847587 | [-0.1595354407364765, -0.01770621578215301] | 0.017982017982017984 | 0.02157842157842158 | bootstrap %95 GA sıfırı içermiyor |
| exploratory | deepsurv | Methylation only | -0.17950164551010817 | [-0.275558226569973, -0.07587852884651938] | 0.005994005994005994 | 0.010275438846867418 | bootstrap %95 GA sıfırı içermiyor |
