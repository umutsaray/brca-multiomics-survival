# Phase 17 Task 7 — Uno versus Harrell C-index sensitivity

**Exploratory sensitivity analysis.** No model was refit; only stored out-of-fold predictions were used.

Design fixed before results: `{"stored_predictions_only": true, "participant_bootstrap_key": "first 12 characters of sample_key", "cohorts": ["mixed", "primary"], "configurations": ["Clinical-only", "Full"], "models": ["coxph", "rsf", "deepsurv"], "taus_days": {"existing_none": null, "3_year": 1095.0, "5_year": 1825.0, "7_year": 2555.0}, "weight_constructions": ["other_four_oof", "outer_training_outcomes"], "delta": "Clinical-only minus Full", "bootstrap_iterations": 9999, "bootstrap_seed": 42}`

Analysis runtime seconds: `8039.219540800201`.

## Per-fold Harrell and Uno C-index

The Uno columns below use the existing censoring-weight construction (outcomes from the other four OOF folds). Values for the independently assembled outer-training-outcome construction are retained in the JSON.

| Cohort | Fold | Configuration | Model | Harrell | Uno existing tau | Uno 3-year | Uno 5-year | Uno 7-year |
|---|---:|---|---|---:|---:|---:|---:|---:|
| mixed | 1 | Clinical-only | coxph | 0.667704280155642 | 0.4169198217081488 | 0.7727242645008323 | 0.6126151101975731 | 0.49094164597344264 |
| mixed | 1 | Clinical-only | rsf | 0.6673151750972762 | 0.5219792253933844 | 0.7208966381589219 | 0.6067117077662906 | 0.5446042108960689 |
| mixed | 1 | Clinical-only | deepsurv | 0.6832684824902724 | 0.527618565592039 | 0.7116275334523979 | 0.6338110874597486 | 0.5562277469553951 |
| mixed | 1 | Full | coxph | 0.621011673151751 | 0.5686073716343351 | 0.7187908047023986 | 0.5424056243648877 | 0.5712793739254963 |
| mixed | 1 | Full | rsf | 0.5587548638132296 | 0.6603013440178999 | 0.5345180517568557 | 0.5073859427002545 | 0.644957932341626 |
| mixed | 1 | Full | deepsurv | 0.5789883268482491 | 0.66814535602882 | 0.612491051708731 | 0.5110201019740617 | 0.5958174872471266 |
| mixed | 2 | Clinical-only | coxph | 0.7494632889652211 | 0.6873610994411892 | 0.736460924806234 | 0.7014657289206586 | 0.6577669911696629 |
| mixed | 2 | Clinical-only | rsf | 0.6983683984542722 | 0.6284165774162562 | 0.6671737942693284 | 0.6821500207198125 | 0.5805017746335573 |
| mixed | 2 | Clinical-only | deepsurv | 0.7052382996994419 | 0.6715843942432957 | 0.6911191851129213 | 0.6825187772557857 | 0.6349189915936274 |
| mixed | 2 | Full | coxph | 0.7548303993130099 | 0.7559099485359396 | 0.7997840615632792 | 0.786918060974005 | 0.7835656519744277 |
| mixed | 2 | Full | rsf | 0.7805925289823958 | 0.7183149068012263 | 0.808514077506858 | 0.7891753816241123 | 0.7812352059232682 |
| mixed | 2 | Full | deepsurv | 0.69944182052383 | 0.6836228555440127 | 0.7137076239494895 | 0.6964278046913708 | 0.6991224859432263 |
| mixed | 3 | Clinical-only | coxph | 0.7144927536231884 | 0.6214474665897076 | 0.7081876782486777 | 0.7729619927565541 | 0.7150789495944204 |
| mixed | 3 | Clinical-only | rsf | 0.7840579710144927 | 0.5968958336782427 | 0.7878913098745682 | 0.8260202773276369 | 0.7318606689434779 |
| mixed | 3 | Clinical-only | deepsurv | 0.736231884057971 | 0.4381350411816114 | 0.7274745605125069 | 0.7769033721863197 | 0.6748283126015268 |
| mixed | 3 | Full | coxph | 0.6601449275362319 | 0.5494104539656565 | 0.7085583546065217 | 0.6508128632774463 | 0.5689765889506072 |
| mixed | 3 | Full | rsf | 0.6086956521739131 | 0.5202553872345652 | 0.6085038171778281 | 0.5932706782412787 | 0.5917491328149632 |
| mixed | 3 | Full | deepsurv | 0.4326086956521739 | 0.6053122069296103 | 0.4551031692544583 | 0.4583175547277688 | 0.5017138406340865 |
| mixed | 4 | Clinical-only | coxph | 0.792255027813436 | 0.774843517232488 | 0.7673951445466618 | 0.7254738917119418 | 0.7661404938840143 |
| mixed | 4 | Clinical-only | rsf | 0.7849807445442876 | 0.7449839363366868 | 0.772996689158941 | 0.7037919668122642 | 0.725950827017395 |
| mixed | 4 | Clinical-only | deepsurv | 0.7811296534017972 | 0.7165187058437184 | 0.7319316813786589 | 0.7173830065540868 | 0.7568289919383007 |
| mixed | 4 | Full | coxph | 0.7145913564398801 | 0.7199623037888033 | 0.7381924489743451 | 0.6969411181121805 | 0.6922728232949851 |
| mixed | 4 | Full | rsf | 0.7163029525032092 | 0.6844722499030806 | 0.7598059475682188 | 0.6963759159631571 | 0.6791003394936781 |
| mixed | 4 | Full | deepsurv | 0.6302952503209243 | 0.6396103293263061 | 0.6876643085067426 | 0.6224643336142581 | 0.6187181911551247 |
| mixed | 5 | Clinical-only | coxph | 0.7548906789413119 | 0.6380407912345885 | 0.6339573928647696 | 0.6849897073088664 | 0.7437861632847014 |
| mixed | 5 | Clinical-only | rsf | 0.7589182968929804 | 0.6542042108429301 | 0.7011450694050659 | 0.6539298988795842 | 0.718523646884606 |
| mixed | 5 | Clinical-only | deepsurv | 0.7917146144994246 | 0.6542371999798763 | 0.7507193679059172 | 0.7062873303780965 | 0.7611085995545347 |
| mixed | 5 | Full | coxph | 0.6006904487917146 | 0.6414912217614086 | 0.6096907503167855 | 0.49224408831153776 | 0.5545556316451821 |
| mixed | 5 | Full | rsf | 0.6737629459148446 | 0.6978204565837756 | 0.7007768126406276 | 0.665537849364232 | 0.7117346059491384 |
| mixed | 5 | Full | deepsurv | 0.619102416570771 | 0.6603284400919418 | 0.6589859764546809 | 0.6393023252753407 | 0.6498199860360961 |
| primary | 1 | Clinical-only | coxph | 0.8271312584573748 | 0.28342436109989066 | 0.8218429127059051 | 0.7858156089286885 | 0.8056826678975486 |
| primary | 1 | Clinical-only | rsf | 0.804127198917456 | 0.28176162257952453 | 0.7712584696848617 | 0.7710411443439747 | 0.8001639450146913 |
| primary | 1 | Clinical-only | deepsurv | 0.8298376184032477 | 0.27698732727417624 | 0.7985671317664206 | 0.7854345613277577 | 0.7939619871906759 |
| primary | 1 | Full | coxph | 0.7401894451962111 | 0.26611520415564355 | 0.7026017595309348 | 0.6637704930351553 | 0.7385884027920022 |
| primary | 1 | Full | rsf | 0.6427604871447903 | 0.2516392892160576 | 0.6480606441393125 | 0.6618714388496545 | 0.7098304199277305 |
| primary | 1 | Full | deepsurv | 0.5690121786197564 | 0.21213073487632986 | 0.6011380579147628 | 0.5953931814428437 | 0.646208682908189 |
| primary | 2 | Clinical-only | coxph | 0.6362637362637362 | 0.40240012339219766 | 0.7516852791355434 | 0.6749967022431779 | 0.6016719453409748 |
| primary | 2 | Clinical-only | rsf | 0.6532967032967033 | 0.4669887392384994 | 0.7368591028683364 | 0.6830197318422061 | 0.6426663581284108 |
| primary | 2 | Clinical-only | deepsurv | 0.589010989010989 | 0.42227608496174557 | 0.6663545439949442 | 0.5917519557591624 | 0.6108009129483146 |
| primary | 2 | Full | coxph | 0.6615384615384615 | 0.6425247016023758 | 0.7067889290057632 | 0.6802422713601716 | 0.5749558961436029 |
| primary | 2 | Full | rsf | 0.6538461538461539 | 0.6359154143801716 | 0.7989745852383208 | 0.45210061182562666 | 0.535334194759015 |
| primary | 2 | Full | deepsurv | 0.4857142857142857 | 0.5405963770880057 | 0.32128359002993195 | 0.5294854732887672 | 0.5397395303955169 |
| primary | 3 | Clinical-only | coxph | 0.7958096590909091 | 0.6948094388296805 | 0.7813986133195527 | 0.7893868273589936 | 0.7284461934063837 |
| primary | 3 | Clinical-only | rsf | 0.8178267045454546 | 0.6555201483362569 | 0.7736072018436955 | 0.799716756487499 | 0.714207898674373 |
| primary | 3 | Clinical-only | deepsurv | 0.78515625 | 0.6882873548641986 | 0.759173610291675 | 0.7805193416327241 | 0.7204726665095298 |
| primary | 3 | Full | coxph | 0.6200284090909091 | 0.5558640004188823 | 0.656945518251227 | 0.6568635493392054 | 0.5440471876438906 |
| primary | 3 | Full | rsf | 0.5987215909090909 | 0.5686335125643109 | 0.6864561155165033 | 0.7000547963334667 | 0.6079852663115627 |
| primary | 3 | Full | deepsurv | 0.5426136363636364 | 0.5288461277541388 | 0.5794685336199322 | 0.5672799834105229 | 0.5255489397439963 |
| primary | 4 | Clinical-only | coxph | 0.7361495844875346 | 0.6527986912170504 | 0.7015229012123435 | 0.6929530116567465 | 0.6527986912170504 |
| primary | 4 | Clinical-only | rsf | 0.7347645429362881 | 0.6681417660080842 | 0.7168413306341098 | 0.665088035901347 | 0.6681417660080842 |
| primary | 4 | Clinical-only | deepsurv | 0.7063711911357341 | 0.6051282692920242 | 0.6714888633076184 | 0.6507302890167782 | 0.6051282692920242 |
| primary | 4 | Full | coxph | 0.6814404432132964 | 0.6676603820259045 | 0.6918280009581653 | 0.6645653151176728 | 0.6676603820259045 |
| primary | 4 | Full | rsf | 0.6565096952908587 | 0.6527137528412913 | 0.5697275817498891 | 0.6356136157916205 | 0.6527137528412913 |
| primary | 4 | Full | deepsurv | 0.5457063711911357 | 0.5561019589363643 | 0.4959473172100652 | 0.5720508689375748 | 0.5561019589363643 |
| primary | 5 | Clinical-only | coxph | 0.80375 | 0.623892170691066 | 0.7666906858613757 | 0.6802624027617769 | 0.7350273788668016 |
| primary | 5 | Clinical-only | rsf | 0.77625 | 0.5921226179320089 | 0.7794735760113052 | 0.6812948205396117 | 0.7073461340739867 |
| primary | 5 | Clinical-only | deepsurv | 0.81125 | 0.6147962786133246 | 0.8369335160447932 | 0.6897975130541698 | 0.7324385768149145 |
| primary | 5 | Full | coxph | 0.5591666666666667 | 0.44213383897378833 | 0.5319135606195526 | 0.4330289280620638 | 0.5414576620121844 |
| primary | 5 | Full | rsf | 0.5583333333333333 | 0.4594037325950737 | 0.5241378643322615 | 0.5770066338598514 | 0.6315392737281241 |
| primary | 5 | Full | deepsurv | 0.48833333333333334 | 0.3938444715823263 | 0.48366438288793717 | 0.5863283487661259 | 0.6002679963447767 |

## IPCW diagnostics and participants at risk

`Maximum IPCW` is the maximum base event weight 1/G; `maximum squared Uno weight` is the corresponding weight used by Uno concordance. ESS is computed from the squared Uno weights. Participant counts use the 12-character TCGA barcode.

| Cohort | Fold | Tau | Weight construction | Participants at risk | Rows at risk | Maximum IPCW | Maximum squared Uno weight | Weight ESS | Error |
|---|---:|---|---|---:|---:|---:|---:|---:|---|
| mixed | 1 | existing_none | other_four_oof | None | None | 68.41985879361502 | 4681.277077338219 | 1.451563854349597 | None |
| mixed | 1 | 3_year | other_four_oof | 68 | 79 | 2.1507234935818165 | 4.625611545844774 | 4.966236682035535 | None |
| mixed | 1 | 5_year | other_four_oof | 35 | 41 | 3.268671396260993 | 10.684212696734788 | 7.985964787788984 | None |
| mixed | 1 | 7_year | other_four_oof | 21 | 23 | 5.041926929967179 | 25.421027167128262 | 14.12781602155068 | None |
| mixed | 1 | existing_none | outer_training_outcomes | None | None | 68.41985879361502 | 4681.277077338219 | 1.451563854349597 | None |
| mixed | 1 | 3_year | outer_training_outcomes | 68 | 79 | 2.1507234935818165 | 4.625611545844774 | 4.966236682035535 | None |
| mixed | 1 | 5_year | outer_training_outcomes | 35 | 41 | 3.268671396260993 | 10.684212696734788 | 7.985964787788984 | None |
| mixed | 1 | 7_year | outer_training_outcomes | 21 | 23 | 5.041926929967179 | 25.421027167128262 | 14.12781602155068 | None |
| mixed | 2 | existing_none | other_four_oof | None | None | 6.075954130514331 | 36.91721859611416 | 10.256602486644764 | None |
| mixed | 2 | 3_year | other_four_oof | 62 | 72 | 2.13938166909806 | 4.5769539260728 | 12.77874985336506 | None |
| mixed | 2 | 5_year | other_four_oof | 39 | 43 | 3.342947715048321 | 11.175299425546791 | 12.52378441291785 | None |
| mixed | 2 | 7_year | other_four_oof | 23 | 25 | 5.239291997569176 | 27.45018063579241 | 9.603886729330682 | None |
| mixed | 2 | existing_none | outer_training_outcomes | None | None | 6.075954130514331 | 36.91721859611416 | 10.256602486644764 | None |
| mixed | 2 | 3_year | outer_training_outcomes | 62 | 72 | 2.13938166909806 | 4.5769539260728 | 12.77874985336506 | None |
| mixed | 2 | 5_year | outer_training_outcomes | 39 | 43 | 3.342947715048321 | 11.175299425546791 | 12.52378441291785 | None |
| mixed | 2 | 7_year | outer_training_outcomes | 23 | 25 | 5.239291997569176 | 27.45018063579241 | 9.603886729330682 | None |
| mixed | 3 | existing_none | other_four_oof | None | None | 19.188522802335733 | 368.19940733575834 | 5.763742495216715 | None |
| mixed | 3 | 3_year | other_four_oof | 54 | 61 | 1.8669474812347373 | 3.48549289768873 | 10.107719712097365 | None |
| mixed | 3 | 5_year | other_four_oof | 32 | 34 | 3.2244580561649148 | 10.397129755966821 | 11.783625506218034 | None |
| mixed | 3 | 7_year | other_four_oof | 19 | 21 | 5.172795575623766 | 26.75781406719281 | 9.747722998295275 | None |
| mixed | 3 | existing_none | outer_training_outcomes | None | None | 19.188522802335733 | 368.19940733575834 | 5.763742495216715 | None |
| mixed | 3 | 3_year | outer_training_outcomes | 54 | 61 | 1.8669474812347373 | 3.48549289768873 | 10.107719712097365 | None |
| mixed | 3 | 5_year | outer_training_outcomes | 32 | 34 | 3.2244580561649148 | 10.397129755966821 | 11.783625506218034 | None |
| mixed | 3 | 7_year | outer_training_outcomes | 19 | 21 | 5.172795575623766 | 26.75781406719281 | 9.747722998295275 | None |
| mixed | 4 | existing_none | other_four_oof | None | None | 16.99106994263566 | 288.696457795537 | 3.013512735912174 | None |
| mixed | 4 | 3_year | other_four_oof | 67 | 75 | 2.1319759334512445 | 4.545321380815305 | 15.183510994966188 | None |
| mixed | 4 | 5_year | other_four_oof | 40 | 47 | 3.4953454625489924 | 12.217439902561829 | 14.800497899468091 | None |
| mixed | 4 | 7_year | other_four_oof | 18 | 19 | 4.799846100728245 | 23.038522590676138 | 13.02366934401516 | None |
| mixed | 4 | existing_none | outer_training_outcomes | None | None | 16.99106994263566 | 288.696457795537 | 3.013512735912174 | None |
| mixed | 4 | 3_year | outer_training_outcomes | 67 | 75 | 2.1319759334512445 | 4.545321380815305 | 15.183510994966188 | None |
| mixed | 4 | 5_year | outer_training_outcomes | 40 | 47 | 3.4953454625489924 | 12.217439902561829 | 14.800497899468091 | None |
| mixed | 4 | 7_year | outer_training_outcomes | 18 | 19 | 4.799846100728245 | 23.038522590676138 | 13.02366934401516 | None |
| mixed | 5 | existing_none | other_four_oof | None | None | 10.029883609069955 | 100.59856521149015 | 9.006551566348968 | None |
| mixed | 5 | 3_year | other_four_oof | 58 | 65 | 1.9108443192563394 | 3.651326012434223 | 8.154043908592287 | None |
| mixed | 5 | 5_year | other_four_oof | 30 | 34 | 3.1599277721021597 | 9.985143524902519 | 14.285471011418787 | None |
| mixed | 5 | 7_year | other_four_oof | 21 | 22 | 5.180130025894344 | 26.833747085172135 | 10.586145193347702 | None |
| mixed | 5 | existing_none | outer_training_outcomes | None | None | 10.029883609069955 | 100.59856521149015 | 9.006551566348968 | None |
| mixed | 5 | 3_year | outer_training_outcomes | 58 | 65 | 1.9108443192563394 | 3.651326012434223 | 8.154043908592287 | None |
| mixed | 5 | 5_year | outer_training_outcomes | 30 | 34 | 3.1599277721021597 | 9.985143524902519 | 14.285471011418787 | None |
| mixed | 5 | 7_year | outer_training_outcomes | 21 | 22 | 5.180130025894344 | 26.833747085172135 | 10.586145193347702 | None |
| primary | 1 | existing_none | other_four_oof | None | None | 100.6695433324238 | 10134.356954758752 | 1.0315259503709109 | None |
| primary | 1 | 3_year | other_four_oof | 61 | 61 | 1.9640306975253479 | 3.8574165808219045 | 10.113632740897158 | None |
| primary | 1 | 5_year | other_four_oof | 30 | 30 | 3.3912039863481542 | 11.500264477023613 | 8.933746957813039 | None |
| primary | 1 | 7_year | other_four_oof | 19 | 19 | 5.285829742479845 | 27.939996066484547 | 7.0190643031954005 | None |
| primary | 1 | existing_none | outer_training_outcomes | None | None | 100.6695433324238 | 10134.356954758752 | 1.0315259503709109 | None |
| primary | 1 | 3_year | outer_training_outcomes | 61 | 61 | 1.9640306975253479 | 3.8574165808219045 | 10.113632740897158 | None |
| primary | 1 | 5_year | outer_training_outcomes | 30 | 30 | 3.3912039863481542 | 11.500264477023613 | 8.933746957813039 | None |
| primary | 1 | 7_year | outer_training_outcomes | 19 | 19 | 5.285829742479845 | 27.939996066484547 | 7.0190643031954005 | None |
| primary | 2 | existing_none | other_four_oof | None | None | 20.328121610833378 | 413.232528224831 | 5.386844636015611 | None |
| primary | 2 | 3_year | other_four_oof | 62 | 62 | 1.8311088286576387 | 3.3529595423879495 | 3.233470816426145 | None |
| primary | 2 | 5_year | other_four_oof | 43 | 43 | 3.5626152114240264 | 12.692227144669861 | 5.076121072408272 | None |
| primary | 2 | 7_year | other_four_oof | 24 | 24 | 5.72173675251308 | 32.73827146505893 | 7.856265040264906 | None |
| primary | 2 | existing_none | outer_training_outcomes | None | None | 20.328121610833378 | 413.232528224831 | 5.386844636015611 | None |
| primary | 2 | 3_year | outer_training_outcomes | 62 | 62 | 1.8311088286576387 | 3.3529595423879495 | 3.233470816426145 | None |
| primary | 2 | 5_year | outer_training_outcomes | 43 | 43 | 3.5626152114240264 | 12.692227144669861 | 5.076121072408272 | None |
| primary | 2 | 7_year | outer_training_outcomes | 24 | 24 | 5.72173675251308 | 32.73827146505893 | 7.856265040264906 | None |
| primary | 3 | existing_none | other_four_oof | None | None | 10.583694668845325 | 112.01459284334496 | 4.337205072009526 | None |
| primary | 3 | 3_year | other_four_oof | 54 | 54 | 2.1881567910698667 | 4.7880301423051765 | 9.87544398956631 | None |
| primary | 3 | 5_year | other_four_oof | 34 | 34 | 2.5521850364389844 | 6.51364846022306 | 11.10096843321456 | None |
| primary | 3 | 7_year | other_four_oof | 20 | 20 | 5.272611124203682 | 27.80042806707642 | 7.433564726572515 | None |
| primary | 3 | existing_none | outer_training_outcomes | None | None | 10.583694668845325 | 112.01459284334496 | 4.337205072009526 | None |
| primary | 3 | 3_year | outer_training_outcomes | 54 | 54 | 2.1881567910698667 | 4.7880301423051765 | 9.87544398956631 | None |
| primary | 3 | 5_year | outer_training_outcomes | 34 | 34 | 2.5521850364389844 | 6.51364846022306 | 11.10096843321456 | None |
| primary | 3 | 7_year | outer_training_outcomes | 20 | 20 | 5.272611124203682 | 27.80042806707642 | 7.433564726572515 | None |
| primary | 4 | existing_none | other_four_oof | None | None | 4.491375991130539 | 20.172458293703833 | 12.510896414433068 | None |
| primary | 4 | 3_year | other_four_oof | 67 | 67 | 2.3007844617936453 | 5.293609139631074 | 9.37956477223749 | None |
| primary | 4 | 5_year | other_four_oof | 37 | 37 | 3.6122454185868733 | 13.048316964101856 | 13.0806180945505 | None |
| primary | 4 | 7_year | other_four_oof | 17 | 17 | 4.491375991130539 | 20.172458293703833 | 12.510896414433068 | None |
| primary | 4 | existing_none | outer_training_outcomes | None | None | 4.491375991130539 | 20.172458293703833 | 12.510896414433068 | None |
| primary | 4 | 3_year | outer_training_outcomes | 67 | 67 | 2.3007844617936453 | 5.293609139631074 | 9.37956477223749 | None |
| primary | 4 | 5_year | outer_training_outcomes | 37 | 37 | 3.6122454185868733 | 13.048316964101856 | 13.0806180945505 | None |
| primary | 4 | 7_year | outer_training_outcomes | 17 | 17 | 4.491375991130539 | 20.172458293703833 | 12.510896414433068 | None |
| primary | 5 | existing_none | other_four_oof | None | None | 22.37789990070304 | 500.77040396588507 | 2.8246229690430913 | None |
| primary | 5 | 3_year | other_four_oof | 63 | 63 | 2.1434611161602595 | 4.594425556490985 | 7.318516832739689 | None |
| primary | 5 | 5_year | other_four_oof | 32 | 32 | 3.4652350413784547 | 12.00785389199714 | 9.210670782786876 | None |
| primary | 5 | 7_year | other_four_oof | 22 | 22 | 4.930525279208112 | 24.31007952891023 | 8.938210380936038 | None |
| primary | 5 | existing_none | outer_training_outcomes | None | None | 22.37789990070304 | 500.77040396588507 | 2.8246229690430913 | None |
| primary | 5 | 3_year | outer_training_outcomes | 63 | 63 | 2.1434611161602595 | 4.594425556490985 | 7.318516832739689 | None |
| primary | 5 | 5_year | outer_training_outcomes | 32 | 32 | 3.4652350413784547 | 12.00785389199714 | 9.210670782786876 | None |
| primary | 5 | 7_year | outer_training_outcomes | 22 | 22 | 4.930525279208112 | 24.31007952891023 | 8.938210380936038 | None |

## Clinical-only minus Full Uno differences

The estimand is the unweighted mean of five within-fold C-index differences. Positive values favor Clinical-only. The bootstrap resamples 12-character participants independently within each fold. A valid-iteration count below 9,999 means that at least one resampled fold had no estimable comparable pair; such draws are stored as `null`, not silently replaced. Ordering change compares each tau with the existing-tau sign under the same weight construction.

| Cohort | Model | Tau | Weight construction | Mean difference | 95% CI | Valid / requested | Ordering changed vs existing tau |
|---|---|---|---|---:|---|---:|---|
| mixed | coxph | existing_none | other_four_oof | -0.019353720696004207 | `[-0.10509036106493296, 0.0791708482505417]` | 9999 / 9999 | False |
| mixed | coxph | existing_none | outer_training_outcomes | -0.019353720696004207 | `[-0.10509036106493296, 0.0791708482505417]` | 9999 / 9999 | False |
| mixed | coxph | 3_year | other_four_oof | 0.008741796960769044 | `[-0.08014464719511778, 0.10670251350061215]` | 9972 / 9999 | True |
| mixed | coxph | 3_year | outer_training_outcomes | 0.008741796960769044 | `[-0.08014464719511778, 0.10670251350061215]` | 9972 / 9999 | True |
| mixed | coxph | 5_year | other_four_oof | 0.06563693517110736 | `[-0.020146286451139857, 0.15395506893465352]` | 9998 / 9999 | True |
| mixed | coxph | 5_year | outer_training_outcomes | 0.06563693517110736 | `[-0.020146286451139857, 0.15395506893465352]` | 9998 / 9999 | True |
| mixed | coxph | 7_year | other_four_oof | 0.04061283482310864 | `[-0.036221445121145536, 0.1229057158675549]` | 9999 / 9999 | True |
| mixed | coxph | 7_year | outer_training_outcomes | 0.04061283482310864 | `[-0.036221445121145536, 0.1229057158675549]` | 9999 / 9999 | True |
| mixed | rsf | existing_none | other_four_oof | -0.026936912174609473 | `[-0.12083535227712465, 0.08069741329545473]` | 9999 / 9999 | False |
| mixed | rsf | existing_none | outer_training_outcomes | -0.026936912174609473 | `[-0.12083535227712465, 0.08069741329545473]` | 9999 / 9999 | False |
| mixed | rsf | 3_year | other_four_oof | 0.04759695884328745 | `[-0.06378237886993106, 0.1799927419235375]` | 9972 / 9999 | True |
| mixed | rsf | 3_year | outer_training_outcomes | 0.04759695884328745 | `[-0.06378237886993106, 0.1799927419235375]` | 9972 / 9999 | True |
| mixed | rsf | 5_year | other_four_oof | 0.044171620722510775 | `[-0.04721743397636417, 0.13978343621502895]` | 9998 / 9999 | True |
| mixed | rsf | 5_year | outer_training_outcomes | 0.044171620722510775 | `[-0.04721743397636417, 0.13978343621502895]` | 9998 / 9999 | True |
| mixed | rsf | 7_year | other_four_oof | -0.021467217629513757 | `[-0.09751692674385629, 0.062608988215795]` | 9999 / 9999 | False |
| mixed | rsf | 7_year | outer_training_outcomes | -0.021467217629513757 | `[-0.09751692674385629, 0.062608988215795]` | 9999 / 9999 | False |
| mixed | deepsurv | existing_none | other_four_oof | -0.04978505621603002 | `[-0.1427022153784877, 0.07875228941058408]` | 9999 / 9999 | False |
| mixed | deepsurv | existing_none | outer_training_outcomes | -0.04978505621603002 | `[-0.1427022153784877, 0.07875228941058408]` | 9999 / 9999 | False |
| mixed | deepsurv | 3_year | other_four_oof | 0.09698403969765998 | `[-0.02219020377964897, 0.23863689545199274]` | 9972 / 9999 | True |
| mixed | deepsurv | 3_year | outer_training_outcomes | 0.09698403969765998 | `[-0.02219020377964897, 0.23863689545199274]` | 9972 / 9999 | True |
| mixed | deepsurv | 5_year | other_four_oof | 0.11787429071024744 | `[0.010503247831023691, 0.2258819115691073]` | 9998 / 9999 | True |
| mixed | deepsurv | 5_year | outer_training_outcomes | 0.11787429071024744 | `[0.010503247831023691, 0.2258819115691073]` | 9998 / 9999 | True |
| mixed | deepsurv | 7_year | other_four_oof | 0.06374413032554486 | `[-0.02703193370747662, 0.1568785384278921]` | 9999 / 9999 | True |
| mixed | deepsurv | 7_year | outer_training_outcomes | 0.06374413032554486 | `[-0.02703193370747662, 0.1568785384278921]` | 9999 / 9999 | True |
| primary | coxph | existing_none | other_four_oof | 0.016605331610658126 | `[-0.06592665131201647, 0.11701822670538543]` | 9999 / 9999 | False |
| primary | coxph | existing_none | outer_training_outcomes | 0.016605331610658126 | `[-0.06592665131201647, 0.11701822670538543]` | 9999 / 9999 | False |
| primary | coxph | 3_year | other_four_oof | 0.10661252477381548 | `[0.011698764265270434, 0.19547502439571712]` | 9839 / 9999 | False |
| primary | coxph | 3_year | outer_training_outcomes | 0.10661252477381548 | `[0.011698764265270434, 0.19547502439571712]` | 9839 / 9999 | False |
| primary | coxph | 5_year | other_four_oof | 0.10498879920702286 | `[0.01308690835494376, 0.1929483348617457]` | 9997 / 9999 | False |
| primary | coxph | 5_year | outer_training_outcomes | 0.10498879920702286 | `[0.01308690835494376, 0.1929483348617457]` | 9997 / 9999 | False |
| primary | coxph | 7_year | other_four_oof | 0.0913834692222349 | `[0.004200550730348545, 0.17966957719642002]` | 9999 / 9999 | False |
| primary | coxph | 7_year | outer_training_outcomes | 0.0913834692222349 | `[0.004200550730348545, 0.17966957719642002]` | 9999 / 9999 | False |
| primary | rsf | existing_none | other_four_oof | 0.01924583849949375 | `[-0.06743312306075193, 0.127483906523009]` | 9999 / 9999 | False |
| primary | rsf | existing_none | outer_training_outcomes | 0.01924583849949375 | `[-0.06743312306075193, 0.127483906523009]` | 9999 / 9999 | False |
| primary | rsf | 3_year | other_four_oof | 0.11013657801320427 | `[-0.0031884624808052, 0.21412844018882207]` | 9839 / 9999 | False |
| primary | rsf | 3_year | outer_training_outcomes | 0.11013657801320427 | `[-0.0031884624808052, 0.21412844018882207]` | 9839 / 9999 | False |
| primary | rsf | 5_year | other_four_oof | 0.11470267849088375 | `[0.017453026506140053, 0.20688458163527607]` | 9997 / 9999 | False |
| primary | rsf | 5_year | outer_training_outcomes | 0.11470267849088375 | `[0.017453026506140053, 0.20688458163527607]` | 9997 / 9999 | False |
| primary | rsf | 7_year | other_four_oof | 0.07902463886636446 | `[-0.002827072096000236, 0.16428935127754404]` | 9999 / 9999 | False |
| primary | rsf | 7_year | outer_training_outcomes | 0.07902463886636446 | `[-0.002827072096000236, 0.16428935127754404]` | 9999 / 9999 | False |
| primary | deepsurv | existing_none | other_four_oof | 0.07519112895366087 | `[0.000420027336828892, 0.18665923484060246]` | 9999 / 9999 | False |
| primary | deepsurv | existing_none | outer_training_outcomes | 0.07519112895366087 | `[0.000420027336828892, 0.18665923484060246]` | 9999 / 9999 | False |
| primary | deepsurv | 3_year | other_four_oof | 0.2502031567485644 | `[0.11406675470169332, 0.3799172360990463]` | 9839 / 9999 | False |
| primary | deepsurv | 3_year | outer_training_outcomes | 0.2502031567485644 | `[0.11406675470169332, 0.3799172360990463]` | 9839 / 9999 | False |
| primary | deepsurv | 5_year | other_four_oof | 0.12953916098895152 | `[0.017480171538106014, 0.25496504360720923]` | 9997 / 9999 | False |
| primary | deepsurv | 5_year | outer_training_outcomes | 0.12953916098895152 | `[0.017480171538106014, 0.25496504360720923]` | 9997 / 9999 | False |
| primary | deepsurv | 7_year | other_four_oof | 0.11898706088532318 | `[0.03344335482553859, 0.20938108794278035]` | 9999 / 9999 | False |
| primary | deepsurv | 7_year | outer_training_outcomes | 0.11898706088532318 | `[0.03344335482553859, 0.20938108794278035]` | 9999 / 9999 | False |

## Construction audit

- mixed: outer-training outcomes and other-four-OOF outcomes were identical in all folds = `True`.
- primary: outer-training outcomes and other-four-OOF outcomes were identical in all folds = `True`.

Fold-level IPCW diagnostic errors: `[]`.

All raw bootstrap draws, fold-level values, row/participant counts, and weight diagnostics are retained without rounding in the JSON.

Report finalization: `{"command": "C:\\moenv\\Scripts\\python.exe src/evaluation/phase17_task7_finalize.py --repo-dir . --json results/phase17/task7_uno_harrell_sensitivity.json --markdown results/phase17/task7_uno_harrell_sensitivity.md", "statistics_recomputed": false, "derived_fields": ["at-risk participant counts from stored OOF outcomes", "ordering-change booleans from stored mean differences", "non-finite bootstrap draws serialized as null"], "reason": "Make participant counting explicit and JSON standards-compliant; no C-index, confidence interval, or weight statistic was recomputed."}`.
