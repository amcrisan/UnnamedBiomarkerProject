## A generalizable and improved framework for microbial biomarker discovery.

**Background :** LEfSe is a popular method for biomarker discovery, however its statistical framing has some issues that limits its generalizability. To begin, LEfSe relies on binary class comparisons and cannot be applied to continuous or truly multinomial response types. Furthermore, LEfSe does not include additional quantitative or qualitative sample metadata variable for biomarker discovery largely because of the limitations imposed by the linear discriminate analysis step, which does not handle such data well. The result is that biomarkers are not properly adjusted against other informative, and perhaps more easily obtained, metatdata variables. We propose an alternative method that does generalize to multiple repsonse types and can include metadata. We explore its impact in a previously published study (Schubert, 2014) 

**Methods:**  We compared biomarkers found by LEfSe in a previously published study against our new approach. Discovery and validation efforts used a total of 338 samples: 94 with C. difficile; 89 with Diarrheal but no C. difficle; and 155 healthy controls. Biomarkers were identified by three steps: 1) an abundance distribution filter; 2) a Kruskal-Wallis test; 3) a penalized logistic regression. We explored multinomial, binary and continuous response types, as well as found biomarkers with and without adjust against sample metadata variables. 

 **Results:** 
- Overlapped all but 2 Schubert OTUs when using all biomarkers found. The lack of overlap increased depending on strictness of further filter.
- Metadata variables dropped some OTUs, and added some new ones as well.
- *the important thing that remains to be done is answer the question "So what?"*
- **To Do:** Use OTU 19 as continuous response, see what else correlates to it.
- **To Do:** Binary response [need to dbl check some of the binom code).
- **To Do:** Report AUCs for biomarker models against base model + combined model.

**Conclusions:**  The limitations of the current LEfSE implementation may effect the type and utility of the biomarkers discovered and ultimately impact the adoption to these biomarkers for broader application in the clinic or elsewhere.  A penalized regression approach provides a more generalizable and robust framework for microbial biomarker discovery.
