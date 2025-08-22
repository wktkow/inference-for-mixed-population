# AMT: Homework Assignment

## Inference for mixed populations, 2024-2025
---

### 1 Description of the data

The data considered in this assignment are from a study in which patients with renal deficiency, receiving hemodialysis are followed longitudinally and measured on a monthly basis. Every month, it is checked whether there is indication of iron deficiency defined as serum ferritin below 100µg and/or transferrin saturation below 20%.

The objective of the study is to investigate the distribution of iron deficiency, and how it relates to the age and sex of the patient.

### 2 Data file

*   SAS file Hemodialysismix.csv/rds
*   Variables:
    1.  **ID**: identification number of the subject
    2.  **AGE**: age of the patient at the time of entering the study (years)
    3.  **SEX**: sex of the patient (1: male; 2: female)
    4.  **NR**: the number of times serum ferritin and transferrin were measured
    5.  **NRIRON**: the number of times adequate iron stores were observed

### 3 Assignment

1.  The response of interest is the number of occasions where the iron stores were adequate
2.  Study the distribution of the response. Is there any indication of over- or underdispersion ?
3.  Fit a mixture distribution to the data, and explore the different components in the mixture
4.  Are the components related to the age and/or sex of the patient ? Do those covariates completely explain presence of potential clusters in the outcome ?

---

### 4 Points of Attention

0. The report should be written in markdown inside of a README.md file.
1.  The page limit for the submitted report is 5 pages (not including the appendix, if present).
3.  Your report should specify models fitted, interpretation of the parameters in the models, estimation methodology, software used, and motivations for all choices made.
4.  Furthermore, the results should be presented in a format that can be communicated to a non-statistical audience.
5.  Only report the most relevant tables and figures. Plain software output should not be part of the report. Any syntax considered important should be presented in the appendix.
