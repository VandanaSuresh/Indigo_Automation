**Overview**
This repository provides a hybrid automation pipeline for high-throughput analysis of Sanger sequencing chromatogram files (.ab1) generated from CRISPR/Cas9-edited plant lines.

    Gear Genomics INDIGO Webserver (web-based deconvolution & alignment)

    ICE (Inference of CRISPR Edits) – local fallback analysis engine

    Selenium WebDriver automation

    Biopython-based validation

    Automated CSV reporting

    Advanced logging and structured error handling
This Version 2 hybrid system ensures maximum reliability, accuracy, and scalability by automatically routing failed INDIGO analyses to ICE as a fallback.


**Required Python Packages**

pip install selenium pandas biopython

ICE must be locally cloned and referenced via:

ice_source_path = r"path_to_ice-master"
Citation

****If you use this hybrid automation system, please cite:**

Suresh, V., Girish, C., & Tavva, V.S.S.
Python-based automation of INDIGO webserver using Selenium for high-throughput analysis of CRISPR-induced allelic variations.
bioRxiv (2025).**

**Contact**
Dr. V. S. Sresty Tavva
Tata Institute for Genetics and Society (TIGS)
sresty.tavva@tigs.res.in
